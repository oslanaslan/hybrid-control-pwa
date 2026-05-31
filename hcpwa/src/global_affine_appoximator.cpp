#include "global_affine_approximator.hpp"

#include "util/affine_lp_utils.hpp"
#include "util/assert_utils.hpp"
#include "util/sparse_matrix_utils.hpp"
#include "utility.hpp"
#include <iostream>

#include <algo.hpp>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <format>
#include <fstream>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>
#include "cddwrap/lineareq.hpp"
#include "morph.hpp"
#include "types.hpp"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Highs.h>
#include <spdlog/spdlog.h>
#include "spdlog/sinks/stdout_color_sinks.h"
#include "util/HighsCDouble.h"
#include "thread_pool.hpp"

#include <future>

namespace {

// HiGHS stopping tolerances (tighter than defaults / previous 1e-5).
constexpr double kHighsSolutionTol = 1e-6;
constexpr double kHighsSmallMatrixValue = 1e-9;
constexpr double kHighsPdlpOptimalityTol = 1e-6;

// Collects t_index[i] for all i where t_range[i] lies in the interval
// [low - half_step, high + half_step] (if high_inclusive) or
// [low - half_step, high - half_step) (if !high_inclusive).
// Used to replace NumPy-style boolean indexing over the t grid.
std::vector<int> getTRangeIdsInInterval(const std::vector<double>& t_range,
                                        const std::vector<int>& t_index,
                                        double half_step, double low,
                                        double high, bool high_inclusive) {
  std::vector<int> result;
  for (size_t i = 0; i < t_range.size(); ++i) {
    if (t_range[i] < low - half_step) {
      continue;
    }
    if (high_inclusive) {
      if (t_range[i] <= high + half_step) {
        result.push_back(t_index[i]);
      }
    } else {
      if (t_range[i] < high - half_step) {
        result.push_back(t_index[i]);
      }
    }
  }
  return result;
}

}  // namespace

namespace global_affine_approximator {

const std::vector<int> kInIds = {2 - 1, 3 - 1, 5 - 1, 8 - 1};
const std::vector<int> kOutIds = {1 - 1, 4 - 1, 6 - 1, 7 - 1};
const int kPhases = 2;

GlobalAffineApproximator::GlobalAffineApproximator(
    double t_max, int t_split_count, double tau_min, double tau_max,
    const SystemParams& system_params, bool highs_verbose)
    : t_max_(t_max),
      t_split_count_(t_split_count),
      tau_min_(tau_min),
      tau_max_(tau_max),
      system_params_(system_params),
      highs_verbose_(highs_verbose) {
  if (t_split_count < 1) {
    throw std::invalid_argument(
        "t_split_count must be at least 1 in GlobalAffineApproximator "
        "constructor.");
  }
  if (tau_min >= tau_max || tau_min < 0.0 || tau_max < 0.0) {
    throw std::invalid_argument(
        "tau_min must be less than tau_max and both must be >= 0 in "
        "GlobalAffineApproximator constructor.");
  }
  max_switches_ = static_cast<int>(std::ceil(t_max_ / tau_min_));
  t_range_.resize(t_split_count);
  t_index_.resize(t_split_count);
  t_delta_ = t_split_count > 1 ? std::abs(t_max / (t_split_count - 1)) : 0.0;
  for (int i = 0; i < t_split_count; ++i) {
    t_range_[i] = (t_split_count == 1) ? t_max : i * t_delta_;
    t_index_[i] = i;
  }
  int n_vertices = 1 << kSpaceDim;
  cube_angle_vertices_.clear();
  for (int vert = 0; vert < n_vertices; ++vert) {
    Eigen::VectorXd v(kSpaceDim);
    for (int d = 0; d < kSpaceDim; ++d) {
      v(d) = ((vert & (1 << d)) != 0) ? this->system_params_.N : 0.0;
    }
    cube_angle_vertices_.push_back(v);
  }
  theta_t_index_lists_ = interval_building::buildThetaToTIndexLists(
      t_max_, tau_min_, tau_max_, t_range_, max_switches_);
  logger_ = spdlog::stdout_color_mt("affine_approximator");
  logger_->set_level(spdlog::level::info);

  if (highs_verbose_) {
    interval_building::prettyPrintThetaTLists(theta_t_index_lists_, t_range_);
  }
}

void GlobalAffineApproximator::dumpInitParamsToJson(
    const std::string& filepath) const {
  const SystemParams& p = system_params_;
  std::ostringstream out;
  out << "{\n"
      << "  \"t_max\": " << t_max_ << ",\n"
      << "  \"t_split_count\": " << t_split_count_ << ",\n"
      << "  \"max_switches\": " << max_switches_ << ",\n"
      << "  \"tau_min\": " << tau_min_ << ",\n"
      << "  \"tau_max\": " << tau_max_ << ",\n"
      << "  \"system_params\": {\n"
      << "    \"N\": " << p.N << ",\n"
      << "    \"F\": " << p.F << ",\n"
      << "    \"v\": " << p.v << ",\n"
      << "    \"w\": " << p.w << ",\n"
      << "    \"b51\": " << p.b51 << ",\n"
      << "    \"b57\": " << p.b57 << ",\n"
      << "    \"b84\": " << p.b84 << ",\n"
      << "    \"b86\": " << p.b86 << ",\n"
      << "    \"b31\": " << p.b31 << ",\n"
      << "    \"b36\": " << p.b36 << ",\n"
      << "    \"b24\": " << p.b24 << ",\n"
      << "    \"b27\": " << p.b27 << ",\n"
      << "    \"f2min\": " << p.f2min << ",\n"
      << "    \"f3min\": " << p.f3min << ",\n"
      << "    \"f5min\": " << p.f5min << ",\n"
      << "    \"f8min\": " << p.f8min << ",\n"
      << "    \"f2max\": " << p.f2max << ",\n"
      << "    \"f3max\": " << p.f3max << ",\n"
      << "    \"f5max\": " << p.f5max << ",\n"
      << "    \"f8max\": " << p.f8max << "\n"
      << "  }\n"
      << "}\n";
  std::ofstream f(filepath);
  if (!f) {
    throw std::runtime_error("dumpInitParamsToJson: cannot open file "
                             + filepath + " for writing");
  }
  f << out.str();
}

double GlobalAffineApproximator::getBetaParamForAxis(int i, int j) const {
  if (i == 5 - 1 && j == 1 - 1) {
    return system_params_.b51;
  } else if (i == 5 - 1 && j == 7 - 1) {
    return system_params_.b57;
  } else if (i == 8 - 1 && j == 4 - 1) {
    return system_params_.b84;
  } else if (i == 8 - 1 && j == 6 - 1) {
    return system_params_.b86;
  } else if (i == 3 - 1 && j == 1 - 1) {
    return system_params_.b31;
  } else if (i == 3 - 1 && j == 6 - 1) {
    return system_params_.b36;
  } else if (i == 2 - 1 && j == 4 - 1) {
    return system_params_.b24;
  } else if (i == 2 - 1 && j == 7 - 1) {
    return system_params_.b27;
  } else {
    throw std::invalid_argument(std::string("No such beta params for axis (")
                                + std::to_string(i) + ", " + std::to_string(j)
                                + ")");
  }
}

std::pair<double, double> GlobalAffineApproximator::getFMinMaxForAxis(
    int i) const {
  if (i == 2 - 1) {
    return std::make_pair(system_params_.f2min, system_params_.f2max);
  } else if (i == 3 - 1) {
    return std::make_pair(system_params_.f3min, system_params_.f3max);
  } else if (i == 5 - 1) {
    return std::make_pair(system_params_.f5min, system_params_.f5max);
  } else if (i == 8 - 1) {
    return std::make_pair(system_params_.f8min, system_params_.f8max);
  } else {
    throw std::invalid_argument(std::string("No such f min max for axis ")
                                + std::to_string(i));
  }
}

void GlobalAffineApproximator::getIntersectionPoints() {
  hcpwa::PolygonAreasVerticesResult areas_vertices
      = hcpwa::compute_polygon_areas_vertices(
          system_params_.N, system_params_.F, system_params_.v,
          system_params_.w, system_params_.b51, system_params_.b57,
          system_params_.b84, system_params_.b86, system_params_.b31,
          system_params_.b36, system_params_.b24, system_params_.b27,
          system_params_.f2min, system_params_.f3min, system_params_.f5min,
          system_params_.f8min, system_params_.f2max, system_params_.f3max,
          system_params_.f5max, system_params_.f8max);

  intersection_points_.resize(2);

  auto convert_phase = [](const std::vector<std::vector<hcpwa::Vec<8>>>& src,
                          std::vector<std::vector<Eigen::VectorXd>>& dst) {
    dst.resize(src.size());
    for (int a = 0; a < src.size(); ++a) {
      const auto& area = src[a];
      dst[a].resize(area.size());
      for (int v = 0; v < area.size(); ++v) {
        Eigen::VectorXd vec(kSpaceDim);
        for (int i = 0; i < kSpaceDim; ++i) {
          vec(i) = static_cast<double>(area[v][i]);
        }
        dst[a][v] = std::move(vec);
      }
    }
  };

  convert_phase(areas_vertices.intersection_points_phase0,
                intersection_points_[0]);
  convert_phase(areas_vertices.intersection_points_phase1,
                intersection_points_[1]);
}

std::pair<Eigen::RowVectorXd, Eigen::RowVectorXd>
GlobalAffineApproximator::getFIJMinResolution(int i, int j,
                                              const Eigen::VectorXd& n) const {
  // f_i_j_min_a = min{beta_i_j * F, beta_i_j * v * n_i, w(N − n_j)}
  // Where f_matr_row is row vector [length SPACE_DIM], f_vec_row is 1x1 (just
  // value)
  const double N = system_params_.N;
  const double F = system_params_.F;
  const double v = system_params_.v;
  const double w = system_params_.w;
  double beta_i_j = getBetaParamForAxis(i, j);
  double n_i = n(i);
  double n_j = n(j);
  Eigen::RowVectorXd f_matr_row
      = Eigen::RowVectorXd::Zero(kSpaceDim);                   // 1 x SPACE_DIM
  Eigen::RowVectorXd f_vec_row = Eigen::RowVectorXd::Zero(1);  // 1 x 1

  double a = beta_i_j * v * n_i;
  double b = beta_i_j * F;
  double c = w * (N - n_j);

  if (b < a + kEps && b < c + kEps) {
    // Case: b is smallest
    f_vec_row(0) = beta_i_j * F;
  } else if (c < a + kEps && c < b + kEps) {
    // Case: c is smallest
    f_matr_row(j) = -w;
    f_vec_row(0) = w * N;
  } else if (a < b + kEps && a < c + kEps) {
    // Case: a is smallest
    f_matr_row(i) = beta_i_j * v;
    // f_vec_row remains zero
  } else {
    std::ostringstream oss;
    oss << "No such case: a=" << a << ", b=" << b << ", c=" << c
        << " (getFIJMinResolution)";
    throw std::invalid_argument(oss.str());
  }
  return std::make_pair(f_matr_row, f_vec_row);
}

Eigen::VectorXd GlobalAffineApproximator::areaCentroidCoords(int j,
                                                             int phase) const {
  // Check if intersection_points are computed for given phase
  if (intersection_points_.empty() || phase < 0
      || phase >= static_cast<int>(intersection_points_.size())
      || intersection_points_[phase].empty()) {
    throw std::runtime_error(
        "Intersection points are not computed for the specified phase. "
        "Please "
        "compute intersection points before calling areaCentroidCoords.");
  }
  const std::vector<std::vector<Eigen::VectorXd>>& intersection_points_per_phase
      = intersection_points_[phase];
  if (j < 0 || j >= static_cast<int>(intersection_points_per_phase.size())) {
    throw std::invalid_argument("Invalid area index j in areaCentroidCoords");
  }

  const auto& area = intersection_points_per_phase[j];
  if (area.empty()) {
    throw std::invalid_argument("Area is empty in areaCentroidCoords");
  }

  Eigen::VectorXd centroid = Eigen::VectorXd::Zero(kSpaceDim);
  for (const auto& point : area) {
    centroid += point;
  }
  centroid /= static_cast<double>(area.size());

  hcpwa::util::assertShape(centroid, kSpaceDim);
  return centroid;
}

std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd, double>
GlobalAffineApproximator::getAMatrFVecGVecAndGScalJ(int j, int phase) const {
  Eigen::VectorXd n = areaCentroidCoords(j, phase);

  Eigen::MatrixXd a_matr;
  Eigen::VectorXd b_vec;
  Eigen::VectorXd g_vec;
  double g_scal;

  if (phase == 0) {
    auto [f31_matr_row, f31_vec_row] = getFIJMinResolution(3 - 1, 1 - 1, n);
    auto [f36_matr_row, f36_vec_row] = getFIJMinResolution(3 - 1, 6 - 1, n);
    auto [f24_matr_row, f24_vec_row] = getFIJMinResolution(2 - 1, 4 - 1, n);
    auto [f27_matr_row, f27_vec_row] = getFIJMinResolution(2 - 1, 7 - 1, n);

    // Build A_matr: SPACE_DIM x SPACE_DIM
    a_matr.resize(kSpaceDim, kSpaceDim);
    a_matr.row(0) = f31_matr_row;
    a_matr.row(1) = -f24_matr_row - f27_matr_row;
    a_matr.row(2) = -f31_matr_row - f36_matr_row;
    a_matr.row(3) = f24_matr_row;
    a_matr.row(4) = Eigen::RowVectorXd::Zero(kSpaceDim);
    a_matr.row(5) = f36_matr_row;
    a_matr.row(6) = f27_matr_row;
    a_matr.row(7) = Eigen::RowVectorXd::Zero(kSpaceDim);

    // Build b_vec: SPACE_DIM x 1
    b_vec.resize(kSpaceDim);
    b_vec(0) = f31_vec_row(0);
    b_vec(1) = -f24_vec_row(0) - f27_vec_row(0);
    b_vec(2) = -f31_vec_row(0) - f36_vec_row(0);
    b_vec(3) = f24_vec_row(0);
    b_vec(4) = 0.0;
    b_vec(5) = f36_vec_row(0);
    b_vec(6) = f27_vec_row(0);
    b_vec(7) = 0.0;

    // Build g_vec: SPACE_DIM x 1
    g_vec = (f31_matr_row + f36_matr_row + f24_matr_row + f27_matr_row)
                .transpose();

    // Build g_scal: scalar
    g_scal = f31_vec_row(0) + f36_vec_row(0) + f24_vec_row(0) + f27_vec_row(0);
  } else if (phase == 1) {
    auto [f51_matr_row, f51_vec_row] = getFIJMinResolution(5 - 1, 1 - 1, n);
    auto [f57_matr_row, f57_vec_row] = getFIJMinResolution(5 - 1, 7 - 1, n);
    auto [f84_matr_row, f84_vec_row] = getFIJMinResolution(8 - 1, 4 - 1, n);
    auto [f86_matr_row, f86_vec_row] = getFIJMinResolution(8 - 1, 6 - 1, n);

    // Build A_matr: SPACE_DIM x SPACE_DIM
    a_matr.resize(kSpaceDim, kSpaceDim);
    a_matr.row(0) = f51_matr_row;
    a_matr.row(1) = Eigen::RowVectorXd::Zero(kSpaceDim);
    a_matr.row(2) = Eigen::RowVectorXd::Zero(kSpaceDim);
    a_matr.row(3) = f84_matr_row;
    a_matr.row(4) = -f51_matr_row - f57_matr_row;
    a_matr.row(5) = f86_matr_row;
    a_matr.row(6) = f57_matr_row;
    a_matr.row(7) = -f84_matr_row - f86_matr_row;

    // Build b_vec: SPACE_DIM x 1
    b_vec.resize(kSpaceDim);
    b_vec(0) = f51_vec_row(0);
    b_vec(1) = 0.0;
    b_vec(2) = 0.0;
    b_vec(3) = f84_vec_row(0);
    b_vec(4) = -f51_vec_row(0) - f57_vec_row(0);
    b_vec(5) = f86_vec_row(0);
    b_vec(6) = f57_vec_row(0);
    b_vec(7) = -f84_vec_row(0) - f86_vec_row(0);

    // Build g_vec: SPACE_DIM x 1
    g_vec = (f51_matr_row + f57_matr_row + f84_matr_row + f86_matr_row)
                .transpose();

    // Build g_scal: scalar
    g_scal = f51_vec_row(0) + f57_vec_row(0) + f84_vec_row(0) + f86_vec_row(0);
  } else {
    std::ostringstream oss;
    oss << "No such phase: " << phase;
    throw std::invalid_argument(oss.str());
  }

  hcpwa::util::assertShape(b_vec, kSpaceDim);
  hcpwa::util::assertShape(a_matr, kSpaceDim, kSpaceDim);
  hcpwa::util::assertShape(g_vec, kSpaceDim);
  hcpwa::util::assertScalar(g_scal);

  return std::make_tuple(a_matr, b_vec, g_vec, g_scal);
}

std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd, Eigen::VectorXd>
GlobalAffineApproximator::getQQForArea(int j, int phase) const {
  // Get representative point in Ω(j)
  Eigen::VectorXd n0 = areaCentroidCoords(j, phase);

  Eigen::MatrixXd Q_upper_matr = Eigen::MatrixXd::Zero(kSpaceDim, kSpaceDim);
  Eigen::VectorXd q_upper_vec = Eigen::VectorXd::Zero(kSpaceDim);
  Eigen::MatrixXd Q_lower_matr = Eigen::MatrixXd::Zero(kSpaceDim, kSpaceDim);
  Eigen::VectorXd q_lower_vec = Eigen::VectorXd::Zero(kSpaceDim);

  const double N = system_params_.N;
  const double w = system_params_.w;
  const double v = system_params_.v;
  const double F = system_params_.F;

  // Process IN_IDS
  for (int i : kInIds) {
    auto ni = n0(i);
    auto [f_min, f_max] = getFMinMaxForAxis(i);
    // lower
    if (f_min < w * (N - ni)) {
      q_lower_vec(i) = f_min;
    } else {
      Q_lower_matr(i, i) = -w;
      q_lower_vec(i) = w * N;
    }
    // upper
    if (f_max < w * (N - ni)) {
      q_upper_vec(i) = f_max;
    } else {
      Q_upper_matr(i, i) = -w;
      q_upper_vec(i) = w * N;
    }
  }

  // Process OUT_IDS
  for (int i : kOutIds) {
    auto ni = n0(i);
    // lower
    if (F < v * ni) {
      q_lower_vec(i) = -F;
    } else {
      Q_lower_matr(i, i) = -v;
    }
  }

  // Q_c = 1/2 * (Q_upper + Q_lower)
  // Q_r = 1/2 * (Q_upper - Q_lower)
  Eigen::MatrixXd Q_c = (Q_upper_matr + Q_lower_matr) / 2;
  Eigen::MatrixXd Q_r = (Q_upper_matr - Q_lower_matr) / 2;
  Eigen::VectorXd q_c = (q_upper_vec + q_lower_vec) / 2;
  Eigen::VectorXd q_r = (q_upper_vec - q_lower_vec) / 2;

  return std::make_tuple(Q_c, q_c, Q_r, q_r);
}

//
std::tuple<std::vector<Eigen::MatrixXd>,  // A_j_matrs_lst
           std::vector<Eigen::VectorXd>,  // f_j_vecs_lst
           std::vector<Eigen::MatrixXd>,  // Q_c_j_matrs_lst
           std::vector<Eigen::VectorXd>,  // q_c_j_vecs_lst
           std::vector<Eigen::MatrixXd>,  // Q_r_j_matrs_lst
           std::vector<Eigen::VectorXd>,  // q_r_j_vecs_lst
           std::vector<Eigen::VectorXd>,  // g_j_vecs_lst
           std::vector<double>            // g_j_scals_lst
           >
GlobalAffineApproximator::precomputeSystemMatrices(int phase) {
  auto intersection_points_phase = intersection_points_[phase];
  int n_areas = intersection_points_phase.size();
  std::vector<Eigen::MatrixXd> A_j_matrs_lst;
  std::vector<Eigen::VectorXd> f_j_vecs_lst;
  std::vector<Eigen::MatrixXd> Q_c_j_matrs_lst;
  std::vector<Eigen::VectorXd> q_c_j_vecs_lst;
  std::vector<Eigen::MatrixXd> Q_r_j_matrs_lst;
  std::vector<Eigen::VectorXd> q_r_j_vecs_lst;
  std::vector<Eigen::VectorXd> g_j_vecs_lst;
  std::vector<double> g_j_scals_lst;
  A_j_matrs_lst.reserve(n_areas);
  f_j_vecs_lst.reserve(n_areas);
  Q_c_j_matrs_lst.reserve(n_areas);
  q_c_j_vecs_lst.reserve(n_areas);
  Q_r_j_matrs_lst.reserve(n_areas);
  q_r_j_vecs_lst.reserve(n_areas);
  g_j_vecs_lst.reserve(n_areas);
  g_j_scals_lst.reserve(n_areas);

  for (int j = 0; j < n_areas; ++j) {
    auto [A_j_matr, f_j_vec, g_j_vec, g_j_scal]
        = getAMatrFVecGVecAndGScalJ(j, phase);
    auto [Q_c_j_matr, q_c_j_vec, Q_r_j_matr, q_r_j_vec]
        = getQQForArea(j, phase);
    hcpwa::util::assertShape(A_j_matr, kSpaceDim, kSpaceDim);
    hcpwa::util::assertShape(f_j_vec, kSpaceDim);
    hcpwa::util::assertShape(g_j_vec, kSpaceDim);
    hcpwa::util::assertScalar(g_j_scal);
    A_j_matrs_lst.emplace_back(A_j_matr);
    f_j_vecs_lst.emplace_back(f_j_vec);
    Q_c_j_matrs_lst.emplace_back(Q_c_j_matr);
    q_c_j_vecs_lst.emplace_back(q_c_j_vec);
    Q_r_j_matrs_lst.emplace_back(Q_r_j_matr);
    q_r_j_vecs_lst.emplace_back(q_r_j_vec);
    g_j_vecs_lst.emplace_back(g_j_vec);
    g_j_scals_lst.emplace_back(g_j_scal);
  }

  return std::make_tuple(A_j_matrs_lst, f_j_vecs_lst, Q_c_j_matrs_lst,
                         q_c_j_vecs_lst, Q_r_j_matrs_lst, q_r_j_vecs_lst,
                         g_j_vecs_lst, g_j_scals_lst);
}

// ---------- Main builder: builds the (j)-segment feasibility rows only
// ----------
//
// Decision vector is x = [ V ; v ; s ], size = (m) + 1 + (m) = 2m+1.
//
// For each vertex n in vertices_j, we add ONE feasibility constraint:
//
//    p(n)^T V - 1*v + dt*r(n)^T s <= -kappa(n)
//
// Right-hand side is split:
//   -kappa(n) = -kappa_fix(n) - kappa_t(n)
// where kappa_t(n) = V_prev^T n + v_prev is applied later via RHS updates.
//
void GlobalAffineApproximator::buildLPSegmentForJ(
    const Eigen::MatrixXd& Aj,   // A^{(j)} (m x m)
    const Eigen::VectorXd& fj,   // f^{(j)} (m)
    const Eigen::MatrixXd& Qc,   // Q_c^{(j)} (m x m)
    const Eigen::VectorXd& qc,   // q_c^{(j)} (m)
    const Eigen::MatrixXd& Qr,   // Q_r^{(j)} (m x m)
    const Eigen::VectorXd& qr,   // q_r^{(j)} (m)
    const Eigen::VectorXd& g_n,  // g_{i,n}^{(j)} (m)
    double g0,                   // g_{i,0}^{(j)} (scalar)
    const std::vector<Eigen::VectorXd>&
        vertices_j,  // vertices of Omega^{(j)}, each (m)
    double dt,
    // CSR format vectors for the A (constraint) matrix of the LP
    std::vector<int>& starts,     // Row start indices (CSR)
    std::vector<int>& cols,       // Column indices (CSR)
    std::vector<double>& values,  // Nonzero values (CSR)
    std::vector<double>&
        b_fixed,  // fixed RHS part of Ax <= b_fixed + b_time
    Eigen::VectorXd& obj_p_accum, int& vertex_count) {
  if (vertices_j.empty()) {
    throw std::runtime_error("build_lp_segment_for_j: vertices_j is empty.");
  }

  const int n_vertices = static_cast<int>(vertices_j.size());
  const int n_cols = kLPCols;  // [V(m), v(1), s(m)]

  // Basic dimension checks (lightweight)
  if (Aj.rows() != kSpaceDim || Aj.cols() != kSpaceDim) {
    throw std::runtime_error("Aj has wrong shape.");
  }
  if (fj.size() != kSpaceDim) {
    throw std::runtime_error("fj has wrong size.");
  }
  if (Qc.rows() != kSpaceDim || Qc.cols() != kSpaceDim) {
    throw std::runtime_error("Qc has wrong shape.");
  }
  if (qc.size() != kSpaceDim) {
    throw std::runtime_error("qc has wrong size.");
  }
  if (Qr.rows() != kSpaceDim || Qr.cols() != kSpaceDim) {
    throw std::runtime_error("Qr has wrong shape.");
  }
  if (qr.size() != kSpaceDim) {
    throw std::runtime_error("qr has wrong size.");
  }
  if (g_n.size() != kSpaceDim) {
    throw std::runtime_error("g_n has wrong size.");
  }

  for (int k = 0; k < n_vertices; ++k) {
    const Eigen::VectorXd& n = vertices_j[k];
    if (n.size() != kSpaceDim) {
      throw std::runtime_error("Vertex has inconsistent dimension.");
    }

    const Eigen::VectorXd p = hcpwa::util::coeffP(Aj, fj, Qc, qc, n, dt);
    const Eigen::VectorXd r = hcpwa::util::radiusR(Qr, qr, n);
    double kappa_fix = hcpwa::util::kappaFixed(g_n, g0, n, dt);

    Eigen::RowVectorXd feas_row = Eigen::RowVectorXd::Zero(n_cols);
    feas_row.block(0, 0, 1, kSpaceDim) = p.transpose();
    feas_row(kSpaceDim) = -1.0;
    feas_row.block(0, kSpaceDim + 1, 1, kSpaceDim) = dt * r.transpose();

    if (std::abs(kappa_fix) > kEps) {
      b_fixed.push_back(-kappa_fix);
    } else {
      b_fixed.push_back(0);
    }
    hcpwa::util::csrAppendRow(starts, cols, values, feas_row, kEps);

    obj_p_accum += p;
    ++vertex_count;
  }
}

std::tuple<std::vector<int>, std::vector<int>, std::vector<double>,
           std::vector<double>, Eigen::RowVectorXd>
GlobalAffineApproximator::prepareLpMatrices(int phase) {
  // Gather all per-area matrix/vector data for this phase from member
  // variables
  const std::vector<Eigen::MatrixXd>& a_j_matrs = this->A_j_matrs_[phase];
  const std::vector<Eigen::VectorXd>& f_j_vecs = this->f_j_vecs_[phase];
  const std::vector<Eigen::MatrixXd>& q_c_j_matrs = this->Q_c_j_matrs_[phase];
  const std::vector<Eigen::VectorXd>& q_c_j_vecs = this->q_c_j_vecs_[phase];
  const std::vector<Eigen::MatrixXd>& q_r_j_matrs = this->Q_r_j_matrs_[phase];
  const std::vector<Eigen::VectorXd>& q_r_j_vecs = this->q_r_j_vecs_[phase];
  const std::vector<Eigen::VectorXd>& g_j_vecs = this->g_j_vecs_[phase];
  const std::vector<double>& g_j_scals = this->g_j_scals_[phase];
  const std::vector<Eigen::MatrixXd>& b_j_matrs = this->b_j_matrs_[phase];
  const std::vector<Eigen::VectorXd>& b_j_vecs = this->b_j_vecs_[phase];
  const auto& intersection_points_phase = intersection_points_[phase];
  int n_areas = intersection_points_phase.size();
  std::vector<int> starts = {0};
  std::vector<int> col_index;
  std::vector<double> value;
  std::vector<double> row_upper;
  Eigen::VectorXd obj_p_accum = Eigen::VectorXd::Zero(kSpaceDim);
  int vertex_count = 0;
  // n_areas = 1000000;

  for (int j = 0; j < n_areas; ++j) {
    buildLPSegmentForJ(a_j_matrs[j], f_j_vecs[j], q_c_j_matrs[j], q_c_j_vecs[j],
                       q_r_j_matrs[j], q_r_j_vecs[j], g_j_vecs[j], g_j_scals[j],
                       intersection_points_phase[j], t_delta_, starts,
                       col_index, value, row_upper, obj_p_accum, vertex_count);
  }
  // s >= V and s >= -V
  // [V, v, s], where s \in R^kSpaceDim
  for (int i = 0; i < kSpaceDim; i++) {
    Eigen::RowVectorXd row = Eigen::RowVectorXd::Zero(kLPCols);
    row(i) = 1;
    row(i + kSpaceDim + 1) = -1;
    hcpwa::util::csrAppendRow(starts, col_index, value, row);
  }
  for (int i = 0; i < kSpaceDim; i++) {
    Eigen::RowVectorXd row = Eigen::RowVectorXd::Zero(kLPCols);
    row(i) = -1;
    row(i + kSpaceDim + 1) = -1;
    hcpwa::util::csrAppendRow(starts, col_index, value, row);
  }
  for (int i = 0; i < 2 * kSpaceDim; i++) {
    row_upper.push_back(0);
  }

  Eigen::RowVectorXd c_vec = Eigen::RowVectorXd::Zero(kLPCols);
  c_vec.head(kSpaceDim) = -obj_p_accum.transpose();
  c_vec(kSpaceDim) = static_cast<double>(vertex_count);
  // TODO(experiment): Temporarily disable s-regularization in the objective
  // to test whether it is unnecessary.
  // === s-pinning penalty (Option 1, see plan): pin s_i = |V_i| ===
  c_vec.segment(kSpaceDim + 1, kSpaceDim)
      = Eigen::RowVectorXd::Constant(kSpaceDim, kSPinWeight);
  // ==============================================================

  return std::make_tuple(std::move(starts), std::move(col_index),
                         std::move(value), std::move(row_upper),
                         std::move(c_vec));
}

double GlobalAffineApproximator::getBorderFuncValuesAtN(
    int r, int theta_idx, int theta_end_idx, int phase,
    const Eigen::VectorXd& n) {
  std::vector<double> x;
  try {
    x = value_function_.get(phase, r, theta_idx, theta_end_idx);
  } catch (const std::invalid_argument&) {
    // TDODO: think later about the correct value: 0 or -inf
    // return -std::numeric_limits<double>::infinity();
    return 0.0;
  }
  if (x.size() != kVDeltaDim) {
    logger_->error(
        "getBorderFuncValuesAtN: invalid value_function size {} (expected {}) "
        "for phase={}, r={}, theta_idx={}, theta_end_idx={}",
        x.size(), kVDeltaDim, phase, r, theta_idx, theta_end_idx);
    throw std::runtime_error(
        "getBorderFuncValuesAtN: invalid value_function vector size");
  }
  double result = 0.0;
  for (int i = 0; i < kSpaceDim; ++i) {
    result += x[i] * n(i);
  }
  result += x[kSpaceDim];
  if (!std::isfinite(result)) {
    logger_->error(
        "getBorderFuncValuesAtN: non-finite result for phase={}, r={}, "
        "theta_idx={}, theta_end_idx={}, n_norm={}",
        phase, r, theta_idx, theta_end_idx, n.norm());
    throw std::runtime_error("getBorderFuncValuesAtN: non-finite result");
  }
  return result;
}

double GlobalAffineApproximator::getMaxBorderFuncValuesAtN(
    int theta_idx, const std::vector<int>& theta_end_ids, int max_switches,
    int phase, const Eigen::VectorXd& n) {
  double max_val = -std::numeric_limits<double>::infinity();
  for (int theta_end_idx : theta_end_ids) {
    for (int r = 0; r < max_switches; ++r) {
      double val
          = getBorderFuncValuesAtN(r, theta_idx, theta_end_idx, phase, n);
      max_val = std::max(max_val, val);
    }
  }

  // Value function is \geq 0 when computed so if max_val < 0, it means that no
  // any value function is computed for this theta_idx, theta_end_idx, phase, n
  // and that is an error.
  if (max_val < -kEps) {
    throw std::runtime_error(
        std::format("getMaxBorderFuncValuesAtN: max_val < 0 for theta_idx: {}, "
                    "max_switches: {}, phase: {}, max_val: {}",
                    theta_idx, max_switches, phase, max_val));
  }

  return max_val;
}

std::vector<double> GlobalAffineApproximator::getBorderConditions(
    int switch_phase, int theta_idx, double theta, int switch_cnt) {
  // If switch_cnt == 0, return zeros
  if (switch_cnt == 0) {
    return std::vector<double>(kVDeltaDim, 0.0);
  }
  Highs highs;

  // Calculate theta_min and theta_max
  double theta_min = std::min(theta + this->tau_min_, this->t_max_);
  double theta_max = std::min(theta + this->tau_max_, this->t_max_);

  const double half_step = this->t_delta_ / 2.0;
  std::vector<int> theta_range_ids = getTRangeIdsInInterval(
      this->t_range_, this->t_index_, half_step, theta_min, theta_max, true);

  if (theta_range_ids.empty()) {
    throw std::runtime_error("getBorderConditions: No theta range found");
  }

  // Initialize c_vec as zeros
  std::vector<double> c_vec(kVDeltaDim + 1, 0.0);
  c_vec[kVDeltaDim] = 1.0;  // z

  // Build a_matr and b_vec
  std::vector<std::vector<double>> a_matr_lst;
  std::vector<double> b_vec_lower_lst;
  std::vector<double> b_vec_upper_lst;
  // auto inf = std::numeric_limits<double>::infinity();
  const double inf = highs.getInfinity();

  for (const auto& vertex : this->cube_angle_vertices_) {
    double f_scal = getMaxBorderFuncValuesAtN(theta_idx, theta_range_ids,
                                              switch_cnt, switch_phase, vertex);
    if (!std::isfinite(f_scal)) {
      logger_->error(
          "getBorderConditions: non-finite f_scal for switch_phase={}, "
          "theta_idx={}, theta={}, switch_cnt={}, vertex_norm={}",
          switch_phase, theta_idx, theta, switch_cnt, vertex.norm());
      throw std::runtime_error("getBorderConditions: non-finite f_scal");
    }

    // Build row: [vertex, 1, 0] for Ax >= b
    std::vector<double> row_lower(kVDeltaDim + 1);
    // Build row: [vertex, 1, 0] for z >= Ax -> Ax - z <= 0
    std::vector<double> row_upper(kVDeltaDim + 1);
    for (int i = 0; i < kSpaceDim; ++i) {
      row_lower[i] = vertex(i);  // V
      row_upper[i] = vertex(i);  // V
    }
    row_lower[kVDeltaDim - 1] = 1.0;  // v
    row_lower[kVDeltaDim] = 0.0;      // z
    row_upper[kVDeltaDim - 1] = 1.0;  // v
    row_upper[kVDeltaDim] = -1.0;     // z

    b_vec_lower_lst.push_back(f_scal);
    b_vec_upper_lst.push_back(inf);
    b_vec_lower_lst.push_back(-inf);
    b_vec_upper_lst.push_back(0);
    a_matr_lst.push_back(row_lower);
    a_matr_lst.push_back(row_upper);
  }

  // Convert a_matr_lst to a single matrix (for Highs sparse format)
  const int m = static_cast<int>(a_matr_lst.size());  // 2 * 2^kSpaceDim
  const int n = static_cast<int>(c_vec.size());       // kSpaceDim + 2

  // Build sparse matrix representation for Highs
  std::vector<int> starts(m + 1, 0);
  std::vector<int> col_index;
  std::vector<double> value;

  int nnz = 0;
  for (int i = 0; i < m; ++i) {
    starts[i] = nnz;
    for (int j = 0; j < n; ++j) {
      const double a_val = a_matr_lst[i][j];
      if (std::abs(a_val) > 1e-10) {  // Skip near-zero values
        col_index.push_back(j);
        value.push_back(a_val);
        ++nnz;
      }
    }
  }
  starts[m] = nnz;

  // Set up and solve LP with Highs
  highs.setOptionValue("solver", "simplex");
  highs.setOptionValue("presolve", "on");
  highs.setOptionValue("log_to_console", highs_verbose_);
  highs.changeObjectiveSense(ObjSense::kMinimize);
  highs.setOptionValue("kkt_tolerance", kHighsSolutionTol);
  highs.setOptionValue("primal_feasibility_tolerance", kHighsSolutionTol);
  highs.setOptionValue("dual_feasibility_tolerance", kHighsSolutionTol);
  highs.setOptionValue("primal_residual_tolerance", kHighsSolutionTol);
  highs.setOptionValue("dual_residual_tolerance", kHighsSolutionTol);
  highs.setOptionValue("optimality_tolerance", kHighsSolutionTol);
  highs.setOptionValue("small_matrix_value", kHighsSmallMatrixValue);

  // Add columns (variables)
  std::vector<double> col_lower(n, -inf);
  std::vector<double> col_upper(n, inf);

  HighsStatus st
      = highs.addCols(n, c_vec.data(), col_lower.data(), col_upper.data(), 0,
                      nullptr, nullptr, nullptr);
  if (st != HighsStatus::kOk) {
    throw std::runtime_error(
        "getBorderConditions: highs.addCols failed for "
        + std::to_string(switch_phase) + ", " + std::to_string(theta_idx) + ", "
        + std::to_string(theta) + ", " + std::to_string(switch_cnt));
  }

  // Add rows (constraints)
  st = highs.addRows(m, b_vec_lower_lst.data(), b_vec_upper_lst.data(), nnz,
                     starts.data(), col_index.data(), value.data());
  if (st != HighsStatus::kOk) {
    throw std::runtime_error(
        "getBorderConditions: highs.addRows failed for "
        + std::to_string(switch_phase) + ", " + std::to_string(theta_idx) + ", "
        + std::to_string(theta) + ", " + std::to_string(switch_cnt));
  }

  // Solve
  st = highs.run();
  if (st != HighsStatus::kOk) {
    throw std::runtime_error(
        "getBorderConditions: highs.run() failed for "
        + std::to_string(switch_phase) + ", " + std::to_string(theta_idx) + ", "
        + std::to_string(theta) + ", " + std::to_string(switch_cnt));
  }

  if (highs.getModelStatus() != HighsModelStatus::kOptimal) {
    throw std::runtime_error(
        "getBorderConditions: LP solution not found for "
        + std::to_string(switch_phase) + ", " + std::to_string(theta_idx) + ", "
        + std::to_string(theta) + ", " + std::to_string(switch_cnt)
        + ", model status: "
        + std::to_string(static_cast<int>(highs.getModelStatus())));
  }

  // Extract solution
  const auto& solution = highs.getSolution();
  std::vector<double> result(solution.col_value.begin(),
                             solution.col_value.end() - 1);  // exclude z

  return result;
}

std::tuple<std::unique_ptr<Highs>, std::vector<double>, std::vector<double>>
GlobalAffineApproximator::initializeHighs(int phase) {
  auto [starts, col_index, value, row_upper, c_vec] = prepareLpMatrices(phase);
  const int m = row_upper.size();
  const int n = c_vec.size();

  std::unique_ptr<Highs> highs = std::make_unique<Highs>();

  // highs->setOptionValue("solver", "pdlp");
  highs->setOptionValue("solver", "simplex");
  highs->setOptionValue("presolve", "on");
  // Simplex conf
  highs->setOptionValue("simplex_strategy", 2);
  // PDLP conf
  highs->setOptionValue("pdlp_optimality_tolerance", kHighsPdlpOptimalityTol);
  highs->setOptionValue("kkt_tolerance", kHighsSolutionTol);
  highs->setOptionValue("primal_feasibility_tolerance", kHighsSolutionTol);
  highs->setOptionValue("dual_feasibility_tolerance", kHighsSolutionTol);
  highs->setOptionValue("primal_residual_tolerance", kHighsSolutionTol);
  highs->setOptionValue("dual_residual_tolerance", kHighsSolutionTol);
  highs->setOptionValue("optimality_tolerance", kHighsSolutionTol);
  highs->setOptionValue("small_matrix_value", kHighsSmallMatrixValue);
  // Logging and etc
  highs->setOptionValue("log_to_console", highs_verbose_);

  // L1-residual inner LP objective is minimized.
  highs->changeObjectiveSense(ObjSense::kMinimize);

  const double inf = highs->getInfinity();
  std::vector<double> col_lower(n, -inf);
  std::vector<double> col_upper(n, inf);
  std::vector<double> row_lower(m, -inf);

  // Add all columns at once, with no matrix coefficients yet (we'll add rows
  // next). Signature: addCols(num_new_col, cost, lower, upper, num_nz, start,
  // index, value)
  HighsStatus st
      = highs->addCols(n, c_vec.data(), col_lower.data(), col_upper.data(),
                       /*num_nz=*/0,
                       /*start=*/nullptr,
                       /*index=*/nullptr,
                       /*value=*/nullptr);
  if (st != HighsStatus::kOk) {
    throw std::runtime_error("InitializeHighs: highs.addCols failed.");
  }

  // ---- Add rows with their coefficients ----
  // Signature: addRows(num_new_row, lower, upper, num_nz, start, index,
  // value)
  st = highs->addRows(m, row_lower.data(), row_upper.data(), value.size(),
                      starts.data(), col_index.data(), value.data());
  if (st != HighsStatus::kOk) {
    throw std::runtime_error("InitializeHighs: highs.addRows failed.");
  }

  // Warm-start the main LP with x0 = ones.
  // Layout: [V(kSpaceDim), v, s(kSpaceDim)].
  HighsSolution initial_solution;
  initial_solution.value_valid = true;
  initial_solution.col_value.assign(n, 1.0);
  st = highs->setSolution(initial_solution);
  if (st != HighsStatus::kOk) {
    throw std::runtime_error("InitializeHighs: highs.setSolution failed.");
  }

  return std::make_tuple<std::unique_ptr<Highs>, std::vector<double>,
                         std::vector<double>>(
      std::move(highs), std::move(row_lower), std::move(row_upper));
}

void GlobalAffineApproximator::updateHighsRhsUpperBounds(
    int phase, int solver_index, const std::vector<double>& v_prev_vec) {
  if (v_prev_vec.size() != kVDeltaDim) {
    throw std::invalid_argument("v_prev_vec size must be "
                                + std::to_string(kVDeltaDim));
  }
  const auto& row_lower = row_lowers_[solver_index];
  const auto& row_upper = row_uppers_[solver_index];
  auto& highs_solver = highs_solvers_[solver_index];
  std::vector<double> new_row_upper(row_lower.size());
  const auto& intersection_points_phase = intersection_points_[phase];
  std::vector<int> row_ids(row_upper.size());
  std::iota(row_ids.begin(), row_ids.end(), 0);

  int j = 0;
  for (auto& area : intersection_points_phase) {
    for (auto& vertex : area) {
      double b_upd = 0.0;
      for (int i = 0; i < kSpaceDim; ++i) {
        b_upd += v_prev_vec[i] * vertex(i);
      }
      b_upd += v_prev_vec[kSpaceDim];
      new_row_upper[j] = row_upper[j] - b_upd;
      if (std::abs(new_row_upper[j]) <= kEps) {
        new_row_upper[j] = 0.0;
      }
      j++;
    }
  }

  auto st = highs_solver->changeRowsBounds(
      row_upper.size(), row_ids.data(), row_lower.data(), new_row_upper.data());
  if (st != HighsStatus::kOk) {
    throw std::runtime_error(
        "updateHighsRhsUpperBounds: highs.changeRowsUpperBounds failed.");
  }
}

std::vector<double> GlobalAffineApproximator::solveLp(int solver_index) {
  auto& highs_solver = highs_solvers_[solver_index];

  // Run the solver
  HighsStatus run_status = highs_solver->run();
  if (run_status != HighsStatus::kOk) {
    // logger_->warn("solveLp: highs_solver.run() failed with status {}",
    //               static_cast<int>(run_status));
    throw std::runtime_error("solveLp: highs_solver.run() failed with status "
                             + std::to_string(static_cast<int>(run_status)));
  }

  // Check if solution is optimal
  if (highs_solver->getModelStatus() != HighsModelStatus::kOptimal) {
    // logger_->warn(
    // "solveLp: LP solution not found for solver_index {}, model status: {}",
    // solver_index, static_cast<int>(highs_solver->getModelStatus()));
    throw std::runtime_error(
        "solveLp: LP solution not found for solver_index "
        + std::to_string(solver_index) + ", model status: "
        + std::to_string(static_cast<int>(highs_solver->getModelStatus())));
  }

  // Extract solution
  const auto& solution = highs_solver->getSolution();
  std::vector<double> result(solution.col_value.begin(),
                             solution.col_value.end());

  // x_next = [V, v, s], extract only v_next = [V, v]
  std::vector<double> v_next(kVDeltaDim);
  std::copy(result.begin(), result.begin() + kVDeltaDim, v_next.begin());

  return v_next;
}

void GlobalAffineApproximator::precomputeMatrices() {
  logger_->info("Starting precomputeMatrices");
  if (intersection_points_.empty()) {
    throw std::runtime_error(
        "precomputeMatrices: intersection points not computed. Call "
        "getIntersectionPoints() before precomputeMatrices().");
  }
  if (!b_j_matrs_.empty() || !b_j_vecs_.empty() || !g_j_vecs_.empty()
      || !g_j_scals_.empty()) {
    throw std::runtime_error(
        "precomputeMatrices: b_j_matrs, b_j_vecs, g_j_vecs, or g_j_scals "
        "are "
        "not empty.");
  }

  for (int phase = 0; phase < kPhases; ++phase) {
    auto [A_j_matrs_lst, f_j_vecs_lst, Q_c_j_matrs_lst, q_c_j_vecs_lst,
          Q_r_j_matrs_lst, q_r_j_vecs_lst, g_j_vecs_lst, g_j_scals_lst]
        = precomputeSystemMatrices(phase);
    this->A_j_matrs_.push_back(std::move(A_j_matrs_lst));
    this->f_j_vecs_.push_back(std::move(f_j_vecs_lst));
    this->Q_c_j_matrs_.push_back(std::move(Q_c_j_matrs_lst));
    this->q_c_j_vecs_.push_back(std::move(q_c_j_vecs_lst));
    this->Q_r_j_matrs_.push_back(std::move(Q_r_j_matrs_lst));
    this->q_r_j_vecs_.push_back(std::move(q_r_j_vecs_lst));
    this->g_j_vecs_.push_back(std::move(g_j_vecs_lst));
    this->g_j_scals_.push_back(std::move(g_j_scals_lst));
  }
  logger_->info("Finished precomputeMatrices");
}

void GlobalAffineApproximator::run(const std::string& output_folder_path,
                                   int n_threads) {
  const std::filesystem::path out_path(output_folder_path);
  if (!std::filesystem::exists(out_path)
      || !std::filesystem::is_directory(out_path)) {
    throw std::runtime_error(
        "output_folder_path does not exist or is not a directory: "
        + output_folder_path);
  }
  if (n_threads < 2 || n_threads % 2 != 0) {
    throw std::runtime_error(
        "n_threads must be at least 2 and even (so each phase can be run in "
        "parallel). Provided n_threads: "
        + std::to_string(n_threads));
  }
  logger_->info("Starting affine approximator");
  getIntersectionPoints();
  precomputeMatrices();

  // Initialize n_threads HiGHS solvers (n_threads/2 per phase)
  logger_->info("Start initializing Highs solvers");
  highs_solvers_.clear();
  row_lowers_.clear();
  row_uppers_.clear();
  solver_mutexes_.clear();
  const int solvers_per_phase = n_threads / 2;
  const int total_solvers = kPhases * solvers_per_phase;
  highs_solvers_.resize(total_solvers);
  row_lowers_.resize(total_solvers);
  row_uppers_.resize(total_solvers);
  solver_mutexes_.reserve(n_threads);
  for (int i = 0; i < n_threads; ++i) {
    solver_mutexes_.push_back(std::make_unique<std::mutex>());
  }
  ThreadPool init_pool(n_threads);
  std::vector<std::future<void>> init_futures;
  init_futures.reserve(static_cast<size_t>(total_solvers));
  for (int phase = 0; phase < kPhases; ++phase) {
    for (int s = 0; s < solvers_per_phase; ++s) {
      init_futures.push_back(
          init_pool.enqueue([this, phase, s, solvers_per_phase]() {
            const int solver_index = phase * solvers_per_phase + s;
            auto [highs_solver, row_lower, row_upper] = initializeHighs(phase);
            highs_solvers_[solver_index] = std::move(highs_solver);
            row_lowers_[solver_index] = std::move(row_lower);
            row_uppers_[solver_index] = std::move(row_upper);
          }));
    }
  }
  for (auto& f : init_futures) {
    f.get();
  }
  logger_->info("Done initializing Highs solvers");

  ThreadPool pool(n_threads);
  const double half_step = this->t_delta_ / 2;

  // Update HiGHS solvers for each switch count
  for (int switch_cnt = 0; switch_cnt <= max_switches_; ++switch_cnt) {
    logger_->info("Computing value function for switch count {}/{}", switch_cnt,
                  max_switches_);

    auto theta_range_ids
        = theta_t_index_lists_.expanded_t_by_k_theta[switch_cnt];
    const auto n = static_cast<std::ptrdiff_t>(theta_range_ids.size());
    std::vector<std::future<void>> futures;
    futures.reserve(static_cast<size_t>(n * kPhases));

    for (auto [theta_idx, t_range_ids] : theta_range_ids) {
      for (int phase = 0; phase < kPhases; ++phase) {
        const int solver_index
            = phase * (n_threads / 2)
              + (static_cast<int>(theta_idx) % (n_threads / 2));
        futures.push_back(pool.enqueue([this, switch_cnt, phase, theta_idx,
                                        t_range_ids, solver_index]() {
          const double theta = t_range_[theta_idx];
          logger_->info(
              "Starting computation for theta_idx: {}, phase: {}, switch_cnt: "
              "{}",
              theta_idx, phase, switch_cnt);

          std::lock_guard<std::mutex> lock(*solver_mutexes_[solver_index]);

          // int switch_phase = phase == 0 ? 1 : 0;
          // std::vector<double> v_prev
          //     = getBorderConditions(switch_phase, theta_idx, theta,
          //     switch_cnt);
          // if (v_prev.size() != kVDeltaDim) {
          //   throw std::invalid_argument("v_prev size must be "
          //                               + std::to_string(kVDeltaDim));
          // }
          // value_function_.set(phase, switch_cnt, theta_idx, theta_idx,
          // v_prev,
          //                     kVDeltaDim);

          // auto [min_val, max_val, free_val] = MinMaxVAndAffineTerm(v_prev);
          // logger_->debug(
          //     "Border conditions min and max at theta_idx: {}, phase: {}, "
          //     "switch_cnt: {}: \t{:.4f}\t{:.4f}\t{:.4f}",
          //     theta_idx, phase, switch_cnt, min_val, max_val, free_val);

          // if (t_theta_range_ids.empty()) {
          //   throw std::runtime_error(
          //       "t_theta_range_ids is empty: No previous timesteps to compute
          //       " "value function for. " "t_theta_min = "
          //       + std::to_string(t_theta_min)
          //       + ", theta = " + std::to_string(theta));
          // }

          std::vector<double> v_next, v_prev;
          const auto n_t = static_cast<std::ptrdiff_t>(t_range_ids.size());
          if (n_t == 0) {
            throw std::runtime_error(
                "run: empty t_range_ids for theta_idx=" + std::to_string(theta_idx)
                + ", phase=" + std::to_string(phase)
                + ", switch_cnt=" + std::to_string(switch_cnt));
          }

          for (std::ptrdiff_t i_t_idx = n_t - 1; i_t_idx >= 0; --i_t_idx) {
            int t_idx = t_range_ids[i_t_idx];

            if (i_t_idx == n_t - 1) {
              // Last element of the range must be a switch time at which we
              // compute the value border conditions
              if (t_idx != theta_idx) {
                throw std::runtime_error("t_idx != theta_idx: "
                                         + std::to_string(t_idx)
                                         + " != " + std::to_string(theta_idx));
              }
              int switch_phase = phase == 0 ? 1 : 0;
              v_next = getBorderConditions(switch_phase, theta_idx, theta,
                                           switch_cnt);
            } else {
              updateHighsRhsUpperBounds(phase, solver_index, v_prev);
              v_next = solveLp(solver_index);
            }

            if (v_next.size() != kVDeltaDim) {
              throw std::invalid_argument("v_next size must be "
                                          + std::to_string(kVDeltaDim));
            }
            value_function_.set(phase, switch_cnt, t_idx, theta_idx, v_next,
                                kVDeltaDim);
            auto [min_val, max_val, free_val] =
                MinMaxVAndAffineTerm(v_next, kSpaceDim);
            logger_->info(
                "Value function min and max at t_idx: {}, theta_idx: {}, "
                "phase: {}, "
                "switch_cnt: {}: \t{:.4f}\t{:.4f}\t{:.4f}",
                t_idx, theta_idx, phase, switch_cnt, min_val, max_val,
                free_val);
            v_prev = std::move(v_next);
          }
        }));
      }
    }

    for (auto& f : futures) {
      f.get();
    }
  }

  // Save results
  const std::filesystem::path base(output_folder_path);
  value_function_.dumpToJson((base / "value_function.json").string());
  hcpwa::util::dumpVectorToJson(t_range_, (base / "t_range.json").string());
  dumpInitParamsToJson((base / "init_params.json").string());
}
}  // namespace global_affine_approximator