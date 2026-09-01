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
#include <random>
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
    const SystemParams& system_params, bool highs_verbose,
    ApproximationMode mode)
    : t_max_(t_max),
      t_split_count_(t_split_count),
      tau_min_(tau_min),
      tau_max_(tau_max),
      system_params_(system_params),
      highs_verbose_(highs_verbose),
      approximation_mode_(mode) {
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
  const char* mode_str =
      approximation_mode_ == ApproximationMode::Upper ? "upper" : "lower";
  std::ostringstream out;
  out << "{\n"
      << "  \"t_max\": " << t_max_ << ",\n"
      << "  \"approximation_mode\": \"" << mode_str << "\",\n"
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
          system_params_.f5max, system_params_.f8max,
          // The 8D area vertex lists are the truncated ones -- 12 of an area's
          // 108-192 vertices -- and this path no longer uses them: it carries
          // the three block vertex sets instead. Materialising them correctly
          // would cost about 13 GB; the block data is under 1 MB.
          hcpwa::TriangleAreasOptions{.build_8d_regions = false});

  auto ingest = [](int phase, int block_id, const hcpwa::BlockRegions& src) {
    GlobalBlockGeometry block;
    block.coords = src.coords;
    block.coord_count = src.coord_count;
    block.num_regions = static_cast<int>(src.vertices.size());
    if (block.coord_count < 2 || block.coord_count > 3) {
      throw std::runtime_error(
          "getIntersectionPoints: block has an unexpected coord_count");
    }
    for (int d = 0; d < block.coord_count; ++d) {
      if (block.coords[d] < 0 || block.coords[d] >= kSpaceDim) {
        throw std::runtime_error(
            "getIntersectionPoints: block coordinate out of range");
      }
      if (d > 0 && block.coords[d] <= block.coords[d - 1]) {
        throw std::runtime_error(
            "getIntersectionPoints: block coordinates are not ascending");
      }
    }
    block.vertices.reserve(block.num_regions);
    for (int j = 0; j < block.num_regions; ++j) {
      // A stored block polytope is full dimensional, so it has at least three
      // vertices. Block C is a plane cell, which may be a general polygon on
      // this path rather than a triangle.
      if (src.vertices[j].size() < 3) {
        throw std::runtime_error(
            "getIntersectionPoints: block region has fewer than 3 vertices");
      }
      std::vector<Eigen::VectorXd> vertices;
      vertices.reserve(src.vertices[j].size());
      for (const auto& vertex : src.vertices[j]) {
        Eigen::VectorXd point(block.coord_count);
        for (int d = 0; d < block.coord_count; ++d) {
          point(d) = static_cast<double>(vertex[d]);
        }
        vertices.push_back(std::move(point));
      }
      block.vertices.push_back(std::move(vertices));
    }
    return block;
  };

  for (int b = 0; b < kBlockCount; ++b) {
    blocks_[0][b] = ingest(0, b, areas_vertices.blocks_phase0[b]);
    blocks_[1][b] = ingest(1, b, areas_vertices.blocks_phase1[b]);
  }

  for (int phase = 0; phase < 2; ++phase) {
    // The three blocks must partition the eight state coordinates; otherwise
    // the LP rows would not separate and the reduction would be invalid.
    std::array<int, kSpaceDim> coordinate_uses{};
    long long product = 1;
    for (int b = 0; b < kBlockCount; ++b) {
      const auto& block = blocks_[phase][b];
      product *= block.num_regions;
      for (int d = 0; d < block.coord_count; ++d) {
        ++coordinate_uses[block.coords[d]];
      }
    }
    for (int r = 0; r < kSpaceDim; ++r) {
      if (coordinate_uses[r] != 1) {
        throw std::runtime_error(std::format(
            "getIntersectionPoints: phase {} coordinate {} is claimed by {} "
            "blocks, expected exactly 1",
            phase, r, coordinate_uses[r]));
      }
    }
    const std::size_t area_count
        = phase == 0 ? areas_vertices.intersection_prism_indices_phase0.size()
                     : areas_vertices.intersection_prism_indices_phase1.size();
    if (product != static_cast<long long>(area_count)) {
      throw std::runtime_error(std::format(
          "getIntersectionPoints: phase {} has {} areas but M_A*M_B*M_C = {}",
          phase, area_count, product));
    }
    num_regions_[phase] = static_cast<int>(area_count);

    std::size_t block_rows = 0;
    for (int b = 0; b < kBlockCount; ++b) {
      for (const auto& vertices : blocks_[phase][b].vertices) {
        block_rows += vertices.size();
      }
    }
    logger_->info(
        "getIntersectionPoints: phase={}, M_A={}, M_B={}, M_C={}, areas={}, "
        "block rows R_A+R_B+R_C={}",
        phase, blocks_[phase][0].num_regions, blocks_[phase][1].num_regions,
        blocks_[phase][2].num_regions, num_regions_[phase], block_rows);
  }
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

Eigen::VectorXd GlobalAffineApproximator::blockCentroidCoords(
    int phase, int block_id, int j_block) const {
  if (phase < 0 || phase >= 2 || block_id < 0 || block_id >= kBlockCount) {
    throw std::invalid_argument("blockCentroidCoords: invalid phase or block");
  }
  const auto& block = blocks_[phase][block_id];
  if (j_block < 0 || j_block >= block.num_regions) {
    throw std::invalid_argument("blockCentroidCoords: invalid block region id");
  }
  const auto& vertices = block.vertices[j_block];
  if (vertices.empty()) {
    throw std::runtime_error("blockCentroidCoords: block region has no "
                             "vertices");
  }
  Eigen::VectorXd centroid = Eigen::VectorXd::Zero(block.coord_count);
  for (const auto& vertex : vertices) {
    centroid += vertex;
  }
  centroid /= static_cast<double>(vertices.size());
  return centroid;
}

// Expands a block centroid into a full 8-vector for getFIJMinResolution().
//
// Coordinates outside the block are NaN on purpose. Each block's flows read
// only that block's own coordinates -- phase 0 block A resolves f31 and f36
// over {0,2,5}, block B resolves f24 and f27 over {1,3,6} -- so the padding is
// never read. If it ever were, every comparison in getFIJMinResolution() would
// evaluate false and the function would throw rather than silently resolve a
// branch against a made-up value.
Eigen::VectorXd GlobalAffineApproximator::blockPointToFullState(
    int phase, int block_id, const Eigen::VectorXd& nu) const {
  const auto& block = blocks_[phase][block_id];
  Eigen::VectorXd full = Eigen::VectorXd::Constant(
      kSpaceDim, std::numeric_limits<double>::quiet_NaN());
  for (int d = 0; d < block.coord_count; ++d) {
    full(block.coords[d]) = nu(d);
  }
  return full;
}

std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd, double>
GlobalAffineApproximator::getBlockAMatrFVecGVecAndGScalJ(int phase,
                                                         int block_id,
                                                         int j_block) const {
  // The same CTM branch resolution as before, restricted to one block. A^(j) is
  // block diagonal with respect to the coordinate partition, so the restriction
  // loses nothing; restrict_row() asserts that by refusing any entry outside
  // the block. Block C is purely exogenous in both phases, so its A, f and g
  // are zero.
  const auto& block = blocks_[phase][block_id];
  const int c = block.coord_count;
  const Eigen::VectorXd nu = blockCentroidCoords(phase, block_id, j_block);
  const Eigen::VectorXd n = blockPointToFullState(phase, block_id, nu);

  Eigen::MatrixXd a_matr = Eigen::MatrixXd::Zero(c, c);
  Eigen::VectorXd b_vec = Eigen::VectorXd::Zero(c);
  Eigen::VectorXd g_vec = Eigen::VectorXd::Zero(c);
  double g_scal = 0.0;
  if (block_id == 2) {
    return std::make_tuple(a_matr, b_vec, g_vec, g_scal);
  }

  auto restrict_row = [&](const Eigen::RowVectorXd& row) {
    Eigen::RowVectorXd out = Eigen::RowVectorXd::Zero(c);
    std::array<bool, kSpaceDim> in_block{};
    for (int d = 0; d < c; ++d) {
      in_block[block.coords[d]] = true;
      out(d) = row(block.coords[d]);
    }
    for (int r = 0; r < kSpaceDim; ++r) {
      if (!in_block[r] && row(r) != 0.0) {
        throw std::runtime_error(std::format(
            "getBlockAMatrFVecGVecAndGScalJ: phase {} block {} region {} has a "
            "flow coefficient {} on coordinate {}, outside the block. A^(j) is "
            "supposed to be block diagonal.",
            phase, block_id, j_block, row(r), r));
      }
    }
    return out;
  };

  // Phase 0: block A holds n1, n3, n6 with flows f31, f36; block B holds
  // n2, n4, n7 with f24, f27. Phase 1: block A' holds n1, n5, n7 with f51,
  // f57; block B' holds n4, n6, n8 with f84, f86. Both flows of a block share
  // its source cell.
  int first_i = 0, first_j = 0, second_i = 0, second_j = 0;
  if (phase == 0 && block_id == 0) {
    first_i = 3 - 1;  first_j = 1 - 1;
    second_i = 3 - 1; second_j = 6 - 1;
  } else if (phase == 0 && block_id == 1) {
    first_i = 2 - 1;  first_j = 4 - 1;
    second_i = 2 - 1; second_j = 7 - 1;
  } else if (phase == 1 && block_id == 0) {
    first_i = 5 - 1;  first_j = 1 - 1;
    second_i = 5 - 1; second_j = 7 - 1;
  } else if (phase == 1 && block_id == 1) {
    first_i = 8 - 1;  first_j = 4 - 1;
    second_i = 8 - 1; second_j = 6 - 1;
  } else {
    throw std::invalid_argument(
        "getBlockAMatrFVecGVecAndGScalJ: invalid phase/block combination");
  }

  const auto [first_row, first_scal] = getFIJMinResolution(first_i, first_j, n);
  const auto [second_row, second_scal]
      = getFIJMinResolution(second_i, second_j, n);
  const Eigen::RowVectorXd first = restrict_row(first_row);
  const Eigen::RowVectorXd second = restrict_row(second_row);

  auto local = [&](int state_axis) {
    for (int d = 0; d < c; ++d) {
      if (block.coords[d] == state_axis) {
        return d;
      }
    }
    throw std::runtime_error("getBlockAMatrFVecGVecAndGScalJ: coordinate is "
                             "not in this block");
  };

  const int source = local(first_i);
  const int first_dest = local(first_j);
  const int second_dest = local(second_j);
  a_matr.row(first_dest) += first;
  b_vec(first_dest) += first_scal(0);
  a_matr.row(second_dest) += second;
  b_vec(second_dest) += second_scal(0);
  a_matr.row(source) -= first + second;
  b_vec(source) -= first_scal(0) + second_scal(0);

  g_vec = (first + second).transpose();
  g_scal = first_scal(0) + second_scal(0);

  // assertScalar is the only non-finite guard on this data, and it matters more
  // now that the representative point is NaN-padded outside the block: a NaN
  // reaching g_scal would otherwise flow through kappaFixed into every row of
  // this block and surface as an unrelated HiGHS status much later.
  hcpwa::util::assertShape(a_matr, c, c);
  hcpwa::util::assertShape(b_vec, c);
  hcpwa::util::assertShape(g_vec, c);
  hcpwa::util::assertScalar(g_scal);
  return std::make_tuple(a_matr, b_vec, g_vec, g_scal);
}

std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd, Eigen::VectorXd>
GlobalAffineApproximator::getBlockQQ(int phase, int block_id,
                                     int j_block) const {
  // The disturbance box is componentwise -- each component's bounds depend on
  // its own coordinate only -- so restricting it to a block is exact. The loop
  // runs over the block's own coordinates rather than all eight, which is what
  // makes the NaN padding safe.
  const auto& block = blocks_[phase][block_id];
  const int c = block.coord_count;
  const Eigen::VectorXd n0 = blockCentroidCoords(phase, block_id, j_block);

  Eigen::MatrixXd Q_upper = Eigen::MatrixXd::Zero(c, c);
  Eigen::VectorXd q_upper = Eigen::VectorXd::Zero(c);
  Eigen::MatrixXd Q_lower = Eigen::MatrixXd::Zero(c, c);
  Eigen::VectorXd q_lower = Eigen::VectorXd::Zero(c);

  const double N = system_params_.N;
  const double w = system_params_.w;
  const double v = system_params_.v;
  const double F = system_params_.F;

  for (int d = 0; d < c; ++d) {
    const int axis = block.coords[d];
    const double ni = n0(d);
    const bool is_in
        = std::find(kInIds.begin(), kInIds.end(), axis) != kInIds.end();
    const bool is_out
        = std::find(kOutIds.begin(), kOutIds.end(), axis) != kOutIds.end();
    if (is_in == is_out) {
      throw std::runtime_error(std::format(
          "getBlockQQ: coordinate {} is {} an inflow and an outflow cell", axis,
          is_in ? "both" : "neither"));
    }
    if (is_in) {
      auto [f_min, f_max] = getFMinMaxForAxis(axis);
      if (f_min < w * (N - ni)) {
        q_lower(d) = f_min;
      } else {
        Q_lower(d, d) = -w;
        q_lower(d) = w * N;
      }
      if (f_max < w * (N - ni)) {
        q_upper(d) = f_max;
      } else {
        Q_upper(d, d) = -w;
        q_upper(d) = w * N;
      }
    } else {
      if (F < v * ni) {
        q_lower(d) = -F;
      } else {
        Q_lower(d, d) = -v;
      }
    }
  }

  return std::make_tuple((Q_upper + Q_lower) / 2.0, (q_upper + q_lower) / 2.0,
                         (Q_upper - Q_lower) / 2.0, (q_upper - q_lower) / 2.0);
}

void GlobalAffineApproximator::precomputeSystemMatrices(int phase) {
  // About 509 small objects per phase, against 1.36M dense 8x8 blocks before.
  for (int b = 0; b < kBlockCount; ++b) {
    auto& block = blocks_[phase][b];
    block.a_matr.assign(block.num_regions, {});
    block.f_vec.assign(block.num_regions, {});
    block.q_c_matr.assign(block.num_regions, {});
    block.q_c_vec.assign(block.num_regions, {});
    block.q_r_matr.assign(block.num_regions, {});
    block.q_r_vec.assign(block.num_regions, {});
    block.g_vec.assign(block.num_regions, {});
    block.g_scal.assign(block.num_regions, 0.0);

    for (int j = 0; j < block.num_regions; ++j) {
      auto [a_matr, f_vec, g_vec, g_scal]
          = getBlockAMatrFVecGVecAndGScalJ(phase, b, j);
      auto [q_c_matr, q_c_vec, q_r_matr, q_r_vec] = getBlockQQ(phase, b, j);
      block.a_matr[j] = std::move(a_matr);
      block.f_vec[j] = std::move(f_vec);
      block.g_vec[j] = std::move(g_vec);
      block.g_scal[j] = g_scal;
      block.q_c_matr[j] = std::move(q_c_matr);
      block.q_c_vec[j] = std::move(q_c_vec);
      block.q_r_matr[j] = std::move(q_r_matr);
      block.q_r_vec[j] = std::move(q_r_vec);
    }
  }
}

hcpwa::util::global_block_lp::GlobalReducedLp
GlobalAffineApproximator::prepareLpMatrices(int phase) {
  // Builds the block-reduced LP. Three independent block row sets plus one
  // scalar coupling row, in place of one row per element of the Cartesian
  // product R_A x R_B x R_C.
  //
  // The reduction is exact: the feasible sets have identical projections onto
  // the decision variables (notes/step7_block_reduction.md, Lemmas 3 and 4).
  // Unlike the barycentric path there is no per-region y to identify, because
  // the modulus lift s is already a single global vector, so no analogue of
  // Lemma 5 is needed here.
  //
  // Row count on the production geometry: about 3.5e3 + 17, against 1.6e7 for
  // the truncated product form or roughly 2e8 for the product form built
  // correctly. The column count is unchanged apart from the three mu columns.
  namespace gblp = hcpwa::util::global_block_lp;

  const bool is_upper = (approximation_mode_ == ApproximationMode::Upper);
  const double sigma = is_upper ? 1.0 : -1.0;

  std::array<gblp::GlobalBlockInput, kBlockCount> inputs;
  for (int b = 0; b < kBlockCount; ++b) {
    const auto& block = blocks_[phase][b];
    gblp::GlobalBlockInput& input = inputs[b];
    input.coords = block.coords;
    input.coord_count = block.coord_count;

    // Normalised weights. Summing over the Cartesian product would weight each
    // block by the row count of the other two, which is an artefact of how the
    // sum is taken rather than a design decision, and the product objective was
    // never computable at this scale anyway. The resulting optimum is not the
    // one the product objective would give; validity of the bound does not
    // depend on the objective, only tightness does.
    input.objective_weight = 1.0;

    std::size_t total_rows = 0;
    for (const auto& vertices : block.vertices) {
      total_rows += vertices.size();
    }
    input.rows.reserve(total_rows);

    for (int j = 0; j < block.num_regions; ++j) {
      for (const Eigen::VectorXd& nu : block.vertices[j]) {
        gblp::GlobalBlockRow row;
        row.nu = nu;
        // p, r and kappa restricted to this block. Because A is block diagonal
        // and Q_c, Q_r are diagonal, these are exactly the block's slice of the
        // full-dimensional quantities: no approximation.
        row.p = hcpwa::util::coeffP(block.a_matr[j], block.f_vec[j],
                                    block.q_c_matr[j], block.q_c_vec[j], nu,
                                    t_delta_);
        row.r = hcpwa::util::radiusR(block.q_r_matr[j], block.q_r_vec[j], nu);
        row.kappa = hcpwa::util::kappaFixed(block.g_vec[j], block.g_scal[j], nu,
                                            t_delta_);
        input.rows.push_back(std::move(row));
      }
    }
  }

  gblp::GlobalReducedLpOptions options;
  options.dt = t_delta_;
  options.sigma = sigma;
  options.s_pin_weight = kSPinWeight;
  options.coefficient_eps = kEps;

  gblp::GlobalReducedLp lp = assembleGlobalReducedLp(inputs, options);

  logger_->info(
      "prepareLpMatrices: phase={}, cols={}, block regions={}/{}/{}, block "
      "rows={}/{}/{}, rows={}, nnz={}",
      phase, lp.num_cols, blocks_[phase][0].num_regions,
      blocks_[phase][1].num_regions, blocks_[phase][2].num_regions,
      lp.num_block_rows[0], lp.num_block_rows[1], lp.num_block_rows[2],
      lp.row_upper.size(), lp.values.size());

  return lp;
}

void GlobalAffineApproximator::validateStepResiduals(
    int phase, const std::vector<double>& v_prev,
    const std::vector<double>& v_cur) const {
  // Checks the ORIGINAL condition on the original index set:
  //   sigma * F(a,b,c) <= 0  for (a,b,c) in R_A x R_B x R_C,
  // with F additively separable across the three blocks. The reduced rows are
  // what the LP already enforced; the triples are what the bound property is
  // actually about, and evaluating one costs three lookups and two additions.
  //
  // Under the old truncation only 12 of an area's 108-192 vertices were
  // constrained, so almost any sampled triple would show sigma*F > 0. That is
  // exactly what this is here to catch.
  namespace gblp = hcpwa::util::global_block_lp;
  if (v_prev.size() != kVDeltaDim || v_cur.size() != kVDeltaDim) {
    throw std::invalid_argument("validateStepResiduals: bad value size");
  }
  const bool is_upper = (approximation_mode_ == ApproximationMode::Upper);
  const double sigma = is_upper ? 1.0 : -1.0;

  // Per block, the residual contribution of each of its rows. The constant
  // terms (v of the current step, and the previous step's constant) belong to
  // exactly one block, matching the row assembly.
  std::array<std::vector<double>, kBlockCount> contributions;
  for (int b = 0; b < kBlockCount; ++b) {
    const auto& block = blocks_[phase][b];
    const bool carries_const = (b == gblp::kVCarryingBlock);
    for (int j = 0; j < block.num_regions; ++j) {
      for (const Eigen::VectorXd& nu : block.vertices[j]) {
        const Eigen::VectorXd p
            = hcpwa::util::coeffP(block.a_matr[j], block.f_vec[j],
                                  block.q_c_matr[j], block.q_c_vec[j], nu,
                                  t_delta_);
        const Eigen::VectorXd r
            = hcpwa::util::radiusR(block.q_r_matr[j], block.q_r_vec[j], nu);
        const double kappa = hcpwa::util::kappaFixed(
            block.g_vec[j], block.g_scal[j], nu, t_delta_);

        // The row is
        //   sigma * (p^T V - v) + dt * r^T s + sigma * (kappa + b_upd) <= 0
        // with s = |V| at the optimum and b_upd = V_prev^T n + v_prev. Every
        // term carries sigma except the r term, whose sign is the same for both
        // bound directions.
        double value = sigma * kappa;
        for (int d = 0; d < block.coord_count; ++d) {
          const int axis = block.coords[d];
          value += sigma * p(d) * v_cur[axis];
          value += t_delta_ * r(d) * std::abs(v_cur[axis]);
          value += sigma * v_prev[axis] * nu(d);
        }
        if (carries_const) {
          value += -sigma * v_cur[kSpaceDim] + sigma * v_prev[kSpaceDim];
        }
        contributions[b].push_back(value);
      }
    }
    if (contributions[b].empty()) {
      throw std::runtime_error("validateStepResiduals: block has no rows");
    }
  }

  // EXACT, not sampled. The residual is additively separable across the three
  // blocks, so max over all R_A*R_B*R_C triples is the sum of the three
  // per-block maxima -- O(R_A+R_B+R_C) instead of a sample of ~2e8, and it
  // cannot miss a violation confined to one block row. contributions[] already
  // carry sigma, so maximising them directly is correct for both directions.
  const std::size_t total
      = contributions[0].size() * contributions[1].size()
        * contributions[2].size();
  double worst = 0.0;
  std::array<std::size_t, kBlockCount> worst_ids{};
  for (int b = 0; b < kBlockCount; ++b) {
    double best = -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < contributions[b].size(); ++i) {
      if (contributions[b][i] > best) {
        best = contributions[b][i];
        worst_ids[b] = i;
      }
    }
    worst += best;
  }

  if (worst > kResidualValidationTol) {
    throw std::runtime_error(std::format(
        "validateStepResiduals: phase {} is not a bound. worst residual = {}, "
        "witness triple (block rows {}, {}, {}) out of {} triples checked "
        "exhaustively.",
        phase, worst, worst_ids[0], worst_ids[1], worst_ids[2], total));
  }
  logger_->info(
      "validateStepResiduals: phase={}, {} triples covered exhaustively, "
      "worst={}",
      phase, total, worst);
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

double GlobalAffineApproximator::getMinBorderFuncValuesAtN(
    int theta_idx, const std::vector<int>& theta_end_ids, int max_switches,
    int phase, const Eigen::VectorXd& n) {
  double min_val = std::numeric_limits<double>::infinity();
  for (int theta_end_idx : theta_end_ids) {
    for (int r = 0; r < max_switches; ++r) {
      double val
          = getBorderFuncValuesAtN(r, theta_idx, theta_end_idx, phase, n);
      min_val = std::min(min_val, val);
    }
  }

  if (!std::isfinite(min_val)) {
    throw std::runtime_error(
        std::format("getMinBorderFuncValuesAtN: no finite border value for "
                    "theta_idx: {}, max_switches: {}, phase: {}",
                    theta_idx, max_switches, phase));
  }

  return min_val;
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

  const bool is_upper = (approximation_mode_ == ApproximationMode::Upper);

  Eigen::VectorXd n_bar = Eigen::VectorXd::Zero(kSpaceDim);
  for (const auto& vertex : cube_angle_vertices_) {
    n_bar += vertex;
  }
  const int n_corner = static_cast<int>(cube_angle_vertices_.size());
  if (n_corner == 0) {
    throw std::runtime_error("getBorderConditions: cube_angle_vertices_ is empty");
  }

  std::vector<double> c_vec(kVDeltaDim, 0.0);
  for (int i = 0; i < kSpaceDim; ++i) {
    c_vec[i] = is_upper ? n_bar(i) : -n_bar(i);
  }
  c_vec[kSpaceDim] = is_upper ? static_cast<double>(n_corner)
                              : -static_cast<double>(n_corner);

  std::vector<std::vector<double>> a_matr_lst;
  std::vector<double> b_vec_lower_lst;
  std::vector<double> b_vec_upper_lst;
  std::vector<double> f_scals;
  f_scals.reserve(n_corner);
  const double inf = highs.getInfinity();

  for (const auto& vertex : this->cube_angle_vertices_) {
    double f_scal = is_upper
                        ? getMaxBorderFuncValuesAtN(theta_idx, theta_range_ids,
                                                    switch_cnt, switch_phase,
                                                    vertex)
                        : getMinBorderFuncValuesAtN(theta_idx, theta_range_ids,
                                                    switch_cnt, switch_phase,
                                                    vertex);
    if (!std::isfinite(f_scal)) {
      logger_->error(
          "getBorderConditions: non-finite f_scal for switch_phase={}, "
          "theta_idx={}, theta={}, switch_cnt={}, vertex_norm={}",
          switch_phase, theta_idx, theta, switch_cnt, vertex.norm());
      throw std::runtime_error("getBorderConditions: non-finite f_scal");
    }
    f_scals.push_back(f_scal);

    std::vector<double> row(kVDeltaDim);
    if (is_upper) {
      // [-n^T, -1] x <= -f_max  <=>  V^T n + v >= f_max
      for (int i = 0; i < kSpaceDim; ++i) {
        row[i] = -vertex(i);
      }
      row[kSpaceDim] = -1.0;
      b_vec_lower_lst.push_back(-inf);
      b_vec_upper_lst.push_back(-f_scal);
    } else {
      // [n^T, 1] x <= f_min  <=>  V^T n + v <= f_min
      for (int i = 0; i < kSpaceDim; ++i) {
        row[i] = vertex(i);
      }
      row[kSpaceDim] = 1.0;
      b_vec_lower_lst.push_back(-inf);
      b_vec_upper_lst.push_back(f_scal);
    }
    a_matr_lst.push_back(std::move(row));
  }

  const int m = static_cast<int>(a_matr_lst.size());
  const int n = kVDeltaDim;

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
  highs.setOptionValue("simplex_strategy", 2);
  highs.setOptionValue("pdlp_optimality_tolerance", kHighsPdlpOptimalityTol);
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

  const auto& solution = highs.getSolution();
  if (solution.col_value.size() < static_cast<std::size_t>(kVDeltaDim)) {
    throw std::runtime_error("getBorderConditions: solution is too short");
  }
  std::vector<double> result(solution.col_value.begin(),
                             solution.col_value.begin() + kVDeltaDim);

  double min_slack = std::numeric_limits<double>::infinity();
  double max_slack = -std::numeric_limits<double>::infinity();
  double sum_slack = 0.0;
  for (std::size_t corner = 0; corner < cube_angle_vertices_.size(); ++corner) {
    const auto& vertex = cube_angle_vertices_[corner];
    double val = 0.0;
    for (int i = 0; i < kSpaceDim; ++i) {
      val += result[i] * vertex(i);
    }
    val += result[kSpaceDim];
    const double f_scal = f_scals[corner];
    const double slack = is_upper ? (val - f_scal) : (f_scal - val);
    if (slack < -10.0 * kHighsSolutionTol) {
      throw std::runtime_error(
          is_upper
              ? "getBorderConditions: border LP violates majorization constraint"
              : "getBorderConditions: border LP violates minorization constraint");
    }
    min_slack = std::min(min_slack, slack);
    max_slack = std::max(max_slack, slack);
    sum_slack += slack;
  }
  logger_->info(
      "Border LP switch_phase={} theta_idx={} switch_cnt={} mode={} "
      "slack_min={:.4f} slack_max={:.4f} slack_sum={:.4f}",
      switch_phase, theta_idx, switch_cnt, is_upper ? "upper" : "lower",
      min_slack, max_slack, sum_slack);

  return result;
}

std::tuple<std::unique_ptr<Highs>, std::vector<double>, std::vector<double>>
GlobalAffineApproximator::initializeHighs(int phase) {
  // Reads the per-phase assembly prepared by precomputeMatrices(). It must NOT
  // assemble here: run() builds the solvers from a thread pool, so assigning
  // reduced_lps_[phase] from several threads would race, and each thread holds
  // a reference into that object while handing its buffers to HiGHS.
  const hcpwa::util::global_block_lp::GlobalReducedLp& lp = reduced_lps_[phase];
  if (lp.row_upper.empty()) {
    throw std::runtime_error(
        "initializeHighs: LP not assembled. Call precomputeMatrices() first.");
  }
  std::vector<double> row_upper = lp.row_upper;
  const int m = static_cast<int>(row_upper.size());
  const int n = lp.num_cols;

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

  std::vector<double> col_lower = lp.col_lower;
  std::vector<double> col_upper = lp.col_upper;
  std::vector<double> row_lower = lp.row_lower;

  // Add all columns at once, with no matrix coefficients yet (we'll add rows
  // next). Signature: addCols(num_new_col, cost, lower, upper, num_nz, start,
  // index, value)
  HighsStatus st
      = highs->addCols(n, lp.objective.data(), col_lower.data(),
                       col_upper.data(),
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
  st = highs->addRows(m, row_lower.data(), row_upper.data(), lp.values.size(),
                      lp.starts.data(), lp.cols.data(), lp.values.data());
  if (st != HighsStatus::kOk) {
    throw std::runtime_error("InitializeHighs: highs.addRows failed.");
  }

  // Warm-start at the origin. Layout is now [V(8), v, s(8), mu(3)], and the
  // coupling row is sum_b mu_b <= 0, so the old all-ones start would hand the
  // solver a point that violates a constraint of the model by 3.
  HighsSolution initial_solution;
  initial_solution.value_valid = true;
  initial_solution.col_value.assign(n, 0.0);
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
  // Only the residual row right-hand sides depend on the previous step's value
  // function; the s rows, the coupling row and the column bounds are static.
  //
  // b_upd is that value function evaluated at a vertex, V_prev^T n + v_prev. It
  // splits the same way the rows do: the linear part restricted to the block,
  // and the constant only on the block that carries v. Cost drops from O(R) to
  // O(R_A + R_B + R_C).
  if (v_prev_vec.size() != kVDeltaDim) {
    throw std::invalid_argument("v_prev_vec size must be "
                                + std::to_string(kVDeltaDim));
  }
  const auto& row_lower = row_lowers_[solver_index];
  auto& highs_solver = highs_solvers_[solver_index];

  std::vector<double> new_row_upper
      = hcpwa::util::global_block_lp::globalReducedLpRowUpper(
          reduced_lps_[phase], v_prev_vec);
  for (double& upper : new_row_upper) {
    if (std::abs(upper) <= kEps) {
      upper = 0.0;
    }
  }

  std::vector<int> row_ids(new_row_upper.size());
  std::iota(row_ids.begin(), row_ids.end(), 0);
  auto st = highs_solver->changeRowsBounds(
      static_cast<int>(row_ids.size()), row_ids.data(), row_lower.data(),
      new_row_upper.data());
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
  for (int phase = 0; phase < 2; ++phase) {
    if (blocks_[phase][0].num_regions == 0) {
      throw std::runtime_error(
          "precomputeMatrices: geometry not computed. Call "
          "getIntersectionPoints() before precomputeMatrices().");
    }
  }

  for (int phase = 0; phase < kPhases; ++phase) {
    precomputeSystemMatrices(phase);
    std::size_t block_regions = 0;
    std::size_t block_rows = 0;
    for (int b = 0; b < kBlockCount; ++b) {
      block_regions += blocks_[phase][b].num_regions;
      for (const auto& vertices : blocks_[phase][b].vertices) {
        block_rows += vertices.size();
      }
    }
    logger_->info(
        "precomputeMatrices: phase={}, block regions={}, block rows={}, "
        "implied full areas={}",
        phase, block_regions, block_rows, num_regions_[phase]);

    // Assemble once per phase, here, while still single-threaded. The
    // constraint matrix is fixed for the whole backward march, so every solver
    // instance of this phase shares this one assembly.
    reduced_lps_[phase] = prepareLpMatrices(phase);
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
  solver_mutexes_.clear();
  const int solvers_per_phase = n_threads / 2;
  const int total_solvers = kPhases * solvers_per_phase;
  highs_solvers_.resize(total_solvers);
  row_lowers_.resize(total_solvers);
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
              if (validate_) {
                // Checks that the solved step really is a bound, on the
                // original product row set rather than the reduced rows.
                validateStepResiduals(phase, v_prev, v_next);
              }
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