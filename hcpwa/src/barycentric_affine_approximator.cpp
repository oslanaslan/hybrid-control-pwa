#include "barycentric_affine_approximator.hpp"

#include "thread_pool.hpp"
#include "util/assert_utils.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Highs.h>
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <format>
#include <fstream>
#include <future>
#include <limits>
#include <numeric>
#include <span>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

// Keep local helper names and public methods aligned with the existing
// GlobalAffineApproximator naming style so the copied solver structure remains
// easy to audit side-by-side.
// NOLINTBEGIN(readability-identifier-naming)

namespace {

// HiGHS tolerances are copied from the production global affine solver. The new
// barycentric LP is meant to reuse the same solver configuration and update only
// RHS bounds between backward time steps.
constexpr double kHighsSolutionTol = 1e-6;
constexpr double kHighsSmallMatrixValue = 1e-9;
constexpr double kHighsPdlpOptimalityTol = 1e-6;

// Converts hcpwa::Vec<2> to Eigen::Vector2d. Keeping this tiny conversion helper
// avoids mixing two vector APIs inside the indexing-heavy barycentric code.
Eigen::Vector2d toEigen2(const hcpwa::Vec<2>& v) {
  return Eigen::Vector2d(static_cast<double>(v[0]), static_cast<double>(v[1]));
}

// Converts one hcpwa::Vec<8> region vertex into the Eigen representation used by
// the LP formulas. Every LP formula below assumes n is an 8-vector.
Eigen::VectorXd toEigen8(const hcpwa::Vec<8>& v) {
  Eigen::VectorXd out(barycentric_affine_approximator::kSpaceDim);
  for (int i = 0; i < barycentric_affine_approximator::kSpaceDim; ++i) {
    out(i) = static_cast<double>(v[i]);
  }
  return out;
}

// The triangle/prism indices returned by compute_triangle_areas_vertices() are
// ordered by the way compute_intersection_points() combines prisms. These axis
// arrays must stay in that same order because region_triangle_ids[j][s] selects
// a triangle in layers[s].
// Moved to barycentric_geometry_types.hpp so the courier solver and the tests
// share one source of truth for this table.
using barycentric_affine_approximator::projectionAxesForPhase;

// Appends a sparse row to the CSR arrays. This mirrors csrAppendRow() but avoids
// building dense rows with eta_i + 8M columns, which can be large for the
// barycentric basis.
void appendSparseRow(std::vector<int>& starts, std::vector<int>& cols,
                     std::vector<double>& values,
                     const barycentric_affine_approximator::SparseVec& row,
                     double eps = barycentric_affine_approximator::kEps) {
  int nnz = static_cast<int>(values.size());
  for (std::size_t i = 0; i < row.cols.size(); ++i) {
    if (std::abs(row.vals[i]) <= eps) {
      continue;
    }
    cols.push_back(row.cols[i]);
    values.push_back(row.vals[i]);
    ++nnz;
  }
  starts.push_back(nnz);
}

// Computes beta = Psi_j^T * w as a sparse vector over the x block. This is the
// implementation of beta_{j,nu} = Psi_j^T (A_j nu + f_j + c_{j,nu}).
barycentric_affine_approximator::SparseVec psiTransposeTimes(
    const barycentric_affine_approximator::SparsePsi& psi,
    const Eigen::VectorXd& w) {
  barycentric_affine_approximator::SparseVec result;
  for (int r = 0; r < barycentric_affine_approximator::kSpaceDim; ++r) {
    for (std::size_t k = 0; k < psi.rows[r].cols.size(); ++k) {
      result.add(psi.rows[r].cols[k], psi.rows[r].vals[k] * w(r));
    }
  }
  return result;
}

}  // namespace

namespace barycentric_affine_approximator {

const std::vector<int> kInIds = {2 - 1, 3 - 1, 5 - 1, 8 - 1};
const std::vector<int> kOutIds = {1 - 1, 4 - 1, 6 - 1, 7 - 1};

void SparseVec::add(int col, double value, double eps) {
  if (std::abs(value) <= eps) {
    return;
  }
  // Merge duplicate columns immediately. Duplicates are possible when several
  // projected terms contribute to the same barycentric unknown.
  for (std::size_t i = 0; i < cols.size(); ++i) {
    if (cols[i] == col) {
      vals[i] += value;
      if (std::abs(vals[i]) <= eps) {
        cols.erase(cols.begin() + static_cast<std::ptrdiff_t>(i));
        vals.erase(vals.begin() + static_cast<std::ptrdiff_t>(i));
      }
      return;
    }
  }
  cols.push_back(col);
  vals.push_back(value);
}

double SparseVec::dot(const std::vector<double>& x) const {
  double result = 0.0;
  for (std::size_t i = 0; i < cols.size(); ++i) {
    const int col = cols[i];
    if (col < 0 || col >= static_cast<int>(x.size())) {
      throw std::runtime_error("SparseVec::dot: column outside vector size");
    }
    result += vals[i] * x[col];
  }
  return result;
}

int BarycentricVarLayout::idxX(int subsystem, int vertex_id) const {
  if (subsystem < 0 || subsystem >= kSubsystemCount) {
    throw std::invalid_argument("BarycentricVarLayout::idxX: bad subsystem");
  }
  if (vertex_id < 0 || vertex_id >= eta_s[subsystem]) {
    throw std::invalid_argument("BarycentricVarLayout::idxX: bad vertex id");
  }
  return offset_s[subsystem] + vertex_id;
}

int BarycentricVarLayout::idxY(int region, int dim) const {
  if (region < 0 || region >= num_regions) {
    throw std::invalid_argument("BarycentricVarLayout::idxY: bad region");
  }
  if (dim < 0 || dim >= kSpaceDim) {
    throw std::invalid_argument("BarycentricVarLayout::idxY: bad dimension");
  }
  return num_x + region * kSpaceDim + dim;
}

BarycentricAffineApproximator::BarycentricAffineApproximator(
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
  // This constructor mirrors GlobalAffineApproximator because the time-grid and
  // solver-reuse logic should behave the same. Only the LP variable layout is
  // different later.
  if (t_split_count < 1) {
    throw std::invalid_argument(
        "BarycentricAffineApproximator: t_split_count must be at least 1.");
  }
  if (tau_min >= tau_max || tau_min < 0.0 || tau_max < 0.0) {
    throw std::invalid_argument(
        "BarycentricAffineApproximator: tau_min must be less than tau_max and "
        "both must be nonnegative.");
  }

  max_switches_ = static_cast<int>(std::ceil(t_max_ / tau_min_));
  t_range_.resize(t_split_count);
  t_index_.resize(t_split_count);
  t_delta_ = t_split_count > 1 ? std::abs(t_max / (t_split_count - 1)) : 0.0;
  for (int i = 0; i < t_split_count; ++i) {
    t_range_[i] = (t_split_count == 1) ? t_max : i * t_delta_;
    t_index_[i] = i;
  }

  // Cube vertices are not used by the main LP, but we keep them for the future
  // border-condition implementation so the class shape stays close to the
  // global affine solver.
  const int n_vertices = 1 << kSpaceDim;
  cube_angle_vertices_.clear();
  for (int vert = 0; vert < n_vertices; ++vert) {
    Eigen::VectorXd v(kSpaceDim);
    for (int d = 0; d < kSpaceDim; ++d) {
      v(d) = ((vert & (1 << d)) != 0) ? system_params_.N : 0.0;
    }
    cube_angle_vertices_.push_back(v);
  }

  theta_t_index_lists_ = interval_building::buildThetaToTIndexLists(
      t_max_, tau_min_, tau_max_, t_range_, max_switches_);

  logger_ = spdlog::get("barycentric_affine_approximator");
  if (!logger_) {
    logger_ = spdlog::stdout_color_mt("barycentric_affine_approximator");
  }
  logger_->set_level(spdlog::level::info);

  if (highs_verbose_) {
    interval_building::prettyPrintThetaTLists(theta_t_index_lists_, t_range_);
  }
}

void BarycentricAffineApproximator::dumpInitParamsToJson(
    const std::string& filepath) const {
  // This dump mirrors the global affine class so output folders remain
  // inspectable in the same way once the barycentric path is runnable.
  const SystemParams& p = system_params_;
  const char* mode_str
      = approximation_mode_ == ApproximationMode::Upper ? "upper" : "lower";
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
                             + filepath);
  }
  f << out.str();
}

double BarycentricAffineApproximator::getBetaParamForAxis(int i, int j) const {
  // Same beta lookup as the global affine implementation. Arguments are
  // zero-based axis ids even though parameter names use paper notation.
  if (i == 5 - 1 && j == 1 - 1) {
    return system_params_.b51;
  }
  if (i == 5 - 1 && j == 7 - 1) {
    return system_params_.b57;
  }
  if (i == 8 - 1 && j == 4 - 1) {
    return system_params_.b84;
  }
  if (i == 8 - 1 && j == 6 - 1) {
    return system_params_.b86;
  }
  if (i == 3 - 1 && j == 1 - 1) {
    return system_params_.b31;
  }
  if (i == 3 - 1 && j == 6 - 1) {
    return system_params_.b36;
  }
  if (i == 2 - 1 && j == 4 - 1) {
    return system_params_.b24;
  }
  if (i == 2 - 1 && j == 7 - 1) {
    return system_params_.b27;
  }
  throw std::invalid_argument(std::format(
      "No beta parameter for zero-based movement ({}, {})", i, j));
}

std::pair<double, double> BarycentricAffineApproximator::getFMinMaxForAxis(
    int i) const {
  // Incoming-flow uncertainty bounds exist only for incoming axes 2,3,5,8 in
  // paper notation, i.e. zero-based axes 1,2,4,7 here.
  if (i == 2 - 1) {
    return std::make_pair(system_params_.f2min, system_params_.f2max);
  }
  if (i == 3 - 1) {
    return std::make_pair(system_params_.f3min, system_params_.f3max);
  }
  if (i == 5 - 1) {
    return std::make_pair(system_params_.f5min, system_params_.f5max);
  }
  if (i == 8 - 1) {
    return std::make_pair(system_params_.f8min, system_params_.f8max);
  }
  throw std::invalid_argument("No incoming-flow bounds for axis "
                              + std::to_string(i));
}

void BarycentricAffineApproximator::getIntersectionPoints() {
  // Geometry ingestion for the barycentric LP. The triangle variant is required
  // because it returns both full 8D region vertices and the simplex ids that
  // define each region's local barycentric chart.
  hcpwa::TriangleAreasVerticesResult areas_vertices
      = hcpwa::compute_triangle_areas_vertices(
          system_params_.N, system_params_.F, system_params_.v,
          system_params_.w, system_params_.b51, system_params_.b57,
          system_params_.b84, system_params_.b86, system_params_.b31,
          system_params_.b36, system_params_.b24, system_params_.b27,
          system_params_.f2min, system_params_.f3min, system_params_.f5min,
          system_params_.f8min, system_params_.f2max, system_params_.f3max,
          system_params_.f5max, system_params_.f8max, true);

  auto build_projection_layer =
      [](const std::vector<hcpwa::TriangleWithUniqueVertices>& triangles,
         std::array<int, 2> axes) {
        ProjectionLayer layer;
        layer.axes = axes;
        layer.triangles = triangles;
        layer.bases.reserve(triangles.size());

        auto get_or_add_vertex = [&layer](const Eigen::Vector2d& point) {
          // Deliberately use a simple linear tolerance search. The number of 2D
          // arrangement vertices is modest, and this is much easier to audit
          // than a hash whose rounding convention might hide indexing bugs.
          for (int i = 0; i < static_cast<int>(layer.unique_vertices.size());
               ++i) {
            if ((layer.unique_vertices[i] - point).norm() <= kGeomEps) {
              return i;
            }
          }
          layer.unique_vertices.push_back(point);
          return static_cast<int>(layer.unique_vertices.size() - 1);
        };

        for (const auto& triangle : layer.triangles) {
          TriangleBasis basis;
          const std::array<Eigen::Vector2d, 3> vertices = {
              toEigen2(triangle.a), toEigen2(triangle.b), toEigen2(triangle.c)};

          // Each shared 2D vertex must map to one barycentric unknown. The three
          // ids below are therefore deduplicated ids, not per-triangle ids.
          for (int l = 0; l < 3; ++l) {
            basis.vertex_ids[l] = get_or_add_vertex(vertices[l]);
          }

          // Barycentric coordinates are alpha(z) = H z + h. We compute them by
          // inverting G * alpha = [z_x, z_y, 1]^T, where G columns are the three
          // simplex vertices with a final row of ones.
          Eigen::Matrix3d G;
          G << vertices[0](0), vertices[1](0), vertices[2](0), vertices[0](1),
              vertices[1](1), vertices[2](1), 1.0, 1.0, 1.0;
          const double det = G.determinant();
          if (std::abs(det) <= kGeomEps) {
            throw std::runtime_error(
                "Degenerate triangle in barycentric geometry.");
          }
          const Eigen::Matrix3d M = G.inverse();
          basis.H = M.block<3, 2>(0, 0);
          basis.h = M.col(2);

          // Validate alpha(g_l) = e_l. This catches transposed H/h formulas and
          // duplicated vertices before they contaminate LP rows.
          for (int l = 0; l < 3; ++l) {
            const Eigen::Vector3d alpha = basis.H * vertices[l] + basis.h;
            for (int r = 0; r < 3; ++r) {
              const double expected = (r == l) ? 1.0 : 0.0;
              if (std::abs(alpha(r) - expected) > 1e-6) {
                throw std::runtime_error(
                    "Barycentric basis failed vertex reconstruction.");
              }
            }
            if (std::abs(alpha.sum() - 1.0) > 1e-6) {
              throw std::runtime_error(
                  "Barycentric basis coordinates do not sum to one.");
            }
          }

          layer.bases.push_back(std::move(basis));
        }
        return layer;
      };

  auto convert_regions =
      [](const std::vector<std::vector<hcpwa::Vec<8>>>& src) {
        std::vector<std::vector<Eigen::VectorXd>> dst;
        dst.resize(src.size());
        for (int j = 0; j < static_cast<int>(src.size()); ++j) {
          dst[j].reserve(src[j].size());
          for (const auto& vertex : src[j]) {
            dst[j].push_back(toEigen8(vertex));
          }
        }
        return dst;
      };

  auto convert_region_indices =
      [](const std::vector<std::vector<std::size_t>>& src) {
        std::vector<std::array<int, kSubsystemCount>> dst;
        dst.reserve(src.size());
        for (const auto& tuple : src) {
          if (tuple.size() != kSubsystemCount) {
            throw std::runtime_error(
                "Each barycentric region must have exactly five simplex ids.");
          }
          std::array<int, kSubsystemCount> converted{};
          for (int s = 0; s < kSubsystemCount; ++s) {
            converted[s] = static_cast<int>(tuple[s]);
          }
          dst.push_back(converted);
        }
        return dst;
      };

  // Phase 0 must use the same order as intersection_prism_indices_phase0:
  // [31, 36, 24, 27, 58].
  auto phase0_axes = projectionAxesForPhase(0);
  phase_geometries_[0].layers[0]
      = build_projection_layer(areas_vertices.triangles31, phase0_axes[0]);
  phase_geometries_[0].layers[1]
      = build_projection_layer(areas_vertices.triangles36, phase0_axes[1]);
  phase_geometries_[0].layers[2]
      = build_projection_layer(areas_vertices.triangles24, phase0_axes[2]);
  phase_geometries_[0].layers[3]
      = build_projection_layer(areas_vertices.triangles27, phase0_axes[3]);
  phase_geometries_[0].layers[4]
      = build_projection_layer(areas_vertices.triangles58, phase0_axes[4]);
  phase_geometries_[0].region_vertices
      = convert_regions(areas_vertices.intersection_points_phase0);
  phase_geometries_[0].region_triangle_ids
      = convert_region_indices(areas_vertices.intersection_prism_indices_phase0);

  // Phase 1 must use the same order as intersection_prism_indices_phase1:
  // [51, 57, 84, 86, 23].
  auto phase1_axes = projectionAxesForPhase(1);
  phase_geometries_[1].layers[0]
      = build_projection_layer(areas_vertices.triangles51, phase1_axes[0]);
  phase_geometries_[1].layers[1]
      = build_projection_layer(areas_vertices.triangles57, phase1_axes[1]);
  phase_geometries_[1].layers[2]
      = build_projection_layer(areas_vertices.triangles84, phase1_axes[2]);
  phase_geometries_[1].layers[3]
      = build_projection_layer(areas_vertices.triangles86, phase1_axes[3]);
  phase_geometries_[1].layers[4]
      = build_projection_layer(areas_vertices.triangles23, phase1_axes[4]);
  phase_geometries_[1].region_vertices
      = convert_regions(areas_vertices.intersection_points_phase1);
  phase_geometries_[1].region_triangle_ids
      = convert_region_indices(areas_vertices.intersection_prism_indices_phase1);

  for (int phase = 0; phase < kPhases; ++phase) {
    const auto& geometry = phase_geometries_[phase];
    if (geometry.region_vertices.size() != geometry.region_triangle_ids.size()) {
      throw std::runtime_error(
          "Barycentric geometry has mismatched region vertices and ids.");
    }

    BarycentricVarLayout layout;
    layout.num_regions = static_cast<int>(geometry.region_vertices.size());
    for (int s = 0; s < kSubsystemCount; ++s) {
      layout.offset_s[s] = layout.num_x;
      layout.eta_s[s]
          = static_cast<int>(geometry.layers[s].unique_vertices.size());
      layout.num_x += layout.eta_s[s];
    }
    layout.num_cols = layout.num_x + kSpaceDim * layout.num_regions;
    layouts_[phase] = layout;
  }

  computeNodeWeights();
}

void BarycentricAffineApproximator::computeNodeWeights() {
  // w_(s,k) = integral over Omega of the hat function of node (s,k).
  // Omega = [0,N]^8 and P_s keeps two coordinates, so the integral factorizes:
  //   int_Omega alpha^(s)_k(P_s n) dn = N^6 * int_{[0,N]^2} alpha^(s)_k(u) du,
  // and the integral of a barycentric coordinate over its triangle is
  // area/3 (step 2.2, section 7).
  const double n6 = std::pow(system_params_.N, kSpaceDim - 2);
  const double expected = static_cast<double>(kSubsystemCount)
                          * std::pow(system_params_.N, kSpaceDim);

  std::size_t total_triangles = 0;
  for (int phase = 0; phase < kPhases; ++phase) {
    for (int s = 0; s < kSubsystemCount; ++s) {
      total_triangles += phase_geometries_[phase].layers[s].triangles.size();
    }
  }
  logger_->info(
      "computeNodeWeights: starting over {} phases, {} total triangles",
      kPhases, total_triangles);

  std::size_t processed_triangles = 0;
  const auto started_at = std::chrono::steady_clock::now();
  auto next_progress_at = started_at;
  auto report_progress = [&](bool force, int phase, int subsystem) {
    const auto now = std::chrono::steady_clock::now();
    if (!force && now < next_progress_at) {
      return;
    }
    const double percent
        = total_triangles == 0
              ? 100.0
              : 100.0 * static_cast<double>(processed_triangles)
                    / static_cast<double>(total_triangles);
    const double elapsed_seconds
        = std::chrono::duration<double>(now - started_at).count();
    logger_->info(
        "computeNodeWeights progress: {:.1f}% ({}/{} triangles), "
        "phase={}/{}, subsystem={}/{}, elapsed={:.1f}s",
        percent, processed_triangles, total_triangles, phase + 1, kPhases,
        subsystem + 1, kSubsystemCount, elapsed_seconds);
    next_progress_at = now + std::chrono::seconds(5);
  };
  report_progress(true, 0, 0);

  for (int phase = 0; phase < kPhases; ++phase) {
    const auto& geometry = phase_geometries_[phase];
    const auto& layout = layouts_[phase];
    Eigen::VectorXd w = Eigen::VectorXd::Zero(layout.num_x);
    logger_->info(
        "computeNodeWeights: phase {}/{} num_x={}", phase + 1, kPhases,
        layout.num_x);

    for (int s = 0; s < kSubsystemCount; ++s) {
      const auto& layer = geometry.layers[s];
      if (layer.triangles.size() != layer.bases.size()) {
        throw std::runtime_error(
            "computeNodeWeights: triangle and basis counts differ");
      }
      logger_->info(
          "computeNodeWeights: phase {}/{} subsystem {}/{} triangles={}",
          phase + 1, kPhases, s + 1, kSubsystemCount, layer.triangles.size());
      for (std::size_t m = 0; m < layer.triangles.size(); ++m) {
        const Eigen::Vector2d a = toEigen2(layer.triangles[m].a);
        const Eigen::Vector2d b = toEigen2(layer.triangles[m].b);
        const Eigen::Vector2d c = toEigen2(layer.triangles[m].c);
        const double area
            = 0.5
              * std::abs((b(0) - a(0)) * (c(1) - a(1))
                         - (c(0) - a(0)) * (b(1) - a(1)));
        for (int local_vertex = 0; local_vertex < 3; ++local_vertex) {
          const int col
              = layout.idxX(s, layer.bases[m].vertex_ids[local_vertex]);
          w(col) += n6 * area / 3.0;
        }
        ++processed_triangles;
        if ((processed_triangles % 10000) == 0) {
          report_progress(false, phase, s);
        }
      }
      report_progress(true, phase, s);
    }

    // Control identity: every triangle is counted once per vertex, so each layer
    // integrates to |Omega| and the five layers to 5|Omega|. A mismatch means the
    // triangulation does not tile [0,N]^2 or a vertex id is wrong.
    const double w_sum = w.sum();
    const double abs_diff = std::abs(w_sum - expected);
    const double rel_error = abs_diff / expected;
    // if (abs_diff > 1e-6 * expected) {
    //   throw std::runtime_error(
    //       std::format(
    //           "computeNodeWeights: node weights failed the 1^T w = 5 N^8 identity. w.sum = {:.16f}, expected = {:.16f}, abs_diff = {:.16f}, rel_error = {:.16e}",
    //           w_sum, expected, abs_diff, rel_error)
    //       );
    // }
    logger_->info(
        "computeNodeWeights: phase {}/{} identity ok (w.sum={:.6f}, "
        "expected={:.6f}, rel_error={:.3e})",
        phase + 1, kPhases, w_sum, expected, rel_error);
    node_weights_[phase] = std::move(w);
  }
  report_progress(true, kPhases - 1, kSubsystemCount - 1);
  logger_->info("computeNodeWeights: finished");
}

SparseVec BarycentricAffineApproximator::buildPhiRow(
    int phase, int region, const Eigen::VectorXd& point,
    double tolerance) const {
  if (phase < 0 || phase >= kPhases) {
    throw std::invalid_argument("buildPhiRow: invalid phase");
  }
  if (point.size() != kSpaceDim) {
    throw std::invalid_argument("buildPhiRow: point must be 8-dimensional");
  }
  const auto& geometry = phase_geometries_[phase];
  const auto& layout = layouts_[phase];
  if (region < 0
      || region >= static_cast<int>(geometry.region_triangle_ids.size())) {
    throw std::invalid_argument("buildPhiRow: invalid region id");
  }

  SparseVec phi;
  const auto& triangle_ids = geometry.region_triangle_ids[region];
  for (int s = 0; s < kSubsystemCount; ++s) {
    const int triangle_id = triangle_ids[s];
    const auto& layer = geometry.layers[s];
    if (triangle_id < 0
        || triangle_id >= static_cast<int>(layer.bases.size())) {
      throw std::runtime_error("buildPhiRow: triangle id out of range");
    }

    const auto& basis = layer.bases[triangle_id];
    Eigen::Vector2d projected;
    projected(0) = point(layer.axes[0]);
    projected(1) = point(layer.axes[1]);
    Eigen::Vector3d alpha = basis.H * projected + basis.h;

    // A value row is only meaningful in the selected local simplex. The
    // tolerance allows points on shared triangle/area boundaries, but rejects
    // genuine tuple-order mistakes that would project outside the simplex.
    if (std::abs(alpha.sum() - 1.0) > tolerance) {
      throw std::runtime_error("buildPhiRow: barycentric alpha sum is not one");
    }
    for (int local_vertex = 0; local_vertex < 3; ++local_vertex) {
      if (alpha(local_vertex) < -tolerance
          || alpha(local_vertex) > 1.0 + tolerance) {
        throw std::runtime_error(
            "buildPhiRow: point is outside selected simplex");
      }
      const int x_col = layout.idxX(s, basis.vertex_ids[local_vertex]);
      // Preserve near-boundary barycentric coordinates. SparseVec::add defaults
      // to kEps (1e-5), which is larger than the geometry tolerance and can
      // remove enough mass for the five-layer phi row to stop summing to five.
      phi.add(x_col, alpha(local_vertex), kGeomEps);
    }
  }
  return phi;
}

std::vector<int> BarycentricAffineApproximator::locateRegions(
    int phase, const Eigen::VectorXd& point, double tolerance) const {
  if (phase < 0 || phase >= kPhases) {
    throw std::invalid_argument("locateRegions: invalid phase");
  }
  if (point.size() != kSpaceDim) {
    throw std::invalid_argument("locateRegions: point must be 8-dimensional");
  }

  std::vector<int> matches;
  const auto& geometry = phase_geometries_[phase];
  for (int region = 0;
       region < static_cast<int>(geometry.region_triangle_ids.size());
       ++region) {
    bool contains = true;
    const auto& triangle_ids = geometry.region_triangle_ids[region];
    for (int s = 0; s < kSubsystemCount && contains; ++s) {
      const int triangle_id = triangle_ids[s];
      const auto& layer = geometry.layers[s];
      if (triangle_id < 0
          || triangle_id >= static_cast<int>(layer.bases.size())) {
        throw std::runtime_error("locateRegions: triangle id out of range");
      }
      const auto& basis = layer.bases[triangle_id];
      Eigen::Vector2d projected;
      projected(0) = point(layer.axes[0]);
      projected(1) = point(layer.axes[1]);
      const Eigen::Vector3d alpha = basis.H * projected + basis.h;
      if (std::abs(alpha.sum() - 1.0) > tolerance) {
        contains = false;
        break;
      }
      for (int local_vertex = 0; local_vertex < 3; ++local_vertex) {
        if (alpha(local_vertex) < -tolerance
            || alpha(local_vertex) > 1.0 + tolerance) {
          contains = false;
          break;
        }
      }
    }
    if (contains) {
      matches.push_back(region);
    }
  }
  return matches;
}

double BarycentricAffineApproximator::evaluateBarycentricValue(
    int phase, const std::vector<double>& x, const Eigen::VectorXd& point,
    double tolerance) const {
  const auto& layout = layouts_[phase];
  if (x.size() != static_cast<std::size_t>(layout.num_x)) {
    throw std::invalid_argument("evaluateBarycentricValue: x size must be "
                                + std::to_string(layout.num_x));
  }

  const std::vector<int> matches = locateRegions(phase, point, tolerance);
  if (matches.empty()) {
    throw std::runtime_error(
        "evaluateBarycentricValue: point is outside phase partition");
  }

  double min_value = std::numeric_limits<double>::infinity();
  double max_value = -std::numeric_limits<double>::infinity();
  for (int region : matches) {
    const double value = buildPhiRow(phase, region, point, tolerance).dot(x);
    min_value = std::min(min_value, value);
    max_value = std::max(max_value, value);
  }

  if (max_value - min_value > 10.0 * tolerance) {
    throw std::runtime_error(
        "evaluateBarycentricValue: inconsistent boundary-region values");
  }
  return max_value;
}

std::vector<int> BarycentricAffineApproximator::admissibleThetaIds(
    double theta) const {
  const double theta_min = std::min(theta + tau_min_, t_max_);
  const double theta_max = std::min(theta + tau_max_, t_max_);
  const double half_step = t_delta_ / 2.0;

  std::vector<int> result;
  for (std::size_t i = 0; i < t_range_.size(); ++i) {
    if (t_range_[i] >= theta_min - half_step
        && t_range_[i] <= theta_max + half_step) {
      result.push_back(t_index_[i]);
    }
  }
  if (result.empty()) {
    throw std::runtime_error("admissibleThetaIds: no theta range found");
  }
  return result;
}

std::pair<Eigen::RowVectorXd, Eigen::RowVectorXd>
BarycentricAffineApproximator::getFIJMinResolution(
    int i, int j, const Eigen::VectorXd& n) const {
  // Resolves f_ij(n) = min{ beta F, beta v n_i, w(N - n_j) } at representative
  // point n. The returned pair is (row, scalar) such that f_ij(n) = row*n +
  // scalar on this region, matching the global affine implementation.
  const double N = system_params_.N;
  const double F = system_params_.F;
  const double v = system_params_.v;
  const double w = system_params_.w;
  const double beta_i_j = getBetaParamForAxis(i, j);

  const double free_flow = beta_i_j * v * n(i);
  const double capacity = beta_i_j * F;
  const double receiving = w * (N - n(j));

  Eigen::RowVectorXd row = Eigen::RowVectorXd::Zero(kSpaceDim);
  Eigen::RowVectorXd scalar = Eigen::RowVectorXd::Zero(1);

  if (capacity < free_flow + kEps && capacity < receiving + kEps) {
    scalar(0) = capacity;
  } else if (receiving < free_flow + kEps && receiving < capacity + kEps) {
    row(j) = -w;
    scalar(0) = w * N;
  } else if (free_flow < capacity + kEps && free_flow < receiving + kEps) {
    row(i) = beta_i_j * v;
  } else {
    throw std::invalid_argument(std::format(
        "No min branch for f_{}{}: free={}, cap={}, recv={}", i, j, free_flow,
        capacity, receiving));
  }
  return std::make_pair(row, scalar);
}

Eigen::VectorXd BarycentricAffineApproximator::areaCentroidCoords(
    int j, int phase) const {
  // The centroid is only a representative point for branch resolution. The LP
  // still enforces residual constraints at every stored region vertex.
  if (phase < 0 || phase >= kPhases) {
    throw std::invalid_argument("areaCentroidCoords: invalid phase");
  }
  const auto& regions = phase_geometries_[phase].region_vertices;
  if (j < 0 || j >= static_cast<int>(regions.size())) {
    throw std::invalid_argument("areaCentroidCoords: invalid region id");
  }
  if (regions[j].empty()) {
    throw std::runtime_error("areaCentroidCoords: region has no vertices");
  }
  Eigen::VectorXd centroid = Eigen::VectorXd::Zero(kSpaceDim);
  for (const auto& vertex : regions[j]) {
    centroid += vertex;
  }
  centroid /= static_cast<double>(regions[j].size());
  hcpwa::util::assertShape(centroid, kSpaceDim);
  return centroid;
}

std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd, double>
BarycentricAffineApproximator::getAMatrFVecGVecAndGScalJ(
    int j, int phase) const {
  // Copy of the global affine CTM branch-resolution formulas. The barycentric
  // LP changes the value-function representation, not the CTM dynamics.
  Eigen::VectorXd n = areaCentroidCoords(j, phase);

  Eigen::MatrixXd a_matr;
  Eigen::VectorXd b_vec;
  Eigen::VectorXd g_vec;
  double g_scal = 0.0;

  if (phase == 0) {
    auto [f31_matr_row, f31_vec_row] = getFIJMinResolution(3 - 1, 1 - 1, n);
    auto [f36_matr_row, f36_vec_row] = getFIJMinResolution(3 - 1, 6 - 1, n);
    auto [f24_matr_row, f24_vec_row] = getFIJMinResolution(2 - 1, 4 - 1, n);
    auto [f27_matr_row, f27_vec_row] = getFIJMinResolution(2 - 1, 7 - 1, n);

    a_matr.resize(kSpaceDim, kSpaceDim);
    a_matr.row(0) = f31_matr_row;
    a_matr.row(1) = -f24_matr_row - f27_matr_row;
    a_matr.row(2) = -f31_matr_row - f36_matr_row;
    a_matr.row(3) = f24_matr_row;
    a_matr.row(4) = Eigen::RowVectorXd::Zero(kSpaceDim);
    a_matr.row(5) = f36_matr_row;
    a_matr.row(6) = f27_matr_row;
    a_matr.row(7) = Eigen::RowVectorXd::Zero(kSpaceDim);

    b_vec.resize(kSpaceDim);
    b_vec(0) = f31_vec_row(0);
    b_vec(1) = -f24_vec_row(0) - f27_vec_row(0);
    b_vec(2) = -f31_vec_row(0) - f36_vec_row(0);
    b_vec(3) = f24_vec_row(0);
    b_vec(4) = 0.0;
    b_vec(5) = f36_vec_row(0);
    b_vec(6) = f27_vec_row(0);
    b_vec(7) = 0.0;

    g_vec = (f31_matr_row + f36_matr_row + f24_matr_row + f27_matr_row)
                .transpose();
    g_scal = f31_vec_row(0) + f36_vec_row(0) + f24_vec_row(0)
             + f27_vec_row(0);
  } else if (phase == 1) {
    auto [f51_matr_row, f51_vec_row] = getFIJMinResolution(5 - 1, 1 - 1, n);
    auto [f57_matr_row, f57_vec_row] = getFIJMinResolution(5 - 1, 7 - 1, n);
    auto [f84_matr_row, f84_vec_row] = getFIJMinResolution(8 - 1, 4 - 1, n);
    auto [f86_matr_row, f86_vec_row] = getFIJMinResolution(8 - 1, 6 - 1, n);

    a_matr.resize(kSpaceDim, kSpaceDim);
    a_matr.row(0) = f51_matr_row;
    a_matr.row(1) = Eigen::RowVectorXd::Zero(kSpaceDim);
    a_matr.row(2) = Eigen::RowVectorXd::Zero(kSpaceDim);
    a_matr.row(3) = f84_matr_row;
    a_matr.row(4) = -f51_matr_row - f57_matr_row;
    a_matr.row(5) = f86_matr_row;
    a_matr.row(6) = f57_matr_row;
    a_matr.row(7) = -f84_matr_row - f86_matr_row;

    b_vec.resize(kSpaceDim);
    b_vec(0) = f51_vec_row(0);
    b_vec(1) = 0.0;
    b_vec(2) = 0.0;
    b_vec(3) = f84_vec_row(0);
    b_vec(4) = -f51_vec_row(0) - f57_vec_row(0);
    b_vec(5) = f86_vec_row(0);
    b_vec(6) = f57_vec_row(0);
    b_vec(7) = -f84_vec_row(0) - f86_vec_row(0);

    g_vec = (f51_matr_row + f57_matr_row + f84_matr_row + f86_matr_row)
                .transpose();
    g_scal = f51_vec_row(0) + f57_vec_row(0) + f84_vec_row(0)
             + f86_vec_row(0);
  } else {
    throw std::invalid_argument("getAMatrFVecGVecAndGScalJ: invalid phase");
  }

  hcpwa::util::assertShape(a_matr, kSpaceDim, kSpaceDim);
  hcpwa::util::assertShape(b_vec, kSpaceDim);
  hcpwa::util::assertShape(g_vec, kSpaceDim);
  hcpwa::util::assertScalar(g_scal);
  return std::make_tuple(a_matr, b_vec, g_vec, g_scal);
}

std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd, Eigen::VectorXd>
BarycentricAffineApproximator::getQQForArea(int j, int phase) const {
  // Builds affine center/radius maps for the disturbance box:
  //   c_j(n) = Qc_j n + qc_j
  //   rho_j(n) = Qr_j n + qr_j.
  // The branch resolution is intentionally copied from the global affine code.
  Eigen::VectorXd n0 = areaCentroidCoords(j, phase);
  Eigen::MatrixXd Q_upper = Eigen::MatrixXd::Zero(kSpaceDim, kSpaceDim);
  Eigen::VectorXd q_upper = Eigen::VectorXd::Zero(kSpaceDim);
  Eigen::MatrixXd Q_lower = Eigen::MatrixXd::Zero(kSpaceDim, kSpaceDim);
  Eigen::VectorXd q_lower = Eigen::VectorXd::Zero(kSpaceDim);

  const double N = system_params_.N;
  const double w = system_params_.w;
  const double v = system_params_.v;
  const double F = system_params_.F;

  for (int i : kInIds) {
    const double ni = n0(i);
    auto [f_min, f_max] = getFMinMaxForAxis(i);
    if (f_min < w * (N - ni)) {
      q_lower(i) = f_min;
    } else {
      Q_lower(i, i) = -w;
      q_lower(i) = w * N;
    }
    if (f_max < w * (N - ni)) {
      q_upper(i) = f_max;
    } else {
      Q_upper(i, i) = -w;
      q_upper(i) = w * N;
    }
  }

  for (int i : kOutIds) {
    const double ni = n0(i);
    if (F < v * ni) {
      q_lower(i) = -F;
    } else {
      Q_lower(i, i) = -v;
    }
  }

  Eigen::MatrixXd Q_c = (Q_upper + Q_lower) / 2.0;
  Eigen::MatrixXd Q_r = (Q_upper - Q_lower) / 2.0;
  Eigen::VectorXd q_c = (q_upper + q_lower) / 2.0;
  Eigen::VectorXd q_r = (q_upper - q_lower) / 2.0;
  return std::make_tuple(Q_c, q_c, Q_r, q_r);
}

std::tuple<std::vector<Eigen::MatrixXd>, std::vector<Eigen::VectorXd>,
           std::vector<Eigen::MatrixXd>, std::vector<Eigen::VectorXd>,
           std::vector<Eigen::MatrixXd>, std::vector<Eigen::VectorXd>,
           std::vector<Eigen::VectorXd>, std::vector<double>>
BarycentricAffineApproximator::precomputeSystemMatrices(int phase) {
  // Static per-region CTM data. These matrices are independent of time and of
  // x_next, so they are computed once before the reusable LP is built.
  const int n_areas
      = static_cast<int>(phase_geometries_[phase].region_vertices.size());
  std::vector<Eigen::MatrixXd> A_j_matrs;
  std::vector<Eigen::VectorXd> f_j_vecs;
  std::vector<Eigen::MatrixXd> Q_c_j_matrs;
  std::vector<Eigen::VectorXd> q_c_j_vecs;
  std::vector<Eigen::MatrixXd> Q_r_j_matrs;
  std::vector<Eigen::VectorXd> q_r_j_vecs;
  std::vector<Eigen::VectorXd> g_j_vecs;
  std::vector<double> g_j_scals;

  A_j_matrs.reserve(n_areas);
  f_j_vecs.reserve(n_areas);
  Q_c_j_matrs.reserve(n_areas);
  q_c_j_vecs.reserve(n_areas);
  Q_r_j_matrs.reserve(n_areas);
  q_r_j_vecs.reserve(n_areas);
  g_j_vecs.reserve(n_areas);
  g_j_scals.reserve(n_areas);

  for (int j = 0; j < n_areas; ++j) {
    auto [A_j, f_j, g_j, g0_j] = getAMatrFVecGVecAndGScalJ(j, phase);
    auto [Qc_j, qc_j, Qr_j, qr_j] = getQQForArea(j, phase);
    A_j_matrs.push_back(std::move(A_j));
    f_j_vecs.push_back(std::move(f_j));
    Q_c_j_matrs.push_back(std::move(Qc_j));
    q_c_j_vecs.push_back(std::move(qc_j));
    Q_r_j_matrs.push_back(std::move(Qr_j));
    q_r_j_vecs.push_back(std::move(qr_j));
    g_j_vecs.push_back(std::move(g_j));
    g_j_scals.push_back(g0_j);
  }

  return std::make_tuple(A_j_matrs, f_j_vecs, Q_c_j_matrs, q_c_j_vecs,
                         Q_r_j_matrs, q_r_j_vecs, g_j_vecs, g_j_scals);
}

std::tuple<std::vector<int>, std::vector<int>, std::vector<double>,
           std::vector<double>, std::vector<double>, Eigen::RowVectorXd>
BarycentricAffineApproximator::prepareLpMatrices(int phase) {
  // Builds the static LP matrix for one phase. Only residual row upper bounds
  // change with x_next; all left-hand side coefficients and objective costs are
  // fixed after this function returns.
  const auto& geometry = phase_geometries_[phase];
  const auto& layout = layouts_[phase];
  const int n_cols = layout.num_cols;

  std::vector<int> starts = {0};
  std::vector<int> col_index;
  std::vector<double> value;
  std::vector<double> row_lower;
  std::vector<double> row_upper;
  Eigen::RowVectorXd c_vec = Eigen::RowVectorXd::Zero(n_cols);

  // s = +1 for the upper bound, -1 for the lower one. After the corrections of
  // step 2.1 this is the only thing separating the two directions.
  const double s = signS();

  rhs_terms_[phase].clear();
  psi_by_region_[phase].clear();
  psi_by_region_[phase].resize(layout.num_regions);

  const auto& A_j_matrs = A_j_matrs_[phase];
  const auto& f_j_vecs = f_j_vecs_[phase];
  const auto& Q_c_j_matrs = Q_c_j_matrs_[phase];
  const auto& q_c_j_vecs = q_c_j_vecs_[phase];
  const auto& Q_r_j_matrs = Q_r_j_matrs_[phase];
  const auto& q_r_j_vecs = q_r_j_vecs_[phase];
  const auto& g_j_vecs = g_j_vecs_[phase];
  const auto& g_j_scals = g_j_scals_[phase];

  for (int j = 0; j < layout.num_regions; ++j) {
    const auto& triangle_ids = geometry.region_triangle_ids[j];

    // Build Psi_j once per full 8D region. For every local barycentric value
    // v_l, derivative d alpha_l / d n_axis contributes H(l,local_axis).
    SparsePsi psi_j;
    for (int s = 0; s < kSubsystemCount; ++s) {
      const int triangle_id = triangle_ids[s];
      const auto& layer = geometry.layers[s];
      if (triangle_id < 0
          || triangle_id >= static_cast<int>(layer.bases.size())) {
        throw std::runtime_error("prepareLpMatrices: triangle id out of range");
      }
      const auto& basis = layer.bases[triangle_id];
      for (int local_vertex = 0; local_vertex < 3; ++local_vertex) {
        const int x_col = layout.idxX(s, basis.vertex_ids[local_vertex]);
        for (int local_axis = 0; local_axis < 2; ++local_axis) {
          const int state_axis = layer.axes[local_axis];
          psi_j.rows[state_axis].add(x_col,
                                     basis.H(local_vertex, local_axis));
        }
      }
    }
    psi_by_region_[phase][j] = psi_j;

    for (const Eigen::VectorXd& nu : geometry.region_vertices[j]) {
      // phi_{j,nu} is the value selector: V(t,nu) = phi_{j,nu}^T x. Each
      // projection contributes the three barycentric coordinates of P_s nu.
      SparseVec phi = buildPhiRow(phase, j, nu, kEps);

      const Eigen::VectorXd d = A_j_matrs[j] * nu + f_j_vecs[j];
      const Eigen::VectorXd c = Q_c_j_matrs[j] * nu + q_c_j_vecs[j];
      Eigen::VectorXd rho = Q_r_j_matrs[j] * nu + q_r_j_vecs[j];

      // rho is the radius of the uncertainty box. Small negative values can
      // only be numerical noise; larger negatives mean the region did not
      // refine uncertainty branches correctly.
      for (int r = 0; r < kSpaceDim; ++r) {
        if (rho(r) < -kEps) {
          throw std::runtime_error(
              "prepareLpMatrices: negative uncertainty radius.");
        }
        if (rho(r) < 0.0) {
          rho(r) = 0.0;
        }
      }

      // m = A_j nu + f_j + c_j(nu): the drift with the disturbance-box centre
      // folded in. Used by the Right row, whose RHS is recomputed per step.
      const Eigen::VectorXd m = d + c;
      const double g_nu = g_j_vecs[j].dot(nu) + g_j_scals[j];

      // ---- Row (L): the segment endpoint with the unknown value z. ----
      // F^L = a^T z + b + s rho^T |Psi_j z|, with
      //   a = Psi_j^T m - phi / dt,  b = phi^T x_next / dt + g(nu).
      // The modulus is taken at the unknown, so it is lifted to y_j >= |Psi_j z|
      // and the row reads  s a^T z + rho^T y_j <= -s b.  Note that rho enters
      // with a plus for BOTH directions; only a and the RHS flip (step 2.1, 4).
      SparseVec a = psiTransposeTimes(psi_j, m);
      for (std::size_t k = 0; k < phi.cols.size(); ++k) {
        a.add(phi.cols[k], -phi.vals[k] / t_delta_);
      }

      SparseVec left_row;
      for (std::size_t k = 0; k < a.cols.size(); ++k) {
        left_row.add(a.cols[k], s * a.vals[k]);
      }
      for (int r = 0; r < kSpaceDim; ++r) {
        left_row.add(layout.idxY(j, r), rho(r));
      }
      appendSparseRow(starts, col_index, value, left_row, kEps);
      row_lower.push_back(-std::numeric_limits<double>::infinity());
      // Fixed part of -s b; the -s phi^T x_next / dt part is added per step.
      row_upper.push_back(-s * g_nu);
      {
        ResidualRhsTerm term;
        term.row_id = static_cast<int>(row_upper.size() - 1);
        term.kind = RhsKind::Left;
        term.phi = phi;
        rhs_terms_[phase].push_back(std::move(term));
      }

      // ---- Row (R): the segment endpoint with the known value x_next. ----
      // F^R = phi^T (x_next - z) / dt + beta, with
      //   beta = q^T m + s rho^T |q| + g(nu),  q = Psi_j x_next.
      // Multiplying by dt > 0 gives the row  -s phi^T z <= -s (phi^T x_next +
      // dt beta).  There is no y block here: the modulus is taken at the known
      // x_next and therefore evaluates to a number (step 2.1, 3.1).
      SparseVec right_row;
      for (std::size_t k = 0; k < phi.cols.size(); ++k) {
        right_row.add(phi.cols[k], -s * phi.vals[k]);
      }
      appendSparseRow(starts, col_index, value, right_row, kEps);
      row_lower.push_back(-std::numeric_limits<double>::infinity());
      // Entirely dynamic: recomputed from scratch in updateHighsRhsUpperBounds.
      row_upper.push_back(0.0);
      {
        ResidualRhsTerm term;
        term.row_id = static_cast<int>(row_upper.size() - 1);
        term.kind = RhsKind::Right;
        term.phi = phi;
        term.region = j;
        term.m = m;
        term.rho = rho;
        term.g = g_nu;
        rhs_terms_[phase].push_back(std::move(term));
      }

      // Objective: l1 norm of the residuals at the endpoint with the known
      // value. That residual is linear in z (no modulus, no y), so its l1 norm
      // is an exact LP objective; the residual at the other endpoint is concave
      // in z and cannot be minimized by an LP (step 2.1, section 5).
      //   delta^R     = -s F^R = (s/dt) phi^T z + const,
      //   sum delta^R = (s/dt) (sum phi)^T z + const,   minimized.
      for (std::size_t k = 0; k < phi.cols.size(); ++k) {
        c_vec(phi.cols[k]) += (s / t_delta_) * phi.vals[k];
      }
    }

    // Absolute value linearization for this region:
    //   Psi_j x - y_j <= 0
    //  -Psi_j x - y_j <= 0
    // y_j >= 0 is implemented as a column lower bound in initializeHighs().
    for (int r = 0; r < kSpaceDim; ++r) {
      SparseVec row = psi_j.rows[r];
      row.add(layout.idxY(j, r), -1.0);
      appendSparseRow(starts, col_index, value, row, kEps);
      row_lower.push_back(-std::numeric_limits<double>::infinity());
      row_upper.push_back(0.0);
    }
    for (int r = 0; r < kSpaceDim; ++r) {
      SparseVec row;
      for (std::size_t k = 0; k < psi_j.rows[r].cols.size(); ++k) {
        row.add(psi_j.rows[r].cols[k], -psi_j.rows[r].vals[k]);
      }
      row.add(layout.idxY(j, r), -1.0);
      appendSparseRow(starts, col_index, value, row, kEps);
      row_lower.push_back(-std::numeric_limits<double>::infinity());
      row_upper.push_back(0.0);
    }
  }

  const int n_rows = static_cast<int>(row_upper.size());
  const int nnz = static_cast<int>(value.size());
  logger_->info(
      "prepareLpMatrices: phase={}, variables(num_cols)={}, num_x={}, "
      "num_y={}, regions={}, rows={}, nnz={}, rhs_terms={}",
      phase, n_cols, layout.num_x, kSpaceDim * layout.num_regions,
      layout.num_regions, n_rows, nnz, rhs_terms_[phase].size());

  return std::make_tuple(std::move(starts), std::move(col_index),
                         std::move(value), std::move(row_lower),
                         std::move(row_upper), std::move(c_vec));
}

std::vector<double> BarycentricAffineApproximator::getBorderConditions(
    int switch_phase, int theta_idx, double theta, int switch_cnt) const {
  // The run loop keeps the historical "switch_phase" name from the global
  // affine code, where it denotes the phase reached after switching. For the
  // border LP that is the source phase i_src; the vector returned here must be
  // expressed in the target phase basis currently being initialized.
  if (switch_phase < 0 || switch_phase >= kPhases) {
    throw std::invalid_argument("getBorderConditions: invalid switch_phase");
  }
  const int source_phase = switch_phase;
  const int target_phase = 1 - switch_phase;
  const auto& target_layout = layouts_[target_phase];
  const auto& source_layout = layouts_[source_phase];

  if (theta_idx < 0 || theta_idx >= static_cast<int>(t_range_.size())) {
    throw std::invalid_argument("getBorderConditions: invalid theta_idx");
  }
  if (std::abs(t_range_[theta_idx] - theta) > kEps) {
    throw std::runtime_error(
        "getBorderConditions: theta does not match theta_idx");
  }

  // Terminal switching family: by the current mathematical convention the
  // terminal value is zero, so no fitting LP is needed for p = 0.
  if (switch_cnt == 0) {
    return std::vector<double>(target_layout.num_x, 0.0);
  }
  if (switch_cnt < 0) {
    throw std::invalid_argument("getBorderConditions: negative switch_cnt");
  }

  const std::vector<int> theta_end_ids = admissibleThetaIds(theta);
  const bool is_upper = approximation_mode_ == ApproximationMode::Upper;

  // The family of already computed members, fetched once. Two reasons to hoist
  // it out of the vertex loop instead of refetching per vertex:
  //   - ValueFunction::get returns by value under a mutex, so a per-vertex fetch
  //     would take |vertices| * |members| locked copies of an eta-vector and
  //     serialize the border LPs of all worker threads;
  //   - the member index c then denotes the same (level, theta*) pair for every
  //     vertex by construction, which is what the per-cell selection below
  //     relies on.
  std::vector<std::vector<double>> candidates;
  for (int r = 0; r < switch_cnt; ++r) {
    for (int theta_end_idx : theta_end_ids) {
      if (!value_function_.contains(source_phase, r, theta_idx,
                                    theta_end_idx)) {
        // Missing candidates are skipped explicitly and never replaced by zero:
        // a fabricated zero would silently weaken the border condition.
        continue;
      }
      std::vector<double> x_src
          = value_function_.get(source_phase, r, theta_idx, theta_end_idx);
      if (x_src.size() != static_cast<std::size_t>(source_layout.num_x)) {
        throw std::runtime_error(
            "getBorderConditions: source vector has invalid size");
      }
      candidates.push_back(std::move(x_src));
    }
  }
  if (candidates.empty()) {
    throw std::runtime_error(
        "getBorderConditions: no source candidates for the border condition");
  }
  const int n_candidates = static_cast<int>(candidates.size());

  CourierBorderRequest request;
  request.target_phase = target_phase;
  request.source_phase = source_phase;
  request.candidates = std::span<const std::vector<double>>(candidates);
  request.theta_idx = theta_idx;
  request.switch_cnt = switch_cnt;

  CourierBorderStats stats;
  std::vector<double> x_boundary
      = courier_solver_.solve(request, approximation_mode_, &stats);

  logger_->info(
      "Courier border LP target_phase={} source_phase={} theta_idx={} "
      "switch_cnt={} mode={} candidates={} iters={} cuts={} subproblems={} "
      "worst_zeta={:.3e} obj={:.6f} hit_box={}",
      target_phase, source_phase, theta_idx, switch_cnt,
      approximation_mode_ == ApproximationMode::Upper ? "upper" : "lower",
      n_candidates, stats.iterations, stats.cuts_added,
      stats.subproblems_solved, stats.worst_zeta, stats.master_objective,
      stats.master_hit_box);

  // Opt-in independent re-derivation: rebuilds every region's courier from
  // scratch rather than trusting the loop that produced x_boundary.
  if (validate_) {
    const double residual = courier_solver_.worstCertificateResidual(
        request, approximation_mode_, x_boundary);
    if (residual > kResidualValidationTol) {
      throw std::runtime_error(
          "getBorderConditions: independent re-verification failed, worst "
          "certificate residual " + std::to_string(residual));
    }
  }

  return x_boundary;
}

std::tuple<std::unique_ptr<Highs>, std::vector<double>, std::vector<double>>
BarycentricAffineApproximator::initializeHighs(int phase) {
  auto [starts, col_index, value, row_lower, row_upper, c_vec]
      = prepareLpMatrices(phase);
  const auto& layout = layouts_[phase];
  const int m = static_cast<int>(row_upper.size());
  const int n = static_cast<int>(c_vec.size());

  std::unique_ptr<Highs> highs = std::make_unique<Highs>();
  highs->setOptionValue("solver", "simplex");
  highs->setOptionValue("presolve", "on");
  highs->setOptionValue("simplex_strategy", 2);
  highs->setOptionValue("pdlp_optimality_tolerance", kHighsPdlpOptimalityTol);
  highs->setOptionValue("kkt_tolerance", kHighsSolutionTol);
  highs->setOptionValue("primal_feasibility_tolerance", kHighsSolutionTol);
  highs->setOptionValue("dual_feasibility_tolerance", kHighsSolutionTol);
  highs->setOptionValue("primal_residual_tolerance", kHighsSolutionTol);
  highs->setOptionValue("dual_residual_tolerance", kHighsSolutionTol);
  highs->setOptionValue("optimality_tolerance", kHighsSolutionTol);
  highs->setOptionValue("small_matrix_value", kHighsSmallMatrixValue);
  highs->setOptionValue("log_to_console", highs_verbose_);

  // The objective is the l1 norm of the residuals measured at the endpoint with
  // the known value, and it is minimized in both directions: the sign s is
  // already baked into c_vec by prepareLpMatrices (step 2.1, section 5).
  highs->changeObjectiveSense(ObjSense::kMinimize);

  const double inf = highs->getInfinity();
  std::vector<double> col_lower(n, -inf);
  std::vector<double> col_upper(n, inf);

  // y_j represents |Psi_j x| and must be nonnegative. The two abs rows enforce
  // y_j >= +/-Psi_j x; this lower bound supplies y_j >= 0 without extra rows.
  for (int j = 0; j < layout.num_regions; ++j) {
    for (int r = 0; r < kSpaceDim; ++r) {
      col_lower[layout.idxY(j, r)] = 0.0;
    }
  }

  // Gauge fixing removes the additive nullspace between the five projected
  // layers. Do not fix layer 0; layers 1..4 get their first vertex pinned to 0.
  for (int s = 1; s < kSubsystemCount; ++s) {
    if (layout.eta_s[s] == 0) {
      throw std::runtime_error("initializeHighs: empty projection layer");
    }
    const int col = layout.idxX(s, 0);
    col_lower[col] = 0.0;
    col_upper[col] = 0.0;
  }

  HighsStatus st = highs->addCols(
      n, c_vec.data(), col_lower.data(), col_upper.data(), /*num_nz=*/0,
      /*start=*/nullptr, /*index=*/nullptr, /*value=*/nullptr);
  if (st != HighsStatus::kOk) {
    throw std::runtime_error("initializeHighs: highs.addCols failed.");
  }

  st = highs->addRows(m, row_lower.data(), row_upper.data(),
                      static_cast<int>(value.size()), starts.data(),
                      col_index.data(), value.data());
  if (st != HighsStatus::kOk) {
    throw std::runtime_error("initializeHighs: highs.addRows failed.");
  }

  // Zero warm start is compatible with gauge-fixed x columns and y >= 0. The
  // solver may still need to find feasibility, but the initial values obey all
  // simple column bounds.
  HighsSolution initial_solution;
  initial_solution.value_valid = true;
  initial_solution.col_value.assign(n, 0.0);
  st = highs->setSolution(initial_solution);
  if (st != HighsStatus::kOk) {
    throw std::runtime_error("initializeHighs: highs.setSolution failed.");
  }

  return std::make_tuple(std::move(highs), std::move(row_lower),
                         std::move(row_upper));
}

void BarycentricAffineApproximator::updateHighsRhsUpperBounds(
    int phase, int solver_index, const std::vector<double>& x_next) {
  // Only the residual row RHS depends on x_next. Absolute-value rows and
  // gauge-fixed column bounds remain static.
  //   Left  row:  upper = -s g(nu)  -  s phi^T x_next / dt   (base + increment)
  //   Right row:  upper = -s (phi^T x_next + dt beta),
  //               beta  = q^T m + s rho^T |q| + g(nu),  q = Psi_j x_next.
  // The modulus lives here and only here: q is built from the known x_next, so
  // it is a vector of numbers and the LP stays linear (step 2.1, section 3.1).
  const auto& layout = layouts_[phase];
  if (x_next.size() != static_cast<std::size_t>(layout.num_x)) {
    throw std::invalid_argument("updateHighsRhsUpperBounds: x_next size must be "
                                + std::to_string(layout.num_x));
  }

  const double s = signS();
  const auto& row_lower = row_lowers_[solver_index];
  const auto& base_upper = row_uppers_[solver_index];
  auto& highs_solver = highs_solvers_[solver_index];
  std::vector<double> new_row_upper = base_upper;

  // q_j = Psi_j x_next, computed once per region and shared by all Right rows
  // of that region.
  const auto& psi_by_region = psi_by_region_[phase];
  std::vector<Eigen::VectorXd> q(psi_by_region.size());
  for (std::size_t j = 0; j < psi_by_region.size(); ++j) {
    q[j] = Eigen::VectorXd::Zero(kSpaceDim);
    for (int r = 0; r < kSpaceDim; ++r) {
      q[j](r) = psi_by_region[j].rows[r].dot(x_next);
    }
  }

  for (const auto& term : rhs_terms_[phase]) {
    const double phi_x = term.phi.dot(x_next);
    double upper = 0.0;
    if (term.kind == RhsKind::Left) {
      upper = base_upper[term.row_id] - s * phi_x / t_delta_;
    } else {
      if (term.region < 0 || term.region >= static_cast<int>(q.size())) {
        throw std::runtime_error(
            "updateHighsRhsUpperBounds: Right term has an invalid region id");
      }
      const Eigen::VectorXd& qj = q[term.region];
      const double beta
          = qj.dot(term.m) + s * term.rho.dot(qj.cwiseAbs()) + term.g;
      upper = -s * (phi_x + t_delta_ * beta);
    }
    new_row_upper[term.row_id] = (std::abs(upper) <= kEps) ? 0.0 : upper;
  }

  std::vector<int> row_ids(new_row_upper.size());
  std::iota(row_ids.begin(), row_ids.end(), 0);
  HighsStatus st = highs_solver->changeRowsBounds(
      static_cast<int>(row_ids.size()), row_ids.data(), row_lower.data(),
      new_row_upper.data());
  if (st != HighsStatus::kOk) {
    throw std::runtime_error(
        "updateHighsRhsUpperBounds: highs.changeRowsBounds failed.");
  }
}

std::vector<double> BarycentricAffineApproximator::solveLp(
    int phase, int solver_index) {
  // Solves the current reused HiGHS model and extracts only x_k. The auxiliary
  // y_j block is a proof/linearization device and is not part of the value
  // function stored for later time steps.
  auto& highs_solver = highs_solvers_[solver_index];
  HighsStatus run_status = highs_solver->run();
  if (run_status != HighsStatus::kOk) {
    throw std::runtime_error("solveLp: highs_solver.run() failed with status "
                             + std::to_string(static_cast<int>(run_status)));
  }
  if (highs_solver->getModelStatus() != HighsModelStatus::kOptimal) {
    throw std::runtime_error(
        "solveLp: LP solution not found for solver_index "
        + std::to_string(solver_index) + ", model status: "
        + std::to_string(static_cast<int>(highs_solver->getModelStatus())));
  }

  const auto& solution = highs_solver->getSolution();
  const int num_x = layouts_[phase].num_x;
  if (solution.col_value.size() < static_cast<std::size_t>(num_x)) {
    throw std::runtime_error("solveLp: HiGHS solution is too short.");
  }
  return std::vector<double>(solution.col_value.begin(),
                             solution.col_value.begin() + num_x);
}

void BarycentricAffineApproximator::validateStepResiduals(
    int phase, const std::vector<double>& x_next,
    const std::vector<double>& z) const {
  // Recomputes both residuals of the segment straight from the formulas and
  // checks their sign. This is the end-to-end check of the assembly: if the row
  // signs, the RHS updates or the psi/phi tables are wrong, it fires here rather
  // than silently producing a function that is not a bound.
  //
  //   F^L = phi^T d + (Psi_j z)^T m       + s rho^T |Psi_j z|       + g(nu)
  //   F^R = phi^T d + (Psi_j x_next)^T m  + s rho^T |Psi_j x_next|  + g(nu)
  //   d   = (x_next - z) / dt
  // Both must satisfy s * F <= 0, i.e. F >= 0 for the lower bound (s = -1) and
  // F <= 0 for the upper one (s = +1).
  const double s = signS();
  const auto& geometry = phase_geometries_[phase];
  const auto& psi_by_region = psi_by_region_[phase];
  const auto& A_j_matrs = A_j_matrs_[phase];
  const auto& f_j_vecs = f_j_vecs_[phase];
  const auto& Q_c_j_matrs = Q_c_j_matrs_[phase];
  const auto& q_c_j_vecs = q_c_j_vecs_[phase];
  const auto& Q_r_j_matrs = Q_r_j_matrs_[phase];
  const auto& q_r_j_vecs = q_r_j_vecs_[phase];
  const auto& g_j_vecs = g_j_vecs_[phase];
  const auto& g_j_scals = g_j_scals_[phase];

  const int n_regions = static_cast<int>(geometry.region_vertices.size());
  double worst = 0.0;
  for (int j = 0; j < n_regions; ++j) {
    Eigen::VectorXd q_left = Eigen::VectorXd::Zero(kSpaceDim);
    Eigen::VectorXd q_right = Eigen::VectorXd::Zero(kSpaceDim);
    for (int r = 0; r < kSpaceDim; ++r) {
      q_left(r) = psi_by_region[j].rows[r].dot(z);
      q_right(r) = psi_by_region[j].rows[r].dot(x_next);
    }

    for (const Eigen::VectorXd& nu : geometry.region_vertices[j]) {
      const SparseVec phi = buildPhiRow(phase, j, nu, kEps);
      const Eigen::VectorXd m = A_j_matrs[j] * nu + f_j_vecs[j]
                                + Q_c_j_matrs[j] * nu + q_c_j_vecs[j];
      Eigen::VectorXd rho = Q_r_j_matrs[j] * nu + q_r_j_vecs[j];
      rho = rho.cwiseMax(0.0);
      const double g_nu = g_j_vecs[j].dot(nu) + g_j_scals[j];

      const double slope = (phi.dot(x_next) - phi.dot(z)) / t_delta_;
      const double f_left = slope + q_left.dot(m)
                            + s * rho.dot(q_left.cwiseAbs()) + g_nu;
      const double f_right = slope + q_right.dot(m)
                             + s * rho.dot(q_right.cwiseAbs()) + g_nu;

      // s * F <= 0 is the requirement; anything above tolerance is a violation.
      worst = std::max(worst, s * f_left);
      worst = std::max(worst, s * f_right);
    }
  }

  if (worst > kResidualValidationTol) {
    throw std::runtime_error(
        "validateStepResiduals: residual has the wrong sign, worst s*F = "
        + std::to_string(worst));
  }
}

std::vector<double> BarycentricAffineApproximator::solveMainLpStep(
    int phase, int solver_index, const std::vector<double>& x_next) {
  updateHighsRhsUpperBounds(phase, solver_index, x_next);
  std::vector<double> z = solveLp(phase, solver_index);
  if (validate_) {
    validateStepResiduals(phase, x_next, z);
  }
  return z;
}

void BarycentricAffineApproximator::precomputeMatrices() {
  logger_->info("Starting barycentric precomputeMatrices");
  for (int phase = 0; phase < kPhases; ++phase) {
    if (phase_geometries_[phase].region_vertices.empty()) {
      throw std::runtime_error(
          "precomputeMatrices: geometry not initialized. Call "
          "getIntersectionPoints() first.");
    }
  }
  A_j_matrs_.clear();
  f_j_vecs_.clear();
  Q_c_j_matrs_.clear();
  q_c_j_vecs_.clear();
  Q_r_j_matrs_.clear();
  q_r_j_vecs_.clear();
  g_j_vecs_.clear();
  g_j_scals_.clear();

  for (int phase = 0; phase < kPhases; ++phase) {
    auto [A_j, f_j, Qc_j, qc_j, Qr_j, qr_j, g_j, g0_j]
        = precomputeSystemMatrices(phase);
    const int n_regions = static_cast<int>(A_j.size());
    std::size_t n_vertices = 0;
    for (const auto& verts : phase_geometries_[phase].region_vertices) {
      n_vertices += verts.size();
    }
    if (n_regions == 0) {
      logger_->info(
          "precomputeMatrices: phase={}, regions=0, vertices=0 (no system "
          "matrices)",
          phase);
    } else {
      logger_->info(
          "precomputeMatrices: phase={}, regions={}, vertices={}, A_j={}x{}, "
          "f_j={}, Q_c={}x{}, q_c={}, Q_r={}x{}, q_r={}, g_j={}, layout: "
          "num_x={}, num_y={}, num_cols={}",
          phase, n_regions, n_vertices, A_j[0].rows(), A_j[0].cols(),
          f_j[0].size(), Qc_j[0].rows(), Qc_j[0].cols(), qc_j[0].size(),
          Qr_j[0].rows(), Qr_j[0].cols(), qr_j[0].size(), g_j[0].size(),
          layouts_[phase].num_x, kSpaceDim * layouts_[phase].num_regions,
          layouts_[phase].num_cols);
    }
    A_j_matrs_.push_back(std::move(A_j));
    f_j_vecs_.push_back(std::move(f_j));
    Q_c_j_matrs_.push_back(std::move(Qc_j));
    q_c_j_vecs_.push_back(std::move(qc_j));
    Q_r_j_matrs_.push_back(std::move(Qr_j));
    q_r_j_vecs_.push_back(std::move(qr_j));
    g_j_vecs_.push_back(std::move(g_j));
    g_j_scals_.push_back(std::move(g0_j));
  }
  logger_->info("Finished barycentric precomputeMatrices");
}

void BarycentricAffineApproximator::run(const std::string& output_folder_path,
                                        int n_threads) {
  // Follows the global affine structure: geometry once, one HiGHS model per
  // (phase, worker) reused across the whole backward march, and only row bounds
  // updated inside it. Both the main LP and the border LP are implemented.
  const std::filesystem::path out_path(output_folder_path);
  if (!std::filesystem::exists(out_path)
      || !std::filesystem::is_directory(out_path)) {
    throw std::runtime_error(
        "output_folder_path does not exist or is not a directory: "
        + output_folder_path);
  }
  if (n_threads < 2 || n_threads % 2 != 0) {
    throw std::runtime_error(
        "n_threads must be at least 2 and even. Provided n_threads: "
        + std::to_string(n_threads));
  }

  logger_->info("Starting barycentric affine approximator");
  getIntersectionPoints();
  logger_->info(
      "Finished getIntersectionPoints. Got {} areas for phase 0 and {} areas "
      "for phase 1",
      phase_geometries_[0].region_vertices.size(),
      phase_geometries_[1].region_vertices.size());
  courier_options_.highs_verbose = highs_verbose_;
  courier_solver_ = CourierBorderSolver(courier_options_);
  courier_solver_.prepare(phase_geometries_, layouts_, node_weights_,
                          system_params_.N);
  logger_->info("Finished preparing the courier border solver");

  precomputeMatrices();
  logger_->info(
      "Finished precomputeMatrices. Precomputed {} system matrices for phase 0 "
      "and {} system matrices for phase 1",
      A_j_matrs_[0].size(), A_j_matrs_[1].size());

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

  // Initialize solver instances once and then reuse them. This is intentionally
  // serialized: prepareLpMatrices() also refreshes shared per-phase RHS metadata
  // (rhs_terms_[phase]), so parallel initialization of several solvers for the
  // same phase would race on that metadata. The runtime optimization that
  // matters for the algorithm is still preserved: the expensive HiGHS models are
  // built once and only row bounds are updated inside the backward loop.
  for (int phase = 0; phase < kPhases; ++phase) {
    for (int s = 0; s < solvers_per_phase; ++s) {
      const int solver_index = phase * solvers_per_phase + s;
      auto [highs_solver, row_lower, row_upper] = initializeHighs(phase);
      highs_solvers_[solver_index] = std::move(highs_solver);
      row_lowers_[solver_index] = std::move(row_lower);
      row_uppers_[solver_index] = std::move(row_upper);
    }
  }
  logger_->info("Done initializing Highs solvers");

  ThreadPool pool(n_threads);
  for (int switch_cnt = 0; switch_cnt <= max_switches_; ++switch_cnt) {
    logger_->info("Computing value function for switch count {}/{}", switch_cnt,
                  max_switches_);
    auto theta_range_ids
        = theta_t_index_lists_.expanded_t_by_k_theta[switch_cnt];
    const std::size_t total_tasks
        = theta_range_ids.size() * static_cast<std::size_t>(kPhases);
    std::atomic<std::size_t> completed_tasks{0};
    std::vector<std::future<void>> futures;
    futures.reserve(theta_range_ids.size() * kPhases);

    for (auto [theta_idx, t_range_ids] : theta_range_ids) {
      for (int phase = 0; phase < kPhases; ++phase) {
        const int solver_index
            = phase * solvers_per_phase
              + (static_cast<int>(theta_idx) % solvers_per_phase);
        futures.push_back(pool.enqueue([this, switch_cnt, phase, theta_idx,
                                        t_range_ids, solver_index, total_tasks,
                                        &completed_tasks]() {
          const double theta = t_range_[theta_idx];
          logger_->info(
              "Starting computation for theta_idx: {}, phase: {}, switch_cnt: "
              "{}",
              theta_idx, phase, switch_cnt);
          std::lock_guard<std::mutex> lock(*solver_mutexes_[solver_index]);
          std::vector<double> x_next;

          const auto n_t = static_cast<std::ptrdiff_t>(t_range_ids.size());
          if (n_t == 0) {
            throw std::runtime_error("run: empty t_range_ids");
          }

          for (std::ptrdiff_t i_t_idx = n_t - 1; i_t_idx >= 0; --i_t_idx) {
            const int t_idx = t_range_ids[i_t_idx];
            if (i_t_idx == n_t - 1) {
              if (t_idx != theta_idx) {
                throw std::runtime_error("run: last time id is not theta_idx");
              }
              const int switch_phase = phase == 0 ? 1 : 0;
              x_next = getBorderConditions(switch_phase, theta_idx, theta,
                                           switch_cnt);
            } else {
              x_next = solveMainLpStep(phase, solver_index, x_next);
            }
            value_function_.set(phase, switch_cnt, t_idx, theta_idx, x_next,
                                static_cast<std::size_t>(
                                    layouts_[phase].num_x));
            auto [min_it, max_it]
                = std::minmax_element(x_next.begin(), x_next.end());
            const double min_val = min_it == x_next.end() ? 0.0 : *min_it;
            const double max_val = max_it == x_next.end() ? 0.0 : *max_it;
            logger_->info(
                "Value function min/max at t_idx: {}, theta_idx: {}, phase: "
                "{}, switch_cnt: {}: \t{:.4f}\t{:.4f}",
                t_idx, theta_idx, phase, switch_cnt, min_val, max_val);
          }
          const std::size_t done = completed_tasks.fetch_add(1) + 1;
          logger_->info(
              "Finished computation for theta_idx: {}, phase: {}, switch_cnt: "
              "{} ({}/{})",
              theta_idx, phase, switch_cnt, done, total_tasks);
        }));
      }
    }
    for (auto& f : futures) {
      f.get();
    }
    logger_->info("Completed switch count {}/{}", switch_cnt, max_switches_);
  }

  const std::filesystem::path base(output_folder_path);
  logger_->info("Saving barycentric value function artifacts to {}", base.string());
  value_function_.dumpToJson((base / "value_function.json").string());
  hcpwa::util::dumpVectorToJson(t_range_, (base / "t_range.json").string());
  dumpInitParamsToJson((base / "init_params.json").string());
  logger_->info("Finished barycentric affine approximator run");
}

}  // namespace barycentric_affine_approximator

// NOLINTEND(readability-identifier-naming)
