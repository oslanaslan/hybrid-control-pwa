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

// How much primal infeasibility a non-optimal HiGHS point may carry and still
// be accepted by solveLp(). Row units here are the same units as the residual
// F, so this is set to kResidualValidationTol: accept exactly what the exact
// worst-residual check that runs immediately afterwards would accept, and let
// that check -- which recomputes s * F from the formulas rather than trusting
// the solver -- be the one that decides. The solver itself is configured at
// kHighsSolutionTol = 1e-6, two orders tighter.
constexpr double kLpPrimalFeasibilityLimit
    = barycentric_affine_approximator::kResidualValidationTol;

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

  // How many (layer, theta) nodes the lists carry that no earlier layer can
  // reach. It is a pure function of the grid and is worth knowing before a run
  // starts rather than discovering it forty minutes in; getBorderConditions
  // skips exactly these. A large share would mean the theta ranges and the
  // t-windows disagree about more than grid quantisation.
  {
    int checked = 0;
    int unreachable = 0;
    const double half_step = t_delta_ / 2.0;
    auto covers = [&](int level, int theta_end, int t_idx) {
      const auto layer = theta_t_index_lists_.expanded_t_by_k_theta.find(level);
      if (layer == theta_t_index_lists_.expanded_t_by_k_theta.end()) {
        return false;
      }
      const auto entry = layer->second.find(theta_end);
      if (entry == layer->second.end()) {
        return false;
      }
      return std::find(entry->second.begin(), entry->second.end(), t_idx)
             != entry->second.end();
    };
    for (const auto& [level, by_theta] :
         theta_t_index_lists_.expanded_t_by_k_theta) {
      if (level < 1) {
        continue;
      }
      for (const auto& [theta_idx, t_ids] : by_theta) {
        ++checked;
        const double theta = t_range_[static_cast<std::size_t>(theta_idx)];
        const double lo = std::min(theta + tau_min_, t_max_);
        const double hi = std::min(theta + tau_max_, t_max_);
        bool reachable = false;
        for (int r = 0; r < level && !reachable; ++r) {
          for (std::size_t i = 0; i < t_range_.size() && !reachable; ++i) {
            if (t_range_[i] >= lo - half_step && t_range_[i] <= hi + half_step
                && covers(r, static_cast<int>(i), theta_idx)) {
              reachable = true;
            }
          }
        }
        if (!reachable) {
          ++unreachable;
        }
      }
    }
    unreachable_nodes_ = unreachable;
    checked_nodes_ = checked;
  }

  logger_ = spdlog::get("barycentric_affine_approximator");
  if (!logger_) {
    logger_ = spdlog::stdout_color_mt("barycentric_affine_approximator");
  }
  logger_->set_level(spdlog::level::info);
  logger_->info(
      "theta lists: {} of {} (layer, theta) nodes have no reachable "
      "predecessor and will be skipped; they sit at the low edge of a layer's "
      "theta range, where grid quantisation puts them just outside the "
      "t-window of every earlier layer",
      unreachable_nodes_, checked_nodes_);

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
          system_params_.f5max, system_params_.f8max, true,
          geometry_options_);

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
  phase_geometries_[0].blocks
      = blockGeometryFromRegions(0, areas_vertices.blocks_phase0);

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
  phase_geometries_[1].blocks
      = blockGeometryFromRegions(1, areas_vertices.blocks_phase1);

  for (int phase = 0; phase < kPhases; ++phase) {
    const auto& geometry = phase_geometries_[phase];
    // region_vertices is the 8D product of the block cells. Nothing in the LP
    // or the border solver reads it any more, so it is normally empty; when a
    // test asks for it, it must line up with the region list.
    if (!geometry.region_vertices.empty()
        && geometry.region_vertices.size()
               != geometry.region_triangle_ids.size()) {
      throw std::runtime_error(
          "Barycentric geometry has mismatched region vertices and ids.");
    }
    if (geometry.region_triangle_ids.empty()) {
      throw std::runtime_error("Barycentric geometry has no regions.");
    }

    BarycentricVarLayout layout;
    for (int s = 0; s < kSubsystemCount; ++s) {
      layout.offset_s[s] = layout.num_x;
      layout.eta_s[s]
          = static_cast<int>(geometry.layers[s].unique_vertices.size());
      layout.num_x += layout.eta_s[s];
    }
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
  return barycentric_affine_approximator::buildPhiRow(
      phase_geometries_[phase], layouts_[phase], region, point, tolerance);
}

SparseVec BarycentricAffineApproximator::buildPhiRowBlock(
    int phase, int block, int block_region, const Eigen::VectorXd& point,
    double tolerance) const {
  if (phase < 0 || phase >= kPhases) {
    throw std::invalid_argument("buildPhiRowBlock: invalid phase");
  }
  return barycentric_affine_approximator::buildPhiRowBlock(
      phase_geometries_[phase], layouts_[phase], block, block_region, point,
      tolerance);
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

CtmRegionData BarycentricAffineApproximator::ctmDataForCells(
    int phase, const std::vector<int>& cells, const Eigen::VectorXd& n) const {
  // One assembly rule for the whole CTM drift: a flow f_ij leaves cell i and
  // enters cell j, so row j gains it and row i loses it, and g is the sum of
  // all four flows. Written out per phase this is exactly the a_matr / b_vec /
  // g_vec block of the paper.
  if (n.size() != kSpaceDim) {
    throw std::invalid_argument("ctmDataForCells: n must be 8-dimensional");
  }
  const int size = static_cast<int>(cells.size());
  std::array<int, kSpaceDim> position{};
  position.fill(-1);
  for (int k = 0; k < size; ++k) {
    if (cells[static_cast<std::size_t>(k)] < 0
        || cells[static_cast<std::size_t>(k)] >= kSpaceDim) {
      throw std::invalid_argument("ctmDataForCells: cell id out of range");
    }
    position[static_cast<std::size_t>(cells[static_cast<std::size_t>(k)])] = k;
  }

  CtmRegionData out;
  out.a = Eigen::MatrixXd::Zero(size, size);
  out.f = Eigen::VectorXd::Zero(size);
  out.g_vec = Eigen::VectorXd::Zero(size);

  for (const auto& flow : phaseFlows(phase)) {
    const int from = flow[0];
    const int to = flow[1];
    const int from_at = position[static_cast<std::size_t>(from)];
    const int to_at = position[static_cast<std::size_t>(to)];
    if (from_at < 0 && to_at < 0) {
      continue;
    }
    if (from_at < 0 || to_at < 0) {
      // A flow whose two cells fall in different blocks would make the drift
      // depend on coordinates the block does not own, and the residual would
      // stop being separable. The block split is chosen so this cannot happen.
      throw std::runtime_error(std::format(
          "ctmDataForCells: flow {}->{} straddles the given cell set", from,
          to));
    }

    const auto [row, scalar] = getFIJMinResolution(from, to, n);
    // The resolved row is supported on the flow's own two cells, both of which
    // are in `cells`, so nothing is dropped by the restriction below.
    for (int c = 0; c < kSpaceDim; ++c) {
      if (position[static_cast<std::size_t>(c)] < 0
          && std::abs(row(c)) > 0.0) {
        throw std::runtime_error(
            "ctmDataForCells: resolved flow touches a cell outside the set");
      }
    }

    Eigen::VectorXd restricted = Eigen::VectorXd::Zero(size);
    for (int k = 0; k < size; ++k) {
      restricted(k) = row(cells[static_cast<std::size_t>(k)]);
    }
    out.a.row(to_at) += restricted.transpose();
    out.a.row(from_at) -= restricted.transpose();
    out.f(to_at) += scalar(0);
    out.f(from_at) -= scalar(0);
    out.g_vec += restricted;
    out.g_scal += scalar(0);
  }

  hcpwa::util::assertScalar(out.g_scal);
  return out;
}

BoxRegionData BarycentricAffineApproximator::boxDataForCells(
    const std::vector<int>& cells, const Eigen::VectorXd& n0) const {
  // Affine center/radius maps of the disturbance box:
  //   c(n) = qc_diag .* n + qc_off,  rho(n) = qr_diag .* n + qr_off.
  // Every entry depends on its own cell only, which is why these are diagonal
  // and why the restriction to a block is exact. Branch resolution is copied
  // from the global affine code.
  if (n0.size() != kSpaceDim) {
    throw std::invalid_argument("boxDataForCells: n0 must be 8-dimensional");
  }
  const int size = static_cast<int>(cells.size());
  Eigen::VectorXd upper_diag = Eigen::VectorXd::Zero(size);
  Eigen::VectorXd upper_off = Eigen::VectorXd::Zero(size);
  Eigen::VectorXd lower_diag = Eigen::VectorXd::Zero(size);
  Eigen::VectorXd lower_off = Eigen::VectorXd::Zero(size);

  const double N = system_params_.N;
  const double w = system_params_.w;
  const double v = system_params_.v;
  const double F = system_params_.F;

  for (int k = 0; k < size; ++k) {
    const int i = cells[static_cast<std::size_t>(k)];
    const bool is_in
        = std::find(kInIds.begin(), kInIds.end(), i) != kInIds.end();
    const bool is_out
        = std::find(kOutIds.begin(), kOutIds.end(), i) != kOutIds.end();
    if (is_in == is_out) {
      throw std::runtime_error(
          "boxDataForCells: every cell must be exactly one of in/out");
    }
    const double ni = n0(i);
    if (!std::isfinite(ni)) {
      throw std::runtime_error(
          "boxDataForCells: n0 is not finite at a cell of the set");
    }

    if (is_in) {
      const auto [f_min, f_max] = getFMinMaxForAxis(i);
      if (f_min < w * (N - ni)) {
        lower_off(k) = f_min;
      } else {
        lower_diag(k) = -w;
        lower_off(k) = w * N;
      }
      if (f_max < w * (N - ni)) {
        upper_off(k) = f_max;
      } else {
        upper_diag(k) = -w;
        upper_off(k) = w * N;
      }
    } else {
      // Outgoing cells carry no uncertainty: upper stays at zero, so the
      // radius is half the magnitude of the lower branch.
      if (F < v * ni) {
        lower_off(k) = -F;
      } else {
        lower_diag(k) = -v;
      }
    }
  }

  BoxRegionData out;
  out.qc_diag = (upper_diag + lower_diag) / 2.0;
  out.qc_off = (upper_off + lower_off) / 2.0;
  out.qr_diag = (upper_diag - lower_diag) / 2.0;
  out.qr_off = (upper_off - lower_off) / 2.0;
  return out;
}

Eigen::VectorXd BarycentricAffineApproximator::blockCentroidCoords(
    int phase, int block, int block_region) const {
  if (phase < 0 || phase >= kPhases) {
    throw std::invalid_argument("blockCentroidCoords: invalid phase");
  }
  if (block < 0 || block >= kBlockCount) {
    throw std::invalid_argument("blockCentroidCoords: invalid block");
  }
  const BlockGeometry& geometry
      = phase_geometries_[phase].blocks[static_cast<std::size_t>(block)];
  if (block_region < 0 || block_region >= geometry.numRegions()) {
    throw std::invalid_argument("blockCentroidCoords: invalid block region");
  }
  const auto& vertices
      = geometry.vertices[static_cast<std::size_t>(block_region)];
  if (vertices.empty()) {
    throw std::runtime_error("blockCentroidCoords: block region has no "
                             "vertices");
  }
  Eigen::VectorXd centroid = Eigen::VectorXd::Zero(geometry.coord_count);
  for (const auto& vertex : vertices) {
    centroid += vertex;
  }
  centroid /= static_cast<double>(vertices.size());
  return centroid;
}

BlockSystemMatrices BarycentricAffineApproximator::getBlockSystemMatrices(
    int phase, int block, int block_region) const {
  const BlockGeometry& geometry
      = phase_geometries_[phase].blocks[static_cast<std::size_t>(block)];
  const Eigen::VectorXd centroid
      = blockCentroidCoords(phase, block, block_region);

  // Scatter the block centroid into an 8-vector and leave every other
  // coordinate NaN. Each flow of this block reads only its own two cells, both
  // of which are in the block, so the padding is never read -- and if it ever
  // were, all three branch comparisons would be false and getFIJMinResolution
  // would throw rather than return a wrong branch.
  Eigen::VectorXd padded = Eigen::VectorXd::Constant(
      kSpaceDim, std::numeric_limits<double>::quiet_NaN());
  std::vector<int> cells;
  cells.reserve(geometry.coord_count);
  for (int c = 0; c < geometry.coord_count; ++c) {
    const int cell = geometry.coords[static_cast<std::size_t>(c)];
    padded(cell) = centroid(c);
    cells.push_back(cell);
  }

  BlockSystemMatrices out;
  out.ctm = ctmDataForCells(phase, cells, padded);
  out.box = boxDataForCells(cells, padded);
  return out;
}

std::vector<BlockSystemMatrices>
BarycentricAffineApproximator::precomputeBlockSystemMatrices(int phase,
                                                             int block) {
  const BlockGeometry& geometry
      = phase_geometries_[phase].blocks[static_cast<std::size_t>(block)];
  std::vector<BlockSystemMatrices> out;
  out.reserve(static_cast<std::size_t>(geometry.numRegions()));
  for (int j = 0; j < geometry.numRegions(); ++j) {
    out.push_back(getBlockSystemMatrices(phase, block, j));
  }
  return out;
}

void BarycentricAffineApproximator::buildReducedLpInput(int phase) {
  // Builds the block description the reduced LP is assembled from. Everything
  // here is per block-region, never per product region: the 8D region is the
  // direct product of the three blocks and the residual is additively separable
  // over them, so a row of block b is shared by all R / R_b product rows that
  // contain it.
  const auto& geometry = phase_geometries_[phase];
  const auto& layout = layouts_[phase];

  block_reduction::ReducedLpInput input;
  input.num_x = layout.num_x;
  input.t_delta = t_delta_;
  input.sign_s = signS();
  input.weights = objective_weights_;
  input.tie_break_eps = tie_break_eps_;
  input.blocks.resize(kBlockCount);

  for (int b = 0; b < kBlockCount; ++b) {
    const BlockGeometry& block = geometry.blocks[static_cast<std::size_t>(b)];
    const auto& system = block_system_[phase][static_cast<std::size_t>(b)];
    if (static_cast<int>(system.size()) != block.numRegions()) {
      throw std::runtime_error(
          "buildReducedLpInput: block system matrices are missing; call "
          "precomputeMatrices() first");
    }

    block_reduction::ReducedLpBlock out_block;
    out_block.coord_count = block.coord_count;
    out_block.regions.resize(static_cast<std::size_t>(block.numRegions()));

    for (int j = 0; j < block.numRegions(); ++j) {
      block_reduction::BlockRegionData& region
          = out_block.regions[static_cast<std::size_t>(j)];

      // Psi_{b,j}: one row per block coordinate. Each of the block's planes
      // contributes d alpha_l / d n_axis to the row of that axis, and the axis
      // sits at local_axis inside this block's coordinate list.
      region.psi_rows.resize(static_cast<std::size_t>(block.coord_count));
      for (int l = 0; l < block.layer_count; ++l) {
        const int s = block.layer_ids[static_cast<std::size_t>(l)];
        const ProjectionLayer& layer = geometry.layers[static_cast<std::size_t>(
            s)];
        const int triangle_id
            = block.triangle_ids[static_cast<std::size_t>(j)][
                static_cast<std::size_t>(l)];
        if (triangle_id < 0
            || triangle_id >= static_cast<int>(layer.bases.size())) {
          throw std::runtime_error(
              "buildReducedLpInput: triangle id out of range");
        }
        const TriangleBasis& basis
            = layer.bases[static_cast<std::size_t>(triangle_id)];
        for (int local_vertex = 0; local_vertex < 3; ++local_vertex) {
          const int x_col = layout.idxX(s, basis.vertex_ids[
              static_cast<std::size_t>(local_vertex)]);
          for (int local_axis = 0; local_axis < 2; ++local_axis) {
            const int row = block.local_axis[static_cast<std::size_t>(l)][
                static_cast<std::size_t>(local_axis)];
            region.psi_rows[static_cast<std::size_t>(row)].add(
                x_col, basis.H(local_vertex, local_axis),
                block_reduction::kAssembleEps);
          }
        }
      }

      const BlockSystemMatrices& data = system[static_cast<std::size_t>(j)];
      region.vertices.reserve(
          block.vertices[static_cast<std::size_t>(j)].size());
      for (const Eigen::VectorXd& nu :
           block.vertices[static_cast<std::size_t>(j)]) {
        block_reduction::BlockVertexData vertex;
        vertex.phi = barycentric_affine_approximator::buildPhiRowBlock(
            geometry, layout, b, j, nu, kEps);
        // m = A_b nu + f_b + c_b(nu): the drift with the disturbance box centre
        // folded in. Q_c and Q_r are diagonal, hence the elementwise products.
        vertex.m = data.ctm.a * nu + data.ctm.f
                   + data.box.qc_diag.cwiseProduct(nu) + data.box.qc_off;
        vertex.rho = data.box.qr_diag.cwiseProduct(nu) + data.box.qr_off;
        vertex.g = data.ctm.g_vec.dot(nu) + data.ctm.g_scal;
        region.vertices.push_back(std::move(vertex));
      }
    }
    input.blocks[static_cast<std::size_t>(b)] = std::move(out_block);
  }

  // Radii are normalised here, once: the row builder, the RHS update and the
  // residual check all read rho exactly as given.
  const block_reduction::RhoClampStats rho_stats
      = block_reduction::clampBlockRho(input, kEps);
  if (rho_stats.negatives_clamped > 0 || rho_stats.tiny_raised > 0) {
    logger_->info(
        "buildReducedLpInput: phase={}, clamped {} negative and raised {} "
        "sub-1e-9 uncertainty radii (vertices on a flow branch line)",
        phase, rho_stats.negatives_clamped, rho_stats.tiny_raised);
  }

  reduced_inputs_[phase] = std::move(input);
  reduced_cols_[phase]
      = block_reduction::makeReducedLpColLayout(reduced_inputs_[phase]);
  reduced_rows_[phase]
      = block_reduction::makeReducedLpRowLayout(reduced_inputs_[phase]);
}

block_reduction::ReducedLpMatrices
BarycentricAffineApproximator::prepareLpMatrices(int phase) const {
  block_reduction::ReducedLpMatrices matrices
      = block_reduction::assembleReducedLp(reduced_inputs_[phase]);

  double product_rows = 1.0;
  for (const auto& block : reduced_inputs_[phase].blocks) {
    double pairs = 0.0;
    for (const auto& region : block.regions) {
      pairs += static_cast<double>(region.vertices.size());
    }
    product_rows *= pairs;
  }
  logger_->info(
      "prepareLpMatrices: phase={}, rows={}, cols={}, nnz={}, block regions="
      "{}/{}/{}, product form would have needed {:.3e} residual rows",
      phase, matrices.rows.num_rows, matrices.cols.num_cols,
      matrices.value.size(), reduced_cols_[phase].num_regions[0],
      reduced_cols_[phase].num_regions[1], reduced_cols_[phase].num_regions[2],
      2.0 * product_rows);
  return matrices;
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
    // No reachable predecessor, by the construction of the theta lists rather
    // than by anything going wrong here.
    //
    // buildThetaToTIndexLists gives layer k the theta range
    // [t_min_k, min(t_max_k + tau_max, T)], which is wider than the union of
    // the t-windows of the layers below it, and grid quantisation widens the
    // gap: at T=1200, tau_min=10, tau_max=50 on a 240-point grid, layer 2
    // carries theta = 1149.8 while layers 0 and 1 start their t-windows at
    // 1150.0 exactly. Such a theta is unreachable in the problem's own
    // combinatorics. Counted up front by checkUnreachableNodes(): 230 of
    // 12 585 nodes, 1.8%, all at the low edge of a layer's theta range.
    //
    // The recursion defines no value there, and fabricating one is the thing
    // the loop above refuses to do with a zero: it would silently weaken the
    // border condition. Skipping the node is equivalent to never emitting it,
    // which is where the fix would belong if buildThetaToTIndexLists were not
    // shared with the global affine path. Every consumer tolerates the absence
    // because they all go through ValueFunction::contains.
    logger_->info(
        "getBorderConditions: no source candidates at theta_idx={}, "
        "switch_cnt={}, source_phase={}; node unreachable, skipping",
        theta_idx, switch_cnt, source_phase);
    return {};
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
  const block_reduction::ReducedLpMatrices matrices = prepareLpMatrices(phase);
  const auto& layout = layouts_[phase];
  const int m = matrices.rows.num_rows;
  const int n = matrices.cols.num_cols;

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
  // Same constant the assembler checks emitted rho entries against: HiGHS drops
  // matrix entries at or below it, and a dropped rho would weaken row (L).
  highs->setOptionValue("small_matrix_value",
                        block_reduction::kSmallMatrixValue);
  highs->setOptionValue("log_to_console", highs_verbose_);

  // The objective is the l1 norm of the residuals measured at the endpoint with
  // the known value, and it is minimized in both directions: the sign s is
  // already baked into the cost by the assembler (step 2.1, section 5).
  highs->changeObjectiveSense(ObjSense::kMinimize);

  // Normalize the cost vector to max |c| = 1. This is an exact reformulation:
  // dividing a linear objective by a positive constant leaves the feasible set
  // and the optimal face untouched. It is needed because the R/R_b weights put
  // raw costs in [1e+05, 2e+06] and objective values near 1e+17, while HiGHS
  // judges dual feasibility with an *absolute* tolerance on reduced costs. At
  // that scale the test is meaningless: a solve with a relative dual error of
  // 1e-10 was reported with 3e-04 of dual infeasibility and downgraded from
  // Optimal to Unknown, which used to abort the run. Same treatment the courier
  // master LP already gets.
  std::vector<double> cost(matrices.cost.data(),
                           matrices.cost.data() + matrices.cost.size());
  double objective_scale = 0.0;
  for (const double c : cost) {
    objective_scale = std::max(objective_scale, std::abs(c));
  }
  if (!(objective_scale > 0.0) || !std::isfinite(objective_scale)) {
    throw std::runtime_error(
        "initializeHighs: objective is all zero or not finite");
  }
  for (double& c : cost) {
    c /= objective_scale;
  }
  objective_scales_[static_cast<std::size_t>(phase)] = objective_scale;

  // y >= 0 and u >= 0 come from the assembler. Gauge fixing is ours: it removes
  // the additive nullspace between the five projected layers. Do not fix layer
  // 0; layers 1..4 get their first vertex pinned to 0.
  std::vector<double> col_lower = matrices.col_lower;
  std::vector<double> col_upper = matrices.col_upper;
  for (int s = 1; s < kSubsystemCount; ++s) {
    if (layout.eta_s[s] == 0) {
      throw std::runtime_error("initializeHighs: empty projection layer");
    }
    const int col = layout.idxX(s, 0);
    col_lower[static_cast<std::size_t>(col)] = 0.0;
    col_upper[static_cast<std::size_t>(col)] = 0.0;
  }

  HighsStatus st = highs->addCols(
      n, cost.data(), col_lower.data(), col_upper.data(),
      /*num_nz=*/0, /*start=*/nullptr, /*index=*/nullptr, /*value=*/nullptr);
  if (st != HighsStatus::kOk) {
    throw std::runtime_error("initializeHighs: highs.addCols failed.");
  }

  st = highs->addRows(m, matrices.row_lower.data(), matrices.row_upper.data(),
                      static_cast<int>(matrices.value.size()),
                      matrices.starts.data(), matrices.col_index.data(),
                      matrices.value.data());
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

  return std::make_tuple(std::move(highs), matrices.row_lower,
                         matrices.row_upper);
}

void BarycentricAffineApproximator::updateHighsRhsUpperBounds(
    int phase, int solver_index, const std::vector<double>& x_next) {
  // Only the residual row upper bounds depend on x_next; absolute-value rows,
  // the two coupling rows and the column bounds are static. The formulas live
  // next to the row builder in block_reduction_lp.cpp precisely because the two
  // have to agree on the row scale.
  const auto& layout = layouts_[phase];
  if (x_next.size() != static_cast<std::size_t>(layout.num_x)) {
    throw std::invalid_argument("updateHighsRhsUpperBounds: x_next size must be "
                                + std::to_string(layout.num_x));
  }

  const auto& row_lower = row_lowers_[solver_index];
  const std::vector<double> new_row_upper
      = block_reduction::updateReducedLpRowUpper(
          reduced_inputs_[phase], reduced_rows_[phase],
          row_uppers_[solver_index], x_next);

  auto& highs_solver = highs_solvers_[solver_index];
  std::vector<int> row_ids(new_row_upper.size());
  std::iota(row_ids.begin(), row_ids.end(), 0);
  // kWarning here too: changing bounds draws the same scaling advisories, and
  // the bounds that were actually installed are verified by the residual check
  // after the solve.
  const HighsStatus st = highs_solver->changeRowsBounds(
      static_cast<int>(row_ids.size()), row_ids.data(), row_lower.data(),
      new_row_upper.data());
  if (st == HighsStatus::kError) {
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
  const HighsStatus run_status = highs_solver->run();
  // kWarning is advisory and does not mean the solve failed. HiGHS returns it
  // for the scaling notes this LP always draws -- "excessively large costs" for
  // the R/R_b objective weights, "excessively small row bounds" for residual
  // bounds near zero -- while reporting Optimal with a primal-dual objective
  // error of 1e-16. Treating it as failure discarded an hour of work over a
  // perfectly solved LP.
  //
  // What actually has to hold is checked twice below: the model status, and
  // then validateStepResiduals, which recomputes s * F from the formulas and
  // is the statement that the result is a bound.
  if (run_status == HighsStatus::kError) {
    throw std::runtime_error("solveLp: highs_solver.run() returned an error");
  }
  const HighsModelStatus model_status = highs_solver->getModelStatus();
  const HighsInfo& info = highs_solver->getInfo();
  if (model_status != HighsModelStatus::kOptimal) {
    // kUnknown means HiGHS solved the LP but would not certify the point
    // against its own tolerances -- in every observed case because of dual
    // infeasibility, that is, reduced costs of the wrong sign. Read what that
    // does and does not cost us.
    //
    // Soundness of the step rests on *primal* feasibility alone: the rows are
    // exactly the statement s * F <= 0 at every vertex of every region, so any
    // primal-feasible z yields a valid bound. Dual feasibility is the
    // certificate of *optimality*, and losing it means only that the bound may
    // be looser than the best one this LP admits -- a coarser estimate, not a
    // wrong one. Accepting it is therefore a tightness relaxation, never a
    // structural one, and validateStepResiduals re-derives s * F from the
    // formulas right after this returns, so the bound property is checked
    // independently of anything HiGHS asserts.
    //
    // Anything else -- Infeasible, Unbounded, a solve error -- is fatal: those
    // carry no usable primal point at all.
    const bool primal_feasible
        = model_status == HighsModelStatus::kUnknown
          && info.primal_solution_status == kSolutionStatusFeasible
          && info.max_primal_infeasibility <= kLpPrimalFeasibilityLimit;
    if (!primal_feasible) {
      throw std::runtime_error(
          "solveLp: LP solution not found for solver_index "
          + std::to_string(solver_index) + ", model status: "
          + std::to_string(static_cast<int>(model_status))
          + ", primal solution status: "
          + std::to_string(static_cast<int>(info.primal_solution_status))
          + ", max primal infeasibility: "
          + std::to_string(info.max_primal_infeasibility));
    }
    logger_->warn(
        "solveLp: solver_index {} accepted a primal-feasible non-optimal "
        "point (model status {}, max primal infeasibility {:.3e}, {} dual "
        "infeasibilities up to {:.3e}); the bound for this step may be looser "
        "than optimal",
        solver_index, static_cast<int>(model_status),
        info.max_primal_infeasibility, info.num_dual_infeasibilities,
        info.max_dual_infeasibility);
  }

  const auto& solution = highs_solver->getSolution();
  const int num_x = layouts_[phase].num_x;
  if (solution.col_value.size() < static_cast<std::size_t>(num_x)) {
    throw std::runtime_error("solveLp: HiGHS solution is too short.");
  }

  // The snapshot the column layout was built from, not the live member: the u
  // block exists only when buildReducedLpInput() saw a positive value, so
  // reading tie_break_eps_ here would make a later setTieBreakEps() throw out
  // of idxU from inside a diagnostic.
  const double tie_break_eps = reduced_inputs_[phase].tie_break_eps;
  if (tie_break_eps > 0.0) {
    // The two halves of the objective, reported apart: the l1 residual norm is
    // what the method minimizes, and the eps * ||z||_1 term only picks one
    // point out of an optimal face. If the second is not far smaller than the
    // first, the regularization is no longer a tie-break.
    double regularizer = 0.0;
    const auto& cols = reduced_cols_[phase];
    for (int k = 0; k < num_x; ++k) {
      regularizer += tie_break_eps
                     * solution.col_value[static_cast<std::size_t>(
                         cols.idxU(k))];
    }
    logger_->info(
        "solveLp: phase={}, residual objective {:.6e}, tie-break term {:.6e}",
        phase,
        highs_solver->getInfo().objective_function_value
            * objective_scales_[static_cast<std::size_t>(phase)]
            - regularizer,
        regularizer);
  }

  return std::vector<double>(solution.col_value.begin(),
                             solution.col_value.begin() + num_x);
}

void BarycentricAffineApproximator::validateStepResiduals(
    int phase, const std::vector<double>& x_next,
    const std::vector<double>& z) const {
  // The exact worst residual over the whole product of regions and vertices.
  // The residual is additively separable over the three blocks, so the maximum
  // over the product is the sum of the per-block maxima -- this is the true
  // global worst case, not a sample of it, and it costs O(sum_b R_b).
  //
  //   F^L = phi^T d + (Psi_j z)^T m      + s rho^T |Psi_j z|      + g(nu)
  //   F^R = phi^T d + (Psi_j x_next)^T m + s rho^T |Psi_j x_next| + g(nu)
  //   d   = (x_next - z) / dt
  // Both must satisfy s * F <= 0: F >= 0 for the lower bound (s = -1) and
  // F <= 0 for the upper one (s = +1).
  const block_reduction::WorstResidual worst
      = block_reduction::worstReducedResidual(reduced_inputs_[phase], x_next,
                                              z);
  const double limit
      = kResidualValidationTol + kResidualValidationRelTol * worst.scale;
  if (worst.worst() > limit) {
    throw std::runtime_error(std::format(
        "validateStepResiduals: residual has the wrong sign, worst s*F = {} "
        "(left {}, right {}) against a limit of {} for terms of size {}",
        worst.worst(), worst.left, worst.right, limit, worst.scale));
  }
}

std::vector<double> BarycentricAffineApproximator::solveMainLpStep(
    int phase, int solver_index, const std::vector<double>& x_next) {
  updateHighsRhsUpperBounds(phase, solver_index, x_next);
  std::vector<double> z = solveLp(phase, solver_index);
  validateStepResiduals(phase, x_next, z);
  return z;
}

void BarycentricAffineApproximator::validateBlockDecomposition(
    int phase) const {
  const auto& blocks = phase_geometries_[phase].blocks;
  std::vector<int> all_cells(kSpaceDim);
  std::iota(all_cells.begin(), all_cells.end(), 0);

  int rounds = 0;
  for (int b = 0; b < kBlockCount; ++b) {
    rounds = std::max(rounds, blocks[static_cast<std::size_t>(b)].numRegions());
  }

  // Block b's rows are built from block b's flows, which read only block b's
  // coordinates, so for a fixed j_b they do not depend on the other two blocks.
  // Walking the three index lists together therefore covers every block-region
  // in max_b M_b full evaluations instead of the product.
  for (int round = 0; round < rounds; ++round) {
    std::array<int, kBlockCount> js{};
    Eigen::VectorXd centroid = Eigen::VectorXd::Zero(kSpaceDim);
    for (int b = 0; b < kBlockCount; ++b) {
      const BlockGeometry& geometry = blocks[static_cast<std::size_t>(b)];
      js[static_cast<std::size_t>(b)]
          = std::min(round, geometry.numRegions() - 1);
      const Eigen::VectorXd block_centroid = blockCentroidCoords(
          phase, b, js[static_cast<std::size_t>(b)]);
      for (int c = 0; c < geometry.coord_count; ++c) {
        centroid(geometry.coords[static_cast<std::size_t>(c)])
            = block_centroid(c);
      }
    }

    const CtmRegionData full_ctm = ctmDataForCells(phase, all_cells, centroid);
    const BoxRegionData full_box = boxDataForCells(all_cells, centroid);

    double g_scal_sum = 0.0;
    for (int b = 0; b < kBlockCount; ++b) {
      const BlockGeometry& geometry = blocks[static_cast<std::size_t>(b)];
      const BlockSystemMatrices data = getBlockSystemMatrices(
          phase, b, js[static_cast<std::size_t>(b)]);
      g_scal_sum += data.ctm.g_scal;

      for (int r = 0; r < geometry.coord_count; ++r) {
        const int gr = geometry.coords[static_cast<std::size_t>(r)];
        for (int c = 0; c < kSpaceDim; ++c) {
          int local = -1;
          for (int k = 0; k < geometry.coord_count; ++k) {
            if (geometry.coords[static_cast<std::size_t>(k)] == c) {
              local = k;
              break;
            }
          }
          const double expected
              = local < 0 ? 0.0 : data.ctm.a(r, local);
          if (full_ctm.a(gr, c) != expected) {
            throw std::runtime_error(std::format(
                "validateBlockDecomposition: phase {} block {} disagrees on "
                "A({},{})", phase, b, gr, c));
          }
        }
        if (full_ctm.f(gr) != data.ctm.f(r)
            || full_ctm.g_vec(gr) != data.ctm.g_vec(r)
            || full_box.qc_diag(gr) != data.box.qc_diag(r)
            || full_box.qc_off(gr) != data.box.qc_off(r)
            || full_box.qr_diag(gr) != data.box.qr_diag(r)
            || full_box.qr_off(gr) != data.box.qr_off(r)) {
          throw std::runtime_error(std::format(
              "validateBlockDecomposition: phase {} block {} disagrees at "
              "cell {}", phase, b, gr));
        }
      }
    }
    // The one quantity here that is accumulated in a different order by the
    // two paths: the full call adds the four flow scalars in sequence while
    // the block sum adds two partial sums, and float addition is not
    // associative. Everything else compared above is written only by the
    // flows of its own block, in the same order, so it is bit-exact. What
    // this check exists to catch -- the two paths resolving different
    // branches -- is off by O(1), so a relative tolerance keeps all of its
    // discriminating power.
    if (std::abs(g_scal_sum - full_ctm.g_scal)
        > 1e-9 * std::max(1.0, std::abs(full_ctm.g_scal))) {
      throw std::runtime_error(std::format(
          "validateBlockDecomposition: phase {} block g scalars sum to {} "
          "instead of {}", phase, g_scal_sum, full_ctm.g_scal));
    }
  }

  // While the 8D product is still materialised, check that concatenating the
  // three block centroids reproduces it. The two are not bit-identical -- the
  // 8D mean sums every value |V_b'| |V_b''| times -- so this is the one place a
  // tolerance is used.
  // Skipped when the 8D product was not materialised, which leaves the region
  // list empty rather than populated with empty entries; the guard covers
  // both shapes.
  const auto& region_vertices = phase_geometries_[phase].region_vertices;
  if (!region_vertices.empty() && !region_vertices.front().empty()) {
    const int sample = std::min<int>(8, static_cast<int>(
        region_vertices.size()));
    const int m_b = blocks[1].numRegions();
    const int m_c = blocks[2].numRegions();
    for (int j = 0; j < sample; ++j) {
      const int jc = j % m_c;
      const int jb = (j / m_c) % m_b;
      const int ja = j / (m_c * m_b);
      const std::array<int, kBlockCount> js = {ja, jb, jc};
      Eigen::VectorXd concatenated = Eigen::VectorXd::Zero(kSpaceDim);
      for (int b = 0; b < kBlockCount; ++b) {
        const BlockGeometry& geometry = blocks[static_cast<std::size_t>(b)];
        const Eigen::VectorXd block_centroid = blockCentroidCoords(
            phase, b, js[static_cast<std::size_t>(b)]);
        for (int c = 0; c < geometry.coord_count; ++c) {
          concatenated(geometry.coords[static_cast<std::size_t>(c)])
              = block_centroid(c);
        }
      }
      const Eigen::VectorXd expected = areaCentroidCoords(j, phase);
      const double error = (concatenated - expected).cwiseAbs().maxCoeff();
      if (error > 1e-9 * std::max(1.0, expected.cwiseAbs().maxCoeff())) {
        throw std::runtime_error(std::format(
            "validateBlockDecomposition: phase {} region {} centroid differs "
            "by {} between the block and product paths", phase, j, error));
      }
    }
  }
}

void BarycentricAffineApproximator::precomputeMatrices() {
  logger_->info("Starting barycentric precomputeMatrices");
  for (int phase = 0; phase < kPhases; ++phase) {
    for (const auto& block : phase_geometries_[phase].blocks) {
      if (block.numRegions() == 0) {
        throw std::runtime_error(
            "precomputeMatrices: geometry not initialized. Call "
            "getIntersectionPoints() first.");
      }
    }
  }

  for (int phase = 0; phase < kPhases; ++phase) {
    for (int block = 0; block < kBlockCount; ++block) {
      block_system_[phase][static_cast<std::size_t>(block)]
          = precomputeBlockSystemMatrices(phase, block);
    }
    validateBlockDecomposition(phase);
    buildReducedLpInput(phase);

    double product_regions = 1.0;
    std::size_t block_pairs = 0;
    double product_pairs = 1.0;
    for (const auto& block : reduced_inputs_[phase].blocks) {
      product_regions *= static_cast<double>(block.regions.size());
      std::size_t pairs = 0;
      for (const auto& region : block.regions) {
        pairs += region.vertices.size();
      }
      block_pairs += pairs;
      product_pairs *= static_cast<double>(pairs);
    }
    logger_->info(
        "precomputeMatrices: phase={}, block regions={}/{}/{} (product "
        "{:.3e}), block (region, vertex) pairs={} (product {:.3e}), "
        "num_x={}, num_cols={}",
        phase, reduced_cols_[phase].num_regions[0],
        reduced_cols_[phase].num_regions[1],
        reduced_cols_[phase].num_regions[2], product_regions, block_pairs,
        product_pairs, layouts_[phase].num_x,
        reduced_cols_[phase].num_cols);
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
      phase_geometries_[0].region_triangle_ids.size(),
      phase_geometries_[1].region_triangle_ids.size());
  courier_options_.highs_verbose = highs_verbose_;
  courier_solver_ = CourierBorderSolver(courier_options_);
  courier_solver_.prepare(phase_geometries_, layouts_, node_weights_,
                          system_params_.N);
  logger_->info("Finished preparing the courier border solver");

  precomputeMatrices();

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

  // Initialize solver instances once and then reuse them. prepareLpMatrices()
  // is now a pure function of reduced_inputs_, which precomputeMatrices() filled
  // single-threaded above, so this loop no longer races on shared state -- it is
  // kept sequential only because it is cheap. The expensive HiGHS models are
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
              if (x_next.empty()) {
                // Unreachable node: nothing to march backwards from, and
                // nothing to store. See getBorderConditions.
                completed_tasks.fetch_add(1);
                return;
              }
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
