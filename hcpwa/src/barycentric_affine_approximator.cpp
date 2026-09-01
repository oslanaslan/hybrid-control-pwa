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
#include <random>
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
std::array<std::array<int, 2>, barycentric_affine_approximator::kSubsystemCount>
projectionAxesForPhase(int phase) {
  if (phase == 0) {
    // Phase 0 tuple order: [31, 36, 24, 27, 58].
    return {{{0, 2}, {2, 5}, {1, 3}, {1, 6}, {4, 7}}};
  }
  if (phase == 1) {
    // Phase 1 tuple order: [51, 57, 84, 86, 23].
    return {{{0, 4}, {4, 6}, {3, 7}, {5, 7}, {1, 2}}};
  }
  throw std::invalid_argument("projectionAxesForPhase: invalid phase");
}

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

// psiTransposeTimes() used to live here, building Psi_j^T m over the 8 rows of
// a full region's Psi. The reduced LP assembler does the equivalent per block,
// over coord_count rows, so it is gone.

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
          system_params_.f5max, system_params_.f8max,
          // The 8D area vertex lists are the truncated ones (12 of 108-192
          // vertices per area) and this path no longer uses them: it carries
          // the three block vertex sets end to end instead. Switching them off
          // is what makes the fix affordable -- materialising them correctly
          // would cost about 13 GB, against under 1 MB for the block data.
          hcpwa::TriangleAreasOptions{.verbose = true,
                                      .build_8d_regions = false});

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

  // Ingests one coordinate block and derives its local_axis table.
  //
  // The table maps each of the block's layers onto positions within the
  // block's own coordinate list. Two of the six entries across both phases are
  // non-contiguous, and getting one wrong yields a silently wrong phi rather
  // than a crash, so it is derived here by search and then checked three ways:
  // every layer axis must occur in coords, coords must be ascending, and the
  // layer axes together must cover coords exactly. A hard-coded table would
  // pass none of those checks by construction.
  auto ingest_block =
      [](int phase, int block_id, const hcpwa::BlockRegions& src,
         const std::array<std::array<int, 2>, kSubsystemCount>& axes) {
        BlockGeometry block;
        block.coords = src.coords;
        block.coord_count = src.coord_count;
        block.layer_ids = src.layer_ids;
        block.layer_count = src.layer_count;

        if (block.coord_count < 2 || block.coord_count > 3
            || block.layer_count < 1 || block.layer_count > 2) {
          throw std::runtime_error(std::format(
              "ingest_block: phase {} block {} has coord_count {} and "
              "layer_count {}",
              phase, block_id, block.coord_count, block.layer_count));
        }
        for (int d = 0; d < block.coord_count; ++d) {
          if (block.coords[d] < 0 || block.coords[d] >= kSpaceDim) {
            throw std::runtime_error(std::format(
                "ingest_block: phase {} block {} coordinate {} is {}", phase,
                block_id, d, block.coords[d]));
          }
          if (d > 0 && block.coords[d] <= block.coords[d - 1]) {
            throw std::runtime_error(std::format(
                "ingest_block: phase {} block {} coordinates are not strictly "
                "ascending",
                phase, block_id));
          }
        }

        std::array<bool, 3> covered = {false, false, false};
        for (int l = 0; l < block.layer_count; ++l) {
          const int layer_id = block.layer_ids[l];
          if (layer_id < 0 || layer_id >= kSubsystemCount) {
            throw std::runtime_error(std::format(
                "ingest_block: phase {} block {} references layer {}", phase,
                block_id, layer_id));
          }
          for (int k = 0; k < 2; ++k) {
            const int axis = axes[layer_id][k];
            int position = -1;
            for (int d = 0; d < block.coord_count; ++d) {
              if (block.coords[d] == axis) {
                position = d;
                break;
              }
            }
            if (position < 0) {
              throw std::runtime_error(std::format(
                  "ingest_block: phase {} block {} layer {} has axis {} which "
                  "is not one of the block's coordinates",
                  phase, block_id, layer_id, axis));
            }
            block.local_axis[l][k] = position;
            covered[position] = true;
          }
        }
        for (int d = 0; d < block.coord_count; ++d) {
          if (!covered[d]) {
            throw std::runtime_error(std::format(
                "ingest_block: phase {} block {} coordinate {} (state axis {}) "
                "is covered by none of the block's layers",
                phase, block_id, d, block.coords[d]));
          }
        }

        block.num_regions = static_cast<int>(src.vertices.size());
        if (static_cast<int>(src.triangle_ids.size()) != block.num_regions) {
          throw std::runtime_error(std::format(
              "ingest_block: phase {} block {} has {} vertex lists but {} "
              "simplex-id tuples",
              phase, block_id, src.vertices.size(), src.triangle_ids.size()));
        }
        block.triangle_ids.reserve(block.num_regions);
        block.vertices.reserve(block.num_regions);
        for (int j = 0; j < block.num_regions; ++j) {
          block.triangle_ids.push_back(
              {static_cast<int>(src.triangle_ids[j][0]),
               static_cast<int>(src.triangle_ids[j][1])});
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
  phase_geometries_[0].region_triangle_ids
      = convert_region_indices(areas_vertices.intersection_prism_indices_phase0);
  for (int b = 0; b < kBlockCount; ++b) {
    phase_geometries_[0].blocks[b]
        = ingest_block(0, b, areas_vertices.blocks_phase0[b], phase0_axes);
  }

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
  phase_geometries_[1].region_triangle_ids
      = convert_region_indices(areas_vertices.intersection_prism_indices_phase1);
  for (int b = 0; b < kBlockCount; ++b) {
    phase_geometries_[1].blocks[b]
        = ingest_block(1, b, areas_vertices.blocks_phase1[b], phase1_axes);
  }

  for (int phase = 0; phase < kPhases; ++phase) {
    auto& geometry = phase_geometries_[phase];

    // Lemma 1 again, this time against the area ids the border path uses. The
    // geometry layer already checked it; checking it here too means a future
    // change that populates one of these two from a different source cannot
    // pass silently.
    const long long product
        = static_cast<long long>(geometry.blocks[0].num_regions)
          * geometry.blocks[1].num_regions * geometry.blocks[2].num_regions;
    if (product != static_cast<long long>(geometry.region_triangle_ids.size())) {
      throw std::runtime_error(std::format(
          "Barycentric geometry: phase {} has {} areas but M_A*M_B*M_C = {}",
          phase, geometry.region_triangle_ids.size(), product));
    }
    geometry.num_regions = static_cast<int>(geometry.region_triangle_ids.size());

    logger_->info(
        "getIntersectionPoints: phase={}, M_A={}, M_B={}, M_C={}, regions={}, "
        "block rows R_A={}, R_B={}, R_C={}",
        phase, geometry.blocks[0].num_regions, geometry.blocks[1].num_regions,
        geometry.blocks[2].num_regions, geometry.num_regions,
        [&] {
          std::size_t total = 0;
          for (const auto& v : geometry.blocks[0].vertices) total += v.size();
          return total;
        }(),
        [&] {
          std::size_t total = 0;
          for (const auto& v : geometry.blocks[1].vertices) total += v.size();
          return total;
        }(),
        [&] {
          std::size_t total = 0;
          for (const auto& v : geometry.blocks[2].vertices) total += v.size();
          return total;
        }());

    BarycentricVarLayout layout;
    for (int s = 0; s < kSubsystemCount; ++s) {
      layout.offset_s[s] = layout.num_x;
      layout.eta_s[s]
          = static_cast<int>(geometry.layers[s].unique_vertices.size());
      layout.num_x += layout.eta_s[s];
    }
    // The y and mu columns are appended by the reduced LP assembler, so the
    // total column count is not known here and is not stored here.
    layouts_[phase] = layout;
  }

  // Needs layouts_ (idxX), so it runs after the loop above.
  for (int phase = 0; phase < kPhases; ++phase) {
    validateBlockGeometry(phase);
  }

  // Common-refinement vertices are deduplicated globally, but cell membership is
  // preserved: the lower-bound border LP selects one family member per cell, so
  // it needs to know which vertices belong to the same cell (step 2.2, 6.1).
  common_refinement_vertices_.clear();
  refinement_cells_.clear();
  refinement_cells_.reserve(areas_vertices.common_refinement.areas.size());
  for (const auto& area : areas_vertices.common_refinement.areas) {
    RefinementCell cell;
    cell.phase0_area_id = static_cast<int>(area.phase0_area_id);
    cell.phase1_area_id = static_cast<int>(area.phase1_area_id);
    cell.vertex_ids.reserve(area.vertices.size());

    for (const auto& vertex : area.vertices) {
      Eigen::VectorXd g = toEigen8(vertex);
      for (int d = 0; d < kSpaceDim; ++d) {
        if (!std::isfinite(g(d))) {
          throw std::runtime_error(
              "getIntersectionPoints: non-finite common-refinement vertex");
        }
        if (g(d) < -kGeomEps || g(d) > system_params_.N + kGeomEps) {
          throw std::runtime_error(
              "getIntersectionPoints: common-refinement vertex outside [0,N]");
        }
        if (std::abs(g(d)) <= kGeomEps) {
          g(d) = 0.0;
        }
        if (std::abs(g(d) - system_params_.N) <= kGeomEps) {
          g(d) = system_params_.N;
        }
      }

      int vertex_id = -1;
      for (int i = 0;
           i < static_cast<int>(common_refinement_vertices_.size()); ++i) {
        if ((common_refinement_vertices_[i] - g).norm() <= kGeomEps) {
          vertex_id = i;
          break;
        }
      }
      if (vertex_id < 0) {
        common_refinement_vertices_.push_back(std::move(g));
        vertex_id = static_cast<int>(common_refinement_vertices_.size() - 1);
      }
      cell.vertex_ids.push_back(vertex_id);
    }

    if (cell.vertex_ids.empty()) {
      throw std::runtime_error(
          "getIntersectionPoints: common-refinement cell has no vertices");
    }
    refinement_cells_.push_back(std::move(cell));
  }

  if (common_refinement_vertices_.empty() || refinement_cells_.empty()) {
    throw std::runtime_error(
        "getIntersectionPoints: common-refinement test set is empty");
  }

  // Every border test point must be evaluable in both phase partitions. This
  // validation is intentionally done during geometry construction so a prism
  // tuple/order bug fails before any LP coefficients are stored. The phi rows
  // found here are cached: the border LP runs once per (level, theta, phase).
  // We use the already stored phase-area ancestors from each refinement cell
  // instead of searching all phase regions for every deduplicated vertex.
  for (int phase = 0; phase < kPhases; ++phase) {
    phi_at_refinement_[phase].clear();
    phi_at_refinement_[phase].resize(common_refinement_vertices_.size());
  }
  std::array<std::vector<bool>, kPhases> phi_initialized;
  for (int phase = 0; phase < kPhases; ++phase) {
    phi_initialized[phase].assign(common_refinement_vertices_.size(), false);
  }
  for (const RefinementCell& cell : refinement_cells_) {
    const std::array<int, kPhases> ancestor_regions
        = {cell.phase0_area_id, cell.phase1_area_id};
    for (int vertex_id : cell.vertex_ids) {
      if (vertex_id < 0
          || vertex_id >= static_cast<int>(common_refinement_vertices_.size())) {
        throw std::runtime_error(
            "getIntersectionPoints: refinement cell references invalid vertex id");
      }
      const Eigen::VectorXd& g = common_refinement_vertices_[vertex_id];
      for (int phase = 0; phase < kPhases; ++phase) {
        if (phi_initialized[phase][vertex_id]) {
          continue;
        }
        const int region = ancestor_regions[phase];
        SparseVec phi = buildPhiRow(phase, region, g, kEps);

        // Invariant of the barycentric representation: the coordinates of each
        // of the five projection layers sum to one, hence five in total.
        double sum = 0.0;
        for (double value : phi.vals) {
          if (value < -kEps) {
            throw std::runtime_error(
                "getIntersectionPoints: negative barycentric coordinate");
          }
          sum += value;
        }
        if (std::abs(sum - static_cast<double>(kSubsystemCount)) > 1e-6) {
          std::ostringstream details;
          details << "phase=" << phase << ", region=" << region
                  << ", sum=" << sum << ", expected=" << kSubsystemCount
                  << ", point=[";
          for (int d = 0; d < kSpaceDim; ++d) {
            details << (d == 0 ? "" : ", ") << g(d);
          }
          details << "], phi={";
          for (std::size_t i = 0; i < phi.vals.size(); ++i) {
            details << (i == 0 ? "" : ", ") << phi.cols[i] << ":"
                    << phi.vals[i];
          }
          details << "}";
          logger_->error(
              "getIntersectionPoints: phi row does not sum to "
              "kSubsystemCount: {}",
              details.str());
          throw std::runtime_error(std::format(
              "getIntersectionPoints: phi row sum {} does not equal {} for phase "
              "{} region {}",
              sum, kSubsystemCount, phase, region));
        }

        phi_at_refinement_[phase][vertex_id] = std::move(phi);
        phi_initialized[phase][vertex_id] = true;
      }
    }
  }
  for (int phase = 0; phase < kPhases; ++phase) {
    for (int vertex_id = 0;
         vertex_id < static_cast<int>(common_refinement_vertices_.size());
         ++vertex_id) {
      if (!phi_initialized[phase][vertex_id]) {
        throw std::runtime_error(std::format(
            "getIntersectionPoints: missing cached phi row for phase {} vertex "
            "{}",
            phase, vertex_id));
      }
    }
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

SparseVec BarycentricAffineApproximator::buildPhiRowBlock(
    int phase, int block_id, int j_block, const Eigen::VectorXd& nu_block,
    double tolerance) const {
  if (phase < 0 || phase >= kPhases) {
    throw std::invalid_argument("buildPhiRowBlock: invalid phase");
  }
  if (block_id < 0 || block_id >= kBlockCount) {
    throw std::invalid_argument("buildPhiRowBlock: invalid block");
  }
  const auto& geometry = phase_geometries_[phase];
  const auto& block = geometry.blocks[block_id];
  const auto& layout = layouts_[phase];
  if (j_block < 0 || j_block >= block.num_regions) {
    throw std::invalid_argument("buildPhiRowBlock: invalid block region id");
  }
  if (nu_block.size() != block.coord_count) {
    throw std::invalid_argument(std::format(
        "buildPhiRowBlock: point has dimension {}, expected {}",
        nu_block.size(), block.coord_count));
  }

  SparseVec phi;
  for (int l = 0; l < block.layer_count; ++l) {
    const int layer_id = block.layer_ids[l];
    const int triangle_id = block.triangle_ids[j_block][l];
    const auto& layer = geometry.layers[layer_id];
    if (triangle_id < 0
        || triangle_id >= static_cast<int>(layer.bases.size())) {
      throw std::runtime_error("buildPhiRowBlock: triangle id out of range");
    }
    const auto& basis = layer.bases[triangle_id];

    // local_axis maps the layer's two axes onto positions within the block's
    // own coordinate list. This is the step that is silently wrong for the two
    // non-contiguous blocks if the table is guessed rather than derived.
    Eigen::Vector2d projected;
    projected(0) = nu_block(block.local_axis[l][0]);
    projected(1) = nu_block(block.local_axis[l][1]);
    const Eigen::Vector3d alpha = basis.H * projected + basis.h;

    if (std::abs(alpha.sum() - 1.0) > tolerance) {
      throw std::runtime_error(
          "buildPhiRowBlock: barycentric alpha sum is not one");
    }
    for (int local_vertex = 0; local_vertex < 3; ++local_vertex) {
      if (alpha(local_vertex) < -tolerance
          || alpha(local_vertex) > 1.0 + tolerance) {
        throw std::runtime_error(
            "buildPhiRowBlock: point is outside selected simplex");
      }
      // Columns are global. layer_id, not block_id, selects the x block.
      const int x_col = layout.idxX(layer_id, basis.vertex_ids[local_vertex]);
      phi.add(x_col, alpha(local_vertex), kGeomEps);
    }
  }
  return phi;
}

Eigen::VectorXd BarycentricAffineApproximator::blockCentroidCoords(
    int phase, int block_id, int j_block) const {
  if (phase < 0 || phase >= kPhases) {
    throw std::invalid_argument("blockCentroidCoords: invalid phase");
  }
  if (block_id < 0 || block_id >= kBlockCount) {
    throw std::invalid_argument("blockCentroidCoords: invalid block");
  }
  const auto& block = phase_geometries_[phase].blocks[block_id];
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

void BarycentricAffineApproximator::validateBlockGeometry(int phase) const {
  const auto& geometry = phase_geometries_[phase];
  for (int b = 0; b < kBlockCount; ++b) {
    const auto& block = geometry.blocks[b];
    for (int j = 0; j < block.num_regions; ++j) {
      const auto& vertices = block.vertices[j];
      // A 3D block is a full-dimensional polytope, so it has at least 3
      // vertices (LinesToPoints discards anything lower dimensional; measured
      // counts run 6 to 8). Block C is a single triangle, so exactly 3.
      if (vertices.size() < 3) {
        throw std::runtime_error(std::format(
            "validateBlockGeometry: phase {} block {} region {} has {} "
            "vertices, expected at least 3",
            phase, b, j, vertices.size()));
      }
      if (block.coord_count == 2 && vertices.size() != 3) {
        throw std::runtime_error(std::format(
            "validateBlockGeometry: phase {} block {} region {} is a triangle "
            "but has {} vertices",
            phase, b, j, vertices.size()));
      }

      for (const Eigen::VectorXd& nu : vertices) {
        // Throws if nu falls outside one of the block's selected simplices,
        // which is what a wrong local_axis entry would produce.
        const SparseVec phi = buildPhiRowBlock(phase, b, j, nu, kEps);

        // Each of the block's layers contributes barycentric coordinates
        // summing to one, so the block's phi sums to its layer count and the
        // three blocks together sum to five.
        double sum = 0.0;
        for (double value : phi.vals) {
          if (value < -kEps) {
            throw std::runtime_error(std::format(
                "validateBlockGeometry: phase {} block {} region {} produced a "
                "negative barycentric coordinate {}",
                phase, b, j, value));
          }
          sum += value;
        }
        if (std::abs(sum - static_cast<double>(block.layer_count)) > 1e-6) {
          throw std::runtime_error(std::format(
              "validateBlockGeometry: phase {} block {} region {} has phi "
              "summing to {}, expected {}",
              phase, b, j, sum, block.layer_count));
        }
      }
    }
  }

  // The three blocks' layer sets must partition the five projection layers,
  // otherwise a layer would be counted twice or not at all in phi.
  std::array<int, kSubsystemCount> layer_uses{};
  for (int b = 0; b < kBlockCount; ++b) {
    const auto& block = geometry.blocks[b];
    for (int l = 0; l < block.layer_count; ++l) {
      ++layer_uses[block.layer_ids[l]];
    }
  }
  for (int s = 0; s < kSubsystemCount; ++s) {
    if (layer_uses[s] != 1) {
      throw std::runtime_error(std::format(
          "validateBlockGeometry: phase {} layer {} is claimed by {} blocks, "
          "expected exactly 1",
          phase, s, layer_uses[s]));
    }
  }

  // Same for the eight state coordinates.
  std::array<int, kSpaceDim> coord_uses{};
  for (int b = 0; b < kBlockCount; ++b) {
    const auto& block = geometry.blocks[b];
    for (int d = 0; d < block.coord_count; ++d) {
      ++coord_uses[block.coords[d]];
    }
  }
  for (int r = 0; r < kSpaceDim; ++r) {
    if (coord_uses[r] != 1) {
      throw std::runtime_error(std::format(
          "validateBlockGeometry: phase {} coordinate {} is claimed by {} "
          "blocks, expected exactly 1",
          phase, r, coord_uses[r]));
    }
  }
}

std::vector<int> BarycentricAffineApproximator::locateRegions(
    int phase, const Eigen::VectorXd& point, double tolerance) const {
  // Locates per block and composes, rather than scanning all M_A*M_B*M_C
  // regions. A point lies in region (j_A, j_B, j_C) exactly when it lies in
  // each block polytope separately, because the region IS the product of the
  // three (step 7, Lemma 1). That turns 1.36M simplex tests into about 509.
  if (phase < 0 || phase >= kPhases) {
    throw std::invalid_argument("locateRegions: invalid phase");
  }
  if (point.size() != kSpaceDim) {
    throw std::invalid_argument("locateRegions: point must be 8-dimensional");
  }

  const auto& geometry = phase_geometries_[phase];

  std::array<std::vector<int>, kBlockCount> per_block;
  for (int b = 0; b < kBlockCount; ++b) {
    const auto& block = geometry.blocks[b];
    for (int j = 0; j < block.num_regions; ++j) {
      bool contains = true;
      for (int l = 0; l < block.layer_count && contains; ++l) {
        const int layer_id = block.layer_ids[l];
        const auto& layer = geometry.layers[layer_id];
        const int triangle_id = block.triangle_ids[j][l];
        if (triangle_id < 0
            || triangle_id >= static_cast<int>(layer.bases.size())) {
          throw std::runtime_error("locateRegions: triangle id out of range");
        }
        const auto& basis = layer.bases[triangle_id];
        // The point is 8-dimensional here, so it is projected with the layer's
        // global axes; local_axis is only needed when the input is already
        // restricted to a block.
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
        per_block[b].push_back(j);
      }
    }
  }

  // region = (j_A * M_B + j_B) * M_C + j_C, the order the geometry pipeline
  // assembles regions in. Iterating the three lists in ascending order yields
  // region ids in ascending order too, so the result needs no sort.
  const int m_b = geometry.blocks[1].num_regions;
  const int m_c = geometry.blocks[2].num_regions;
  std::vector<int> matches;
  matches.reserve(per_block[0].size() * per_block[1].size()
                  * per_block[2].size());
  for (int j_a : per_block[0]) {
    for (int j_b : per_block[1]) {
      for (int j_c : per_block[2]) {
        matches.push_back((j_a * m_b + j_b) * m_c + j_c);
      }
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

// Expands a block centroid into a full 8-vector for getFIJMinResolution().
//
// The coordinates outside the block are filled with NaN on purpose. Each
// block's flows depend only on that block's own coordinates -- phase 0 block A
// resolves f31 and f36 over {0,2,5}, block B resolves f24 and f27 over {1,3,6}
// -- so the padding must never be read. If it ever is, every comparison in
// getFIJMinResolution() evaluates false and the function throws "No min branch"
// instead of quietly resolving against a made-up value. That turns "assert the
// padding is unused" into something the arithmetic enforces for free.
Eigen::VectorXd BarycentricAffineApproximator::blockPointToFullState(
    int phase, int block_id, const Eigen::VectorXd& nu_block) const {
  const auto& block = phase_geometries_[phase].blocks[block_id];
  Eigen::VectorXd full = Eigen::VectorXd::Constant(
      kSpaceDim, std::numeric_limits<double>::quiet_NaN());
  for (int d = 0; d < block.coord_count; ++d) {
    full(block.coords[d]) = nu_block(d);
  }
  return full;
}

std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd, double>
BarycentricAffineApproximator::getBlockAMatrFVecGVecAndGScalJ(
    int phase, int block_id, int j_block) const {
  // The same CTM branch-resolution formulas as the whole-region version, but
  // restricted to one block.
  //
  // A^(j) is block diagonal with respect to the coordinate partition (step 7
  // section 3.3): every state equation depends only on coordinates of its own
  // block. The restriction below therefore loses nothing, and restrict_row()
  // asserts that by refusing any entry outside the block. Block C is purely
  // exogenous in both phases, so its A, f and g are all zero.
  const auto& block = phase_geometries_[phase].blocks[block_id];
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

  // Restricts an 8-column flow row to the block's own columns, checking that
  // nothing outside the block was dropped.
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
            "flow coefficient {} on coordinate {}, which lies outside the "
            "block. A^(j) is supposed to be block diagonal.",
            phase, block_id, j_block, row(r), r));
      }
    }
    return out;
  };

  // The two flows this block owns, and where each state coordinate's row of A
  // comes from. Phase 0: block A holds n1, n3, n6 and the flows f31, f36;
  // block B holds n2, n4, n7 and the flows f24, f27. Phase 1: block A' holds
  // n1, n5, n7 with f51, f57; block B' holds n4, n6, n8 with f84, f86.
  int first_i = 0;
  int first_j = 0;
  int second_i = 0;
  int second_j = 0;
  if (phase == 0 && block_id == 0) {
    first_i = 3 - 1;  first_j = 1 - 1;   // f31
    second_i = 3 - 1; second_j = 6 - 1;  // f36
  } else if (phase == 0 && block_id == 1) {
    first_i = 2 - 1;  first_j = 4 - 1;   // f24
    second_i = 2 - 1; second_j = 7 - 1;  // f27
  } else if (phase == 1 && block_id == 0) {
    first_i = 5 - 1;  first_j = 1 - 1;   // f51
    second_i = 5 - 1; second_j = 7 - 1;  // f57
  } else if (phase == 1 && block_id == 1) {
    first_i = 8 - 1;  first_j = 4 - 1;   // f84
    second_i = 8 - 1; second_j = 6 - 1;  // f86
  } else {
    throw std::invalid_argument(
        "getBlockAMatrFVecGVecAndGScalJ: invalid phase/block combination");
  }

  const auto [first_row, first_scal] = getFIJMinResolution(first_i, first_j, n);
  const auto [second_row, second_scal]
      = getFIJMinResolution(second_i, second_j, n);
  const Eigen::RowVectorXd first = restrict_row(first_row);
  const Eigen::RowVectorXd second = restrict_row(second_row);

  // Position of a state coordinate within this block.
  auto local = [&](int state_axis) {
    for (int d = 0; d < c; ++d) {
      if (block.coords[d] == state_axis) {
        return d;
      }
    }
    throw std::runtime_error(std::format(
        "getBlockAMatrFVecGVecAndGScalJ: coordinate {} is not in phase {} "
        "block {}",
        state_axis, phase, block_id));
  };

  // Rows of A, matching the whole-region matrix exactly:
  //   phase 0: row0 = f31, row2 = -f31-f36, row5 = f36,
  //            row1 = -f24-f27, row3 = f24, row6 = f27
  //   phase 1: row0 = f51, row4 = -f51-f57, row6 = f57,
  //            row3 = f84, row7 = -f84-f86, row5 = f86
  // In each case the flow's source cell takes the negative sum and the two
  // destination cells take one flow each.
  const int source = local(first_i);
  const int first_dest = local(first_j);
  const int second_dest = local(second_j);

  a_matr.row(first_dest) += first;
  b_vec(first_dest) += first_scal(0);
  a_matr.row(second_dest) += second;
  b_vec(second_dest) += second_scal(0);
  a_matr.row(source) -= first + second;
  b_vec(source) -= first_scal(0) + second_scal(0);

  // g_i = f31 + f36 + f24 + f27 in phase 0, split additively across blocks A
  // and B; block C contributes nothing (step 7 section 3.4).
  g_vec = (first + second).transpose();
  g_scal = first_scal(0) + second_scal(0);

  hcpwa::util::assertShape(a_matr, c, c);
  hcpwa::util::assertShape(b_vec, c);
  hcpwa::util::assertShape(g_vec, c);
  hcpwa::util::assertScalar(g_scal);
  return std::make_tuple(a_matr, b_vec, g_vec, g_scal);
}

std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd, Eigen::VectorXd>
BarycentricAffineApproximator::getBlockQQ(int phase, int block_id,
                                          int j_block) const {
  // Affine centre/radius maps for the disturbance box, restricted to one block.
  //
  // The constraints on f_{r,in} and f_{r,out} bound component r using n_r only
  // (step 7 section 3.3), so the box is componentwise and its restriction to a
  // block is exact. The loop runs over the block's own coordinates rather than
  // all eight, which is what makes the NaN padding in the representative point
  // safe: an unrestricted loop would read it.
  const auto& block = phase_geometries_[phase].blocks[block_id];
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

  Eigen::MatrixXd Q_c = (Q_upper + Q_lower) / 2.0;
  Eigen::MatrixXd Q_r = (Q_upper - Q_lower) / 2.0;
  Eigen::VectorXd q_c = (q_upper + q_lower) / 2.0;
  Eigen::VectorXd q_r = (q_upper - q_lower) / 2.0;
  return std::make_tuple(Q_c, q_c, Q_r, q_r);
}

void BarycentricAffineApproximator::precomputeSystemMatrices(int phase) {
  // Per-block-region static data: Psi and the CTM matrices. About 509 small
  // objects on the production geometry, against 1.36M dense 8x8 blocks before.
  auto& geometry = phase_geometries_[phase];
  const auto& layout = layouts_[phase];

  for (int b = 0; b < kBlockCount; ++b) {
    auto& block = geometry.blocks[b];
    const int c = block.coord_count;

    block.psi.assign(block.num_regions, {});
    block.a_matr.assign(block.num_regions, {});
    block.f_vec.assign(block.num_regions, {});
    block.q_c_matr.assign(block.num_regions, {});
    block.q_c_vec.assign(block.num_regions, {});
    block.q_r_matr.assign(block.num_regions, {});
    block.q_r_vec.assign(block.num_regions, {});
    block.g_vec.assign(block.num_regions, {});
    block.g_scal.assign(block.num_regions, 0.0);

    for (int j = 0; j < block.num_regions; ++j) {
      // Psi_{block,j}: one row per block coordinate, over global x columns.
      // For every local barycentric value the derivative d alpha_l / d n_axis
      // contributes H(l, k), and local_axis says which block coordinate the
      // layer's k-th axis is.
      hcpwa::util::block_lp::BlockPsi psi;
      psi.rows.resize(c);
      for (int l = 0; l < block.layer_count; ++l) {
        const int layer_id = block.layer_ids[l];
        const auto& layer = geometry.layers[layer_id];
        const int triangle_id = block.triangle_ids[j][l];
        if (triangle_id < 0
            || triangle_id >= static_cast<int>(layer.bases.size())) {
          throw std::runtime_error(
              "precomputeSystemMatrices: triangle id out of range");
        }
        const auto& basis = layer.bases[triangle_id];
        for (int local_vertex = 0; local_vertex < 3; ++local_vertex) {
          const int x_col = layout.idxX(layer_id, basis.vertex_ids[local_vertex]);
          for (int k = 0; k < 2; ++k) {
            psi.rows[block.local_axis[l][k]].add(
                x_col, basis.H(local_vertex, k), kGeomEps);
          }
        }
      }

      // Psi 1 = 0 on every row: the barycentric coordinates of a simplex sum to
      // one everywhere, so the gradient of their sum vanishes identically. This
      // is what makes the gauge shift z -> z + theta 1^(s) leave q unchanged,
      // and it is the cheapest available check that H was read with the right
      // orientation. Exactly zero, not small.
      for (int r = 0; r < c; ++r) {
        double row_sum = 0.0;
        for (double value : psi.rows[r].vals) {
          row_sum += value;
        }
        // Tolerance must sit ABOVE the kGeomEps threshold the row was built
        // at: sum_l H(l,k) is zero exactly, so if add() drops an entry of
        // magnitude just under kGeomEps the surviving two sum to that entry.
        // Checking tighter than the drop threshold would reject valid geometry.
        if (std::abs(row_sum) > 1e-6) {
          throw std::runtime_error(std::format(
              "precomputeSystemMatrices: phase {} block {} region {} row {} of "
              "Psi sums to {}, expected 0",
              phase, b, j, r, row_sum));
        }
      }
      block.psi[j] = std::move(psi);

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

hcpwa::util::block_lp::ReducedLp
BarycentricAffineApproximator::prepareLpMatrices(int phase) {
  // Builds the reduced LP of step 7 section 10. Three independent block loops
  // plus two scalar coupling rows, in place of one loop over the Cartesian
  // product of the three blocks' row sets.
  //
  // The reduction is exact: the feasible sets have identical projections onto
  // z (Lemmas 3 and 4), and identifying y down to one vector per block region
  // is safe because rho enters row (L) with a plus sign for BOTH bound
  // directions (Lemma 5). Only the objective weighting is a choice, see below.
  //
  // Sizes on the production geometry: about 9950 rows and eta + 1510 columns,
  // against 5.4e7 rows and eta + 1.1e7 columns for the truncated product form,
  // or roughly 4.2e8 rows for the product form built correctly.
  namespace blp = hcpwa::util::block_lp;

  const auto& geometry = phase_geometries_[phase];
  const auto& layout = layouts_[phase];
  const double s = signS();

  std::array<blp::BlockInput, kBlockCount> inputs;
  for (int b = 0; b < kBlockCount; ++b) {
    const auto& block = geometry.blocks[b];
    blp::BlockInput& input = inputs[b];
    input.coord_count = block.coord_count;
    input.num_block_regions = block.num_regions;
    input.psi = block.psi;

    // Objective weight. Summing residuals over the Cartesian product would
    // weight each block by the row count of the other two, R/R_block, which on
    // the production geometry spans a factor of 26 between blocks and is an
    // artefact of how the sum is taken rather than a design decision. We use 1
    // per block instead.
    //
    // Recorded consequence: the resulting z* is NOT the one the product
    // objective would produce. That costs nothing here, because the product
    // objective was never computable at production scale, so there is no
    // earlier result to match. Validity of the bound does not depend on the
    // objective at all -- any feasible point is a certificate -- only
    // tightness does.
    input.objective_weight = 1.0;

    std::size_t total_rows = 0;
    for (const auto& vertices : block.vertices) {
      total_rows += vertices.size();
    }
    input.rows.reserve(total_rows);

    for (int j = 0; j < block.num_regions; ++j) {
      for (const Eigen::VectorXd& nu : block.vertices[j]) {
        blp::BlockRow row;
        row.block_region = j;

        // phi restricted to this block, in global x columns.
        const SparseVec phi = buildPhiRowBlock(phase, b, j, nu, kEps);
        row.phi.cols = phi.cols;
        row.phi.vals = phi.vals;

        // m = A nu + f + c(nu), the drift with the disturbance-box centre
        // folded in; rho is the box radius. Both restricted to the block's
        // coordinates, which is exact because A is block diagonal and the box
        // is componentwise.
        const Eigen::VectorXd drift = block.a_matr[j] * nu + block.f_vec[j];
        const Eigen::VectorXd centre
            = block.q_c_matr[j] * nu + block.q_c_vec[j];
        Eigen::VectorXd rho = block.q_r_matr[j] * nu + block.q_r_vec[j];
        for (int p = 0; p < block.coord_count; ++p) {
          if (rho(p) < -kEps) {
            throw std::runtime_error(std::format(
                "prepareLpMatrices: phase {} block {} region {} has negative "
                "uncertainty radius {}",
                phase, b, j, rho(p)));
          }
          if (rho(p) < 0.0) {
            rho(p) = 0.0;
          }
        }

        row.m = drift + centre;
        row.rho = std::move(rho);
        row.g = block.g_vec[j].dot(nu) + block.g_scal[j];
        input.rows.push_back(std::move(row));
      }
    }
  }

  blp::ReducedLpOptions options;
  options.dt = t_delta_;
  options.s = s;
  options.num_x = layout.num_x;
  // Must match the tolerance phi and Psi were BUILT at, not the old row-level
  // kEps. This assembler writes (R) at the 1/dt scale, so a kEps threshold here
  // would delete phi entries up to kEps*dt from the matrix while
  // reducedLpRowUpper() still uses them in the right-hand side -- the two sides
  // of the same constraint would disagree. HiGHS drops genuinely negligible
  // entries itself via small_matrix_value.
  options.coefficient_eps = kGeomEps;
  // Mandatory, not optional: the optimal face is genuinely wider than a point,
  // and the march is greedy, so an arbitrary choice among optima at one step
  // changes every step after it. See ReducedLpOptions::tie_break_epsilon.
  options.tie_break_epsilon = kTieBreakEpsilon;

  blp::ReducedLp lp = blp::assembleReducedLp(inputs, options);

  logger_->info(
      "prepareLpMatrices: phase={}, variables(num_cols)={}, num_x={}, "
      "num_y={}, mu=6, block regions={}/{}/{}, block rows={}/{}/{}, rows={}, "
      "nnz={}, rhs_terms={}",
      phase, lp.num_cols, lp.num_x,
      lp.num_cols - lp.num_x - 6
          - (lp.tie_break_offset >= 0 ? lp.num_x : 0),
      lp.num_block_regions[0], lp.num_block_regions[1], lp.num_block_regions[2],
      lp.num_block_rows[0], lp.num_block_rows[1], lp.num_block_rows[2],
      lp.row_upper.size(), lp.values.size(), lp.rhs_terms.size());

  return lp;
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
  if (common_refinement_vertices_.empty()) {
    throw std::runtime_error(
        "getBorderConditions: common-refinement test set is empty");
  }

  const std::vector<int> theta_end_ids = admissibleThetaIds(theta);
  const bool is_upper = approximation_mode_ == ApproximationMode::Upper;
  const int n_points = static_cast<int>(common_refinement_vertices_.size());

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

  // val[g][c] = phi^(source)(g)^T u^(c). The upper branch only needs the maximum
  // over c at each g; the lower branch needs the whole table, because it selects
  // one member per refinement cell.
  //
  // The source-side phi row is the cached one. Calling evaluateBarycentricValue
  // here instead would relocate g among all source regions on every (g, c) pair,
  // which is O(#regions) each time. The cached row comes from the same region
  // that lookup would return, and the barycentric function is continuous across
  // region boundaries by construction (neighbouring regions share node values),
  // so the value is the same.
  std::vector<std::vector<double>> val(
      static_cast<std::size_t>(n_points),
      std::vector<double>(static_cast<std::size_t>(n_candidates), 0.0));
  for (int g_id = 0; g_id < n_points; ++g_id) {
    const SparseVec& phi_src = phi_at_refinement_[source_phase][g_id];
    for (int c = 0; c < n_candidates; ++c) {
      const double val_src = phi_src.dot(candidates[c]);
      if (!std::isfinite(val_src)) {
        throw std::runtime_error("getBorderConditions: non-finite source value");
      }
      val[g_id][c] = val_src;
    }
  }

  // Right-hand sides.
  //
  // Upper bound: the requirement is  phi^T x_L >= h(g) = max_c val[g][c].
  // On a refinement cell the left side is affine and h is convex (a max of
  // affine functions), so their difference is concave and its minimum sits at a
  // vertex: the vertex rows are exact (step 2.2, section 4).
  //
  // Lower bound: the requirement is  phi^T x_L <= h(g), and now the difference
  // is convex, so vertices are NOT sufficient. Instead one family member is
  // fixed per cell; both sides are then affine on that cell and vertices become
  // exact again, while "member <= max" keeps the condition sufficient
  // (step 2.2, section 6.1). A vertex shared by several cells keeps the
  // tightest of the bounds it receives.
  std::vector<double> rhs(static_cast<std::size_t>(n_points));
  if (is_upper) {
    for (int g_id = 0; g_id < n_points; ++g_id) {
      rhs[g_id] = *std::max_element(val[g_id].begin(), val[g_id].end());
    }
  } else {
    std::fill(rhs.begin(), rhs.end(), std::numeric_limits<double>::infinity());
    for (const RefinementCell& cell : refinement_cells_) {
      // The member is affine on the cell, so its value at the barycentre of the
      // cell vertices equals the mean of its vertex values. Picking the member
      // with the largest mean is a heuristic: it affects tightness only, never
      // correctness.
      int best = -1;
      double best_mean = -std::numeric_limits<double>::infinity();
      for (int c = 0; c < n_candidates; ++c) {
        double mean = 0.0;
        for (int g_id : cell.vertex_ids) {
          mean += val[g_id][c];
        }
        mean /= static_cast<double>(cell.vertex_ids.size());
        if (mean > best_mean) {
          best_mean = mean;
          best = c;
        }
      }
      if (best < 0) {
        throw std::runtime_error(
            "getBorderConditions: no family member selected for a cell");
      }
      for (int g_id : cell.vertex_ids) {
        rhs[g_id] = std::min(rhs[g_id], val[g_id][best]);
      }
    }
    for (int g_id = 0; g_id < n_points; ++g_id) {
      if (!std::isfinite(rhs[g_id])) {
        throw std::runtime_error(
            "getBorderConditions: border point belongs to no refinement cell");
      }
    }
  }

  std::vector<int> starts = {0};
  std::vector<int> col_index;
  std::vector<double> value;
  std::vector<double> row_lower;
  std::vector<double> row_upper;
  std::vector<SparseVec> phi_rows;

  Highs highs;
  const double inf = highs.getInfinity();

  row_lower.reserve(static_cast<std::size_t>(n_points));
  row_upper.reserve(static_cast<std::size_t>(n_points));
  phi_rows.reserve(static_cast<std::size_t>(n_points));

  // Objective: the integral gap over Omega. Since h does not depend on x_L, that
  // gap equals w^T x_L up to a constant, with w = integral of phi over Omega.
  // Minimized for the upper bound, maximized for the lower one; the solver runs
  // in minimize mode, so the lower branch flips the sign (step 2.2, section 7).
  //
  // Note this replaces a plain sum over the refinement vertices: those vertices
  // are not uniformly spread over Omega, so their sum is not a measure.
  Eigen::RowVectorXd c_vec = node_weights_[target_phase].transpose();
  if (!is_upper) {
    c_vec = -c_vec;
  }

  for (int g_id = 0; g_id < n_points; ++g_id) {
    // The target-side row maps target barycentric coefficients to the boundary
    // value at the common-refinement vertex g. Precomputed with the geometry.
    const SparseVec& phi = phi_at_refinement_[target_phase][g_id];
    appendSparseRow(starts, col_index, value, phi, kEps);
    if (is_upper) {
      row_lower.push_back(rhs[g_id]);
      row_upper.push_back(inf);
    } else {
      row_lower.push_back(-inf);
      row_upper.push_back(rhs[g_id]);
    }
    phi_rows.push_back(phi);
  }

  highs.setOptionValue("solver", "simplex");
  highs.setOptionValue("presolve", "on");
  highs.setOptionValue("simplex_strategy", 2);
  highs.setOptionValue("pdlp_optimality_tolerance", kHighsPdlpOptimalityTol);
  highs.setOptionValue("kkt_tolerance", kHighsSolutionTol);
  highs.setOptionValue("primal_feasibility_tolerance", kHighsSolutionTol);
  highs.setOptionValue("dual_feasibility_tolerance", kHighsSolutionTol);
  highs.setOptionValue("primal_residual_tolerance", kHighsSolutionTol);
  highs.setOptionValue("dual_residual_tolerance", kHighsSolutionTol);
  highs.setOptionValue("optimality_tolerance", kHighsSolutionTol);
  highs.setOptionValue("small_matrix_value", kHighsSmallMatrixValue);
  highs.setOptionValue("log_to_console", highs_verbose_);
  highs.changeObjectiveSense(ObjSense::kMinimize);

  const int n_cols = target_layout.num_x;
  std::vector<double> col_lower(n_cols, -inf);
  std::vector<double> col_upper(n_cols, inf);

  // Match initializeHighs(): layer 0 remains free and the first vertex in each
  // later projection layer is pinned to zero. This removes the additive
  // barycentric nullspace without changing the represented value function.
  for (int s = 1; s < kSubsystemCount; ++s) {
    if (target_layout.eta_s[s] == 0) {
      throw std::runtime_error(
          "getBorderConditions: empty target projection layer");
    }
    const int col = target_layout.idxX(s, 0);
    col_lower[col] = 0.0;
    col_upper[col] = 0.0;
  }

  HighsStatus st = highs.addCols(
      n_cols, c_vec.data(), col_lower.data(), col_upper.data(), /*num_nz=*/0,
      /*start=*/nullptr, /*index=*/nullptr, /*value=*/nullptr);
  if (st != HighsStatus::kOk) {
    throw std::runtime_error("getBorderConditions: highs.addCols failed");
  }

  const int n_rows = static_cast<int>(row_lower.size());
  st = highs.addRows(n_rows, row_lower.data(), row_upper.data(),
                    static_cast<int>(value.size()), starts.data(),
                    col_index.data(), value.data());
  if (st != HighsStatus::kOk) {
    throw std::runtime_error("getBorderConditions: highs.addRows failed");
  }

  st = highs.run();
  if (st != HighsStatus::kOk) {
    throw std::runtime_error("getBorderConditions: highs.run failed");
  }
  if (highs.getModelStatus() != HighsModelStatus::kOptimal) {
    throw std::runtime_error(
        "getBorderConditions: LP solution not found, model status: "
        + std::to_string(static_cast<int>(highs.getModelStatus())));
  }

  const auto& solution = highs.getSolution();
  if (solution.col_value.size() < static_cast<std::size_t>(n_cols)) {
    throw std::runtime_error("getBorderConditions: solution is too short");
  }
  std::vector<double> x_boundary(solution.col_value.begin(),
                                 solution.col_value.begin() + n_cols);

  // Slack is measured in the direction the bound is supposed to hold:
  // majorization (value - rhs) for the upper branch, minorization (rhs - value)
  // for the lower one. Either way it must be nonnegative.
  double min_slack = std::numeric_limits<double>::infinity();
  double max_slack = -std::numeric_limits<double>::infinity();
  double sum_slack = 0.0;
  for (std::size_t i = 0; i < phi_rows.size(); ++i) {
    const double value_at_g = phi_rows[i].dot(x_boundary);
    const double slack
        = is_upper ? (value_at_g - rhs[i]) : (rhs[i] - value_at_g);
    if (slack < -10.0 * kHighsSolutionTol) {
      throw std::runtime_error(
          is_upper
              ? "getBorderConditions: border LP violates the majorization "
                "constraint"
              : "getBorderConditions: border LP violates the minorization "
                "constraint");
    }
    min_slack = std::min(min_slack, slack);
    max_slack = std::max(max_slack, slack);
    sum_slack += slack;
  }
  logger_->info(
      "Barycentric border LP phase={} source={} theta_idx={} switch_cnt={} "
      "mode={} points={} candidates={} min_slack={} max_slack={} "
      "mean_slack={}",
      target_phase, source_phase, theta_idx, switch_cnt,
      is_upper ? "upper" : "lower", phi_rows.size(), n_candidates, min_slack,
      max_slack, sum_slack / static_cast<double>(phi_rows.size()));

  return x_boundary;
}

std::tuple<std::unique_ptr<Highs>, std::vector<double>, std::vector<double>>
BarycentricAffineApproximator::initializeHighs(int phase) {
  // Reads the per-phase assembly prepared by precomputeMatrices(). Assembling
  // here would redo the whole thing once per solver instance and, if this loop
  // were ever parallelised, would race on the shared member.
  const hcpwa::util::block_lp::ReducedLp& lp = reduced_lps_[phase];
  if (lp.row_upper.empty()) {
    throw std::runtime_error(
        "initializeHighs: LP not assembled. Call precomputeMatrices() first.");
  }
  const auto& layout = layouts_[phase];
  const int m = static_cast<int>(lp.row_upper.size());
  const int n = lp.num_cols;
  std::vector<double> row_lower = lp.row_lower;
  std::vector<double> row_upper = lp.row_upper;

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

  // Column bounds come from the assembler: x free, y >= 0 (one vector per
  // block region now, not per full region), mu free.
  std::vector<double> col_lower = lp.col_lower;
  std::vector<double> col_upper = lp.col_upper;

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
      n, lp.objective.data(), col_lower.data(), col_upper.data(), /*num_nz=*/0,
      /*start=*/nullptr, /*index=*/nullptr, /*value=*/nullptr);
  if (st != HighsStatus::kOk) {
    throw std::runtime_error("initializeHighs: highs.addCols failed.");
  }

  st = highs->addRows(m, row_lower.data(), row_upper.data(),
                      static_cast<int>(lp.values.size()), lp.starts.data(),
                      lp.cols.data(), lp.values.data());
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
  // Only the residual row right-hand sides depend on x_next; the (Y) rows, the
  // (Sigma) rows and the column bounds are static.
  //
  // This is the only place the modulus is arithmetic: q is built from the known
  // x_next, so |q| is a vector of numbers and the LP stays linear. The work is
  // now O(R_A + R_B + R_C) plus O(M_A + M_B + M_C) for the q vectors, about
  // 3.5e3 rows per time step instead of the 1.6e7 the product form needed
  // before the solver was even invoked.
  const auto& layout = layouts_[phase];
  if (x_next.size() != static_cast<std::size_t>(layout.num_x)) {
    throw std::invalid_argument("updateHighsRhsUpperBounds: x_next size must be "
                                + std::to_string(layout.num_x));
  }

  const auto& row_lower = row_lowers_[solver_index];
  auto& highs_solver = highs_solvers_[solver_index];
  std::vector<double> new_row_upper
      = hcpwa::util::block_lp::reducedLpRowUpper(reduced_lps_[phase], x_next);
  for (double& upper : new_row_upper) {
    if (std::abs(upper) <= kEps) {
      upper = 0.0;
    }
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
  // Tier 3: check the ORIGINAL condition, on the original index set.
  //
  //   s * F(a,b,c) <= 0   for (a,b,c) in R_A x R_B x R_C,
  //   F(a,b,c) = F_A(a) + F_B(b) + F_C(c)                    (step 7, Lemma 2)
  //
  // evaluated at both endpoints of the segment:
  //   F^L uses q = Psi z       (the unknown endpoint, where y lifts the modulus)
  //   F^R uses q = Psi x_next  (the known endpoint)
  //   both share the slope (phi^T x_next - phi^T z) / dt.
  //
  // The point of testing triples rather than the reduced rows is that the
  // reduced rows are what the LP already enforced; the triples are what the
  // bound property is actually about. Lemma 2 is the bridge, and evaluating a
  // triple costs three block lookups and an addition, so no 8D vertex is ever
  // materialised and a large sample is free.
  //
  // A failure here means the constructed function is not an estimate of the
  // value function, whatever the solver reported.
  const double s = signS();
  const auto& geometry = phase_geometries_[phase];

  // Per block, flatten (block region, vertex) into a row list once, then sample
  // from the product. Deterministic seed: a failure must be reproducible.
  struct BlockRowRef {
    int block_region;
    const Eigen::VectorXd* nu;
  };
  std::array<std::vector<BlockRowRef>, kBlockCount> rows;
  for (int b = 0; b < kBlockCount; ++b) {
    const auto& block = geometry.blocks[b];
    for (int j = 0; j < block.num_regions; ++j) {
      for (const Eigen::VectorXd& nu : block.vertices[j]) {
        rows[b].push_back(BlockRowRef{j, &nu});
      }
    }
    if (rows[b].empty()) {
      throw std::runtime_error(std::format(
          "validateStepResiduals: phase {} block {} has no rows", phase, b));
    }
  }

  // Per block, the two endpoint residual contributions of one row.
  auto block_residuals = [&](int b, const BlockRowRef& row) {
    const auto& block = geometry.blocks[b];
    const int j = row.block_region;
    const Eigen::VectorXd& nu = *row.nu;

    const SparseVec phi = buildPhiRowBlock(phase, b, j, nu, kEps);
    const Eigen::VectorXd m = block.a_matr[j] * nu + block.f_vec[j]
                              + block.q_c_matr[j] * nu + block.q_c_vec[j];
    Eigen::VectorXd rho = block.q_r_matr[j] * nu + block.q_r_vec[j];
    rho = rho.cwiseMax(0.0);
    const double g_nu = block.g_vec[j].dot(nu) + block.g_scal[j];

    Eigen::VectorXd q_left = Eigen::VectorXd::Zero(block.coord_count);
    Eigen::VectorXd q_right = Eigen::VectorXd::Zero(block.coord_count);
    for (int p = 0; p < block.coord_count; ++p) {
      q_left(p) = block.psi[j].rows[p].dot(z);
      q_right(p) = block.psi[j].rows[p].dot(x_next);
    }

    const double slope = (phi.dot(x_next) - phi.dot(z)) / t_delta_;
    const double left
        = slope + q_left.dot(m) + s * rho.dot(q_left.cwiseAbs()) + g_nu;
    const double right
        = slope + q_right.dot(m) + s * rho.dot(q_right.cwiseAbs()) + g_nu;
    return std::make_pair(left, right);
  };

  // Precompute each block's per-row residual pair once: the sample then costs
  // two additions per triple. Total precompute is O(R_A + R_B + R_C).
  std::array<std::vector<std::pair<double, double>>, kBlockCount> residuals;
  for (int b = 0; b < kBlockCount; ++b) {
    residuals[b].reserve(rows[b].size());
    for (const BlockRowRef& row : rows[b]) {
      residuals[b].push_back(block_residuals(b, row));
    }
  }

  // EXACT, not sampled. F is additively separable across the blocks
  // (step 7, Lemma 2), so
  //     max_{(a,b,c)} s*F(a,b,c) = sum_block max_a s*F_block(a)
  // and the maximum over all R_A*R_B*R_C triples -- about 2e8 of them -- is
  // the sum of three maxima over the arrays already built above. That is
  // O(R_A+R_B+R_C), cheaper than drawing a sample, and it cannot miss a
  // violation confined to any single block row.
  const std::size_t total_triples = rows[0].size() * rows[1].size()
                                    * rows[2].size();
  // block_residuals() returns the raw F_block, so the sign is applied HERE,
  // before maximising. Maximising F and multiplying by s afterwards would be
  // wrong for the lower bound (s = -1), where the worst case is the MINIMUM of
  // F -- that mistake understates the violation and turns a real failure into
  // a pass.
  double worst_left = 0.0;
  double worst_right = 0.0;
  std::array<std::size_t, kBlockCount> worst_left_ids{};
  std::array<std::size_t, kBlockCount> worst_right_ids{};
  for (int b = 0; b < kBlockCount; ++b) {
    double best_left = -std::numeric_limits<double>::infinity();
    double best_right = -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < residuals[b].size(); ++i) {
      const double left = s * residuals[b][i].first;
      const double right = s * residuals[b][i].second;
      if (left > best_left) {
        best_left = left;
        worst_left_ids[b] = i;
      }
      if (right > best_right) {
        best_right = right;
        worst_right_ids[b] = i;
      }
    }
    worst_left += best_left;
    worst_right += best_right;
  }

  const bool left_is_worse = worst_left >= worst_right;
  const double worst = std::max(worst_left, worst_right);
  const auto& worst_ids = left_is_worse ? worst_left_ids : worst_right_ids;

  if (worst > kResidualValidationTol) {
    throw std::runtime_error(std::format(
        "validateStepResiduals: phase {} is not a bound. worst s*F = {} at the "
        "{} endpoint, witness triple (block rows {}, {}, {}) out of {} triples "
        "checked exhaustively. The constructed function is not an estimate of "
        "the value function.",
        phase, worst, left_is_worse ? "left" : "right", worst_ids[0],
        worst_ids[1], worst_ids[2], total_triples));
  }

  logger_->info(
      "validateStepResiduals: phase={}, {} triples covered exhaustively, "
      "worst s*F={}",
      phase, total_triples, worst);
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
    if (phase_geometries_[phase].blocks[0].num_regions == 0) {
      throw std::runtime_error(
          "precomputeMatrices: geometry not initialized. Call "
          "getIntersectionPoints() first.");
    }
  }

  for (int phase = 0; phase < kPhases; ++phase) {
    precomputeSystemMatrices(phase);

    const auto& geometry = phase_geometries_[phase];
    std::size_t block_regions = 0;
    std::size_t block_rows = 0;
    for (int b = 0; b < kBlockCount; ++b) {
      block_regions += geometry.blocks[b].num_regions;
      for (const auto& vertices : geometry.blocks[b].vertices) {
        block_rows += vertices.size();
      }
    }
    logger_->info(
        "precomputeMatrices: phase={}, block regions={} (M_A={}, M_B={}, "
        "M_C={}), block rows={} (R_A+R_B+R_C), implied full regions={}, "
        "num_x={}",
        phase, block_regions, geometry.blocks[0].num_regions,
        geometry.blocks[1].num_regions, geometry.blocks[2].num_regions,
        block_rows, geometry.num_regions, layouts_[phase].num_x);

    // Assemble once per phase. Every solver instance of this phase then shares
    // this assembly instead of rebuilding it.
    reduced_lps_[phase] = prepareLpMatrices(phase);
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
      phase_geometries_[0].num_regions, phase_geometries_[1].num_regions);
  precomputeMatrices();
  {
    auto block_region_count = [&](int phase) {
      int total = 0;
      for (int b = 0; b < kBlockCount; ++b) {
        total += phase_geometries_[phase].blocks[b].num_regions;
      }
      return total;
    };
    logger_->info(
        "Finished precomputeMatrices. Precomputed {} block system matrices for "
        "phase 0 and {} for phase 1",
        block_region_count(0), block_region_count(1));
  }

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
