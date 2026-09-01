#include "barycentric_geometry_types.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <format>
#include <numeric>
#include <stdexcept>
#include <vector>

// NOLINTBEGIN(readability-identifier-naming)

namespace barycentric_affine_approximator {

namespace {

// Shared body of the two phi builders: the barycentric weights of one already
// projected 2D point on one triangle of one layer, written to global columns.
void addLayerPhi(const ProjectionLayer& layer,
                 const BarycentricVarLayout& layout, int layer_id,
                 int triangle_id, const Eigen::Vector2d& projected,
                 double tolerance, const char* who, SparseVec* phi) {
  if (triangle_id < 0 || triangle_id >= static_cast<int>(layer.bases.size())) {
    throw std::runtime_error(std::string(who) + ": triangle id out of range");
  }
  const TriangleBasis& basis = layer.bases[static_cast<std::size_t>(
      triangle_id)];
  const Eigen::Vector3d alpha = basis.H * projected + basis.h;

  // A value row is only meaningful in the selected local simplex. The tolerance
  // allows points on shared triangle/area boundaries, but rejects genuine tuple
  // order mistakes that would project outside the simplex. For the block form
  // this check does double duty: almost every wrong local_axis pair projects
  // the point somewhere else on the plane and lands outside.
  if (std::abs(alpha.sum() - 1.0) > tolerance) {
    throw std::runtime_error(std::string(who)
                             + ": barycentric alpha sum is not one");
  }
  for (int local_vertex = 0; local_vertex < 3; ++local_vertex) {
    if (alpha(local_vertex) < -tolerance
        || alpha(local_vertex) > 1.0 + tolerance) {
      throw std::runtime_error(std::string(who)
                               + ": point is outside selected simplex");
    }
    const int x_col = layout.idxX(layer_id, basis.vertex_ids[
        static_cast<std::size_t>(local_vertex)]);
    // kGeomEps, not the SparseVec default: kEps (1e-5) is larger than the
    // geometry tolerance and can remove enough mass for the five-layer phi row
    // to stop summing to five.
    phi->add(x_col, alpha(local_vertex), kGeomEps);
  }
}

}  // namespace

SparseVec buildPhiRow(const PhaseGeometry& geometry,
                      const BarycentricVarLayout& layout, int region,
                      const Eigen::VectorXd& point, double tolerance) {
  if (point.size() != kSpaceDim) {
    throw std::invalid_argument("buildPhiRow: point must be 8-dimensional");
  }
  if (region < 0
      || region >= static_cast<int>(geometry.region_triangle_ids.size())) {
    throw std::invalid_argument("buildPhiRow: invalid region id");
  }

  SparseVec phi;
  const auto& triangle_ids
      = geometry.region_triangle_ids[static_cast<std::size_t>(region)];
  for (int s = 0; s < kSubsystemCount; ++s) {
    const ProjectionLayer& layer = geometry.layers[static_cast<std::size_t>(s)];
    const Eigen::Vector2d projected(point(layer.axes[0]), point(layer.axes[1]));
    addLayerPhi(layer, layout, s, triangle_ids[static_cast<std::size_t>(s)],
                projected, tolerance, "buildPhiRow", &phi);
  }
  return phi;
}

SparseVec buildPhiRowBlock(const PhaseGeometry& geometry,
                           const BarycentricVarLayout& layout, int block,
                           int block_region, const Eigen::VectorXd& point,
                           double tolerance) {
  if (block < 0 || block >= kBlockCount) {
    throw std::invalid_argument("buildPhiRowBlock: invalid block");
  }
  const BlockGeometry& block_geometry
      = geometry.blocks[static_cast<std::size_t>(block)];
  if (point.size() != block_geometry.coord_count) {
    throw std::invalid_argument(
        "buildPhiRowBlock: point must have the block's coordinate count");
  }
  if (block_region < 0 || block_region >= block_geometry.numRegions()) {
    throw std::invalid_argument("buildPhiRowBlock: invalid block region id");
  }

  SparseVec phi;
  const auto& triangle_ids
      = block_geometry.triangle_ids[static_cast<std::size_t>(block_region)];
  for (int l = 0; l < block_geometry.layer_count; ++l) {
    const int s = block_geometry.layer_ids[static_cast<std::size_t>(l)];
    const ProjectionLayer& layer = geometry.layers[static_cast<std::size_t>(s)];
    const auto& local = block_geometry.local_axis[static_cast<std::size_t>(l)];
    const Eigen::Vector2d projected(point(local[0]), point(local[1]));
    addLayerPhi(layer, layout, s, triangle_ids[static_cast<std::size_t>(l)],
                projected, tolerance, "buildPhiRowBlock", &phi);
  }
  return phi;
}

std::array<BlockGeometry, kBlockCount> blockGeometryFromRegions(
    int phase, const std::array<hcpwa::BlockRegions, kBlockCount>& src) {
  // Ingests one phase's block factorisation and derives, for every projection
  // plane of a block, where its two axes sit inside that block's coordinate
  // list. Everything here is checked rather than assumed: a wrong local axis
  // pair produces a silently wrong phi row, never a crash.
  const auto expected_blocks = coordinateBlocksForPhase(phase);

  std::array<BlockGeometry, kBlockCount> dst;
  std::vector<int> covered_coords;
  std::vector<int> covered_layers;

  for (int b = 0; b < kBlockCount; ++b) {
    const hcpwa::BlockRegions& source = src[static_cast<std::size_t>(b)];
    BlockGeometry block;
    block.coords = source.coords;
    block.coord_count = source.coord_count;
    block.layer_ids = source.layer_ids;
    block.layer_count = source.layer_count;
    if (block.coord_count < 2 || block.layer_count < 1) {
      throw std::runtime_error("blockGeometryFromRegions: degenerate block");
    }

    block.local_axis = localAxesForBlock(phase, block.coords,
                                         block.coord_count, block.layer_ids,
                                         block.layer_count);
    for (int l = 0; l < block.layer_count; ++l) {
      covered_layers.push_back(block.layer_ids[static_cast<std::size_t>(l)]);
    }
    for (int c = 0; c < block.coord_count; ++c) {
      covered_coords.push_back(block.coords[static_cast<std::size_t>(c)]);
    }

    block.triangle_ids.reserve(source.triangle_ids.size());
    for (const auto& ids : source.triangle_ids) {
      std::array<int, kMaxBlockLayers> converted{};
      for (int l = 0; l < block.layer_count; ++l) {
        converted[static_cast<std::size_t>(l)]
            = ids[static_cast<std::size_t>(l)];
      }
      block.triangle_ids.push_back(converted);
    }

    block.vertices.resize(source.vertices.size());
    block.aabb_lower.reserve(source.aabb.size());
    block.aabb_upper.reserve(source.aabb.size());
    for (std::size_t j = 0; j < source.vertices.size(); ++j) {
      block.vertices[j].reserve(source.vertices[j].size());
      for (const auto& vertex : source.vertices[j]) {
        Eigen::VectorXd converted(block.coord_count);
        for (int c = 0; c < block.coord_count; ++c) {
          converted(c) = vertex[static_cast<std::size_t>(c)];
        }
        block.vertices[j].push_back(std::move(converted));
      }
      Eigen::VectorXd lower(block.coord_count);
      Eigen::VectorXd upper(block.coord_count);
      for (int c = 0; c < block.coord_count; ++c) {
        lower(c) = source.aabb[j].lower[static_cast<std::size_t>(c)];
        upper(c) = source.aabb[j].upper[static_cast<std::size_t>(c)];
      }
      block.aabb_lower.push_back(std::move(lower));
      block.aabb_upper.push_back(std::move(upper));
    }
    if (block.triangle_ids.size() != block.vertices.size()) {
      throw std::runtime_error(
          "blockGeometryFromRegions: triangle ids and vertices disagree");
    }
    dst[static_cast<std::size_t>(b)] = std::move(block);
  }

  // The blocks must partition the eight coordinates and the five planes, and
  // they must be the components of the plane incidence graph -- which is
  // derived from projectionAxesForPhase alone and so cannot drift with the
  // geometry code.
  std::sort(covered_coords.begin(), covered_coords.end());
  std::sort(covered_layers.begin(), covered_layers.end());
  std::vector<int> all_coords(kSpaceDim);
  std::iota(all_coords.begin(), all_coords.end(), 0);
  std::vector<int> all_layers(kSubsystemCount);
  std::iota(all_layers.begin(), all_layers.end(), 0);
  if (covered_coords != all_coords || covered_layers != all_layers) {
    throw std::runtime_error(
        "blockGeometryFromRegions: blocks do not partition coordinates and "
        "planes");
  }
  for (int b = 0; b < kBlockCount; ++b) {
    const BlockGeometry& block = dst[static_cast<std::size_t>(b)];
    std::vector<int> coords(block.coords.begin(),
                            block.coords.begin() + block.coord_count);
    std::sort(coords.begin(), coords.end());
    bool found = false;
    for (const auto& expected : expected_blocks) {
      if (expected == coords) {
        found = true;
        break;
      }
    }
    if (!found) {
      throw std::runtime_error(std::format(
          "blockGeometryFromRegions: phase {} block {} is not a connected "
          "component of the plane incidence graph",
          phase, b));
    }
  }
  return dst;
}

std::array<std::array<int, 2>, kMaxBlockLayers> localAxesForBlock(
    int phase, const std::array<int, kMaxBlockCoords>& coords, int coord_count,
    const std::array<int, kMaxBlockLayers>& layer_ids, int layer_count) {
  const auto axes = projectionAxesForPhase(phase);
  std::array<std::array<int, 2>, kMaxBlockLayers> local{};
  for (int l = 0; l < layer_count; ++l) {
    const int layer = layer_ids[static_cast<std::size_t>(l)];
    if (layer < 0 || layer >= kSubsystemCount) {
      throw std::invalid_argument("localAxesForBlock: bad layer id");
    }
    for (int k = 0; k < 2; ++k) {
      const int axis = axes[static_cast<std::size_t>(layer)][
          static_cast<std::size_t>(k)];
      int position = -1;
      for (int c = 0; c < coord_count; ++c) {
        if (coords[static_cast<std::size_t>(c)] == axis) {
          position = c;
          break;
        }
      }
      if (position < 0) {
        // Both axes of a plane living in one block is what makes the region
        // projections exact rectangles and the residual separable. If this
        // ever fires, the reduction does not apply and nothing below it does
        // either.
        throw std::runtime_error(std::format(
            "localAxesForBlock: phase {} layer {} axis {} is not in its own "
            "block", phase, layer, axis));
      }
      local[static_cast<std::size_t>(l)][static_cast<std::size_t>(k)]
          = position;
    }
    if (local[static_cast<std::size_t>(l)][0]
        == local[static_cast<std::size_t>(l)][1]) {
      throw std::runtime_error(
          "localAxesForBlock: a plane maps both axes to one coordinate");
    }
  }
  return local;
}

}  // namespace barycentric_affine_approximator

// NOLINTEND(readability-identifier-naming)
