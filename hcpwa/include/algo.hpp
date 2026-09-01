#pragma once

#include <array>
#include <cstddef>
#include <span>
#include <types.hpp>
#include "uniqie_pool.hpp"
#include <vector>

namespace hcpwa {
// Structure to hold polygon resolution results for each 2D plane
struct PolygonResolutions {
  std::vector<hcpwa::PolygonResolution> resolution_51;
  std::vector<hcpwa::PolygonResolution> resolution_57;
  std::vector<hcpwa::PolygonResolution> resolution_84;
  std::vector<hcpwa::PolygonResolution> resolution_86;
  std::vector<hcpwa::PolygonResolution> resolution_58;
  std::vector<hcpwa::PolygonResolution> resolution_31;
  std::vector<hcpwa::PolygonResolution> resolution_36;
  std::vector<hcpwa::PolygonResolution> resolution_24;
  std::vector<hcpwa::PolygonResolution> resolution_27;
  std::vector<hcpwa::PolygonResolution> resolution_23;
};

// Structure to hold triangulation and prism computation results
struct TriangulationAndPrismsResult {
  // Prisms (Phase 0 and Phase 1)
  std::vector<hcpwa::LineSet<8>> prisms31;
  std::vector<hcpwa::LineSet<8>> prisms36;
  std::vector<hcpwa::LineSet<8>> prisms24;
  std::vector<hcpwa::LineSet<8>> prisms27;
  std::vector<hcpwa::LineSet<8>> prisms58;
  std::vector<hcpwa::LineSet<8>> prisms51;
  std::vector<hcpwa::LineSet<8>> prisms57;
  std::vector<hcpwa::LineSet<8>> prisms84;
  std::vector<hcpwa::LineSet<8>> prisms86;
  std::vector<hcpwa::LineSet<8>> prisms23;
  // Triangles (Phase 0 and Phase 1)
  std::vector<hcpwa::TriangleWithUniqueVertices> triangles31;
  std::vector<hcpwa::TriangleWithUniqueVertices> triangles36;
  std::vector<hcpwa::TriangleWithUniqueVertices> triangles24;
  std::vector<hcpwa::TriangleWithUniqueVertices> triangles27;
  std::vector<hcpwa::TriangleWithUniqueVertices> triangles58;
  std::vector<hcpwa::TriangleWithUniqueVertices> triangles51;
  std::vector<hcpwa::TriangleWithUniqueVertices> triangles57;
  std::vector<hcpwa::TriangleWithUniqueVertices> triangles84;
  std::vector<hcpwa::TriangleWithUniqueVertices> triangles86;
  std::vector<hcpwa::TriangleWithUniqueVertices> triangles23;
};

// Structure to hold polygon prism computation results
struct PolygonPrismsResult {
  // Prisms (Phase 0 and Phase 1)
  std::vector<hcpwa::LineSet<8>> prisms31;
  std::vector<hcpwa::LineSet<8>> prisms36;
  std::vector<hcpwa::LineSet<8>> prisms24;
  std::vector<hcpwa::LineSet<8>> prisms27;
  std::vector<hcpwa::LineSet<8>> prisms58;
  std::vector<hcpwa::LineSet<8>> prisms51;
  std::vector<hcpwa::LineSet<8>> prisms57;
  std::vector<hcpwa::LineSet<8>> prisms84;
  std::vector<hcpwa::LineSet<8>> prisms86;
  std::vector<hcpwa::LineSet<8>> prisms23;
};

// The eight state coordinates of one phase split into three blocks that no
// projection plane straddles, and the 8D region of that phase is the direct
// product of one cell from each block. Phase 0 splits as {0,2,5} / {1,3,6} /
// {4,7}, phase 1 as {0,4,6} / {3,5,7} / {1,2}.
//
// The cells themselves are what compute_intersection_points() already builds
// before it forms the product; these structures only stop them from being
// thrown away.
constexpr int kBlockCount = 3;
constexpr int kMaxBlockCoords = 3;
constexpr int kMaxBlockLayers = 2;

struct BlockAabb {
  std::array<double, kMaxBlockCoords> lower{};
  std::array<double, kMaxBlockCoords> upper{};
};

struct BlockRegions {
  // Global state coordinates of this block, ascending. Only the first
  // coord_count entries are meaningful.
  std::array<int, kMaxBlockCoords> coords{};
  int coord_count = 0;
  // Positions of this block's projection planes in the phase's five-plane
  // tuple. Only the first layer_count entries are meaningful.
  std::array<int, kMaxBlockLayers> layer_ids{};
  int layer_count = 0;

  // Per block-region, one triangle id per plane of this block.
  std::vector<std::array<int, kMaxBlockLayers>> triangle_ids;
  // Per block-region, every vertex in this block's own coordinate order.
  std::vector<std::vector<std::array<double, kMaxBlockCoords>>> vertices;
  // Per block-region, the bounding box of those vertices.
  std::vector<BlockAabb> aabb;

  std::size_t num_regions() const { return vertices.size(); }  // NOLINT
};

// Structure to hold phase intersection computation results
struct PhaseIntersectionResult {
  std::vector<std::vector<size_t>> intersection_prism_indices_phase0;
  std::vector<std::vector<hcpwa::Vec<8>>> intersection_points_phase0;
  std::vector<std::vector<size_t>> intersection_prism_indices_phase1;
  std::vector<std::vector<hcpwa::Vec<8>>> intersection_points_phase1;
  // The block factorisation of the same regions. Always filled, including when
  // the 8D product above is skipped.
  std::array<BlockRegions, kBlockCount> blocks_phase0;
  std::array<BlockRegions, kBlockCount> blocks_phase1;
};

// Options of the triangle geometry path.
struct TriangleGeometryOptions {
  // Materialise intersection_points_phase{0,1}, the full 8D product of the
  // block cells. Nothing reads it any more: the LP and the border solver are
  // both assembled block by block. It costs prod_b |V_b| vertices per region --
  // 90 million on the N=100 arrangement, where the block factors come to 2 639
  // -- so it is off by default and exists for tests that compare the two paths.
  //
  // The prism index lists are always built: they are five integers per region
  // and they carry the triangle ids the barycentric charts are selected by.
  bool build_8d_vertices = false;
};

// Structure to hold computation results for areas vertices
struct TriangleAreasVerticesResult {
  // Phase 0
  std::vector<hcpwa::TriangleWithUniqueVertices> triangles31;
  std::vector<hcpwa::TriangleWithUniqueVertices> triangles36;
  std::vector<hcpwa::TriangleWithUniqueVertices> triangles24;
  std::vector<hcpwa::TriangleWithUniqueVertices> triangles27;
  std::vector<hcpwa::TriangleWithUniqueVertices> triangles58;
  std::vector<std::vector<hcpwa::Vec<8>>> intersection_points_phase0;
  std::vector<std::vector<size_t>> intersection_prism_indices_phase0;
  // Phase 1
  std::vector<hcpwa::TriangleWithUniqueVertices> triangles51;
  std::vector<hcpwa::TriangleWithUniqueVertices> triangles57;
  std::vector<hcpwa::TriangleWithUniqueVertices> triangles84;
  std::vector<hcpwa::TriangleWithUniqueVertices> triangles86;
  std::vector<hcpwa::TriangleWithUniqueVertices> triangles23;
  std::vector<std::vector<hcpwa::Vec<8>>> intersection_points_phase1;
  std::vector<std::vector<size_t>> intersection_prism_indices_phase1;
  std::array<BlockRegions, kBlockCount> blocks_phase0;
  std::array<BlockRegions, kBlockCount> blocks_phase1;
};

struct PolygonAreasVerticesResult {
  // Phase 0
  std::vector<std::vector<hcpwa::Vec<8>>> intersection_points_phase0;
  std::vector<std::vector<size_t>> intersection_prism_indices_phase0;
  // Phase 1
  std::vector<std::vector<hcpwa::Vec<8>>> intersection_points_phase1;
  std::vector<std::vector<size_t>> intersection_prism_indices_phase1;
};

std::vector<hcpwa::PolygonResolution> SplitAABBWithLines(
    AABB<2> aabb, const LineSet<2>& lines);

void NormalizeVertices(std::vector<hcpwa::PolygonResolution>& data,
                       hcpwa::UniquePool<hcpwa::Vec<2>>& pool);

std::vector<std::pair<hcpwa::Triangle, std::size_t>> Triangulate(
    const std::span<hcpwa::PolygonResolution>& data);

hcpwa::LineSet<8> CalcPrism(const Triangle& triangle,
                            const std::array<int, 2>& dims);

hcpwa::LineSet<8> CalcPrism(const Polygon& polygon,
                            const std::array<int, 2>& dims);

template<int N=8>
std::vector<hcpwa::Vec<N>> LinesToPoints(const hcpwa::LineSet<N>& data);

std::vector<hcpwa::TriangleWithUniqueVertices> GetTrianglesWithUniqueVertices(const AABB<2>& aabb, std::vector<hcpwa::PolygonResolution>& polygons);

PolygonResolutions compute_polygon_resolutions(double N, double F, double v, double w,
                                               double b51, double b57, double b84, double b86,
                                               double b31, double b36, double b24, double b27,
                                               double f2min, double f3min, double f5min, double f8min,
                                               double f2max, double f3max, double f5max, double f8max,
                                               bool verbose = false);

TriangulationAndPrismsResult compute_triangulation_and_prisms(
    PolygonResolutions& polygon_resolutions,
    const hcpwa::AABB<2>& aabb2d,
    bool verbose = false);

PolygonPrismsResult compute_prisms_from_polygons(
    PolygonResolutions& polygon_resolutions,
    const hcpwa::AABB<2>& aabb2d,
    bool verbose = false);

PhaseIntersectionResult compute_intersection_points(
    const std::vector<hcpwa::LineSet<8>>& prisms31,
    const std::vector<hcpwa::LineSet<8>>& prisms36,
    const std::vector<hcpwa::LineSet<8>>& prisms24,
    const std::vector<hcpwa::LineSet<8>>& prisms27,
    const std::vector<hcpwa::LineSet<8>>& prisms58,
    const std::vector<hcpwa::LineSet<8>>& prisms51,
    const std::vector<hcpwa::LineSet<8>>& prisms57,
    const std::vector<hcpwa::LineSet<8>>& prisms84,
    const std::vector<hcpwa::LineSet<8>>& prisms86,
    const std::vector<hcpwa::LineSet<8>>& prisms23,
    const std::vector<hcpwa::TriangleWithUniqueVertices>& triangles58,
    const std::vector<hcpwa::TriangleWithUniqueVertices>& triangles23,
    hcpwa::Float N,
    bool verbose = false,
    const TriangleGeometryOptions& options = {});

PhaseIntersectionResult compute_intersection_points(
    const std::vector<hcpwa::LineSet<8>>& prisms31,
    const std::vector<hcpwa::LineSet<8>>& prisms36,
    const std::vector<hcpwa::LineSet<8>>& prisms24,
    const std::vector<hcpwa::LineSet<8>>& prisms27,
    const std::vector<hcpwa::LineSet<8>>& prisms58,
    const std::vector<hcpwa::LineSet<8>>& prisms51,
    const std::vector<hcpwa::LineSet<8>>& prisms57,
    const std::vector<hcpwa::LineSet<8>>& prisms84,
    const std::vector<hcpwa::LineSet<8>>& prisms86,
    const std::vector<hcpwa::LineSet<8>>& prisms23,
    const std::vector<hcpwa::PolygonResolution>& polygons58,
    const std::vector<hcpwa::PolygonResolution>& polygons23,
    hcpwa::Float N,
    bool verbose = false);

// Compute areas vertices - pure C++ computation logic
TriangleAreasVerticesResult compute_triangle_areas_vertices(double N, double F, double v, double w,
                                           double b51, double b57, double b84, double b86,
                                           double b31, double b36, double b24, double b27,
                                           double f2min, double f3min, double f5min, double f8min,
                                           double f2max, double f3max, double f5max, double f8max,
                                           bool verbose = false,
                                           const TriangleGeometryOptions& options = {});

PolygonAreasVerticesResult compute_polygon_areas_vertices(
    double N, double F, double v, double w, double b51, double b57, double b84,
    double b86, double b31, double b36, double b24, double b27, double f2min,
    double f3min, double f5min, double f8min, double f2max, double f3max,
    double f5max, double f8max, bool verbose = false);

}  // namespace hcpwa