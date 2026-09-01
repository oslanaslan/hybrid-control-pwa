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

// Axis-aligned bounding box of one 8D area. Used by the common refinement to
// prune candidate cell pairs. It lives here rather than inside
// compute_common_refinement_area_vertices() because callers now compose these
// boxes from per-block boxes and pass them in.
struct AreaBounds8d {
  std::array<double, 8> min{};
  std::array<double, 8> max{};
};

// Axis-aligned bounding box of one block region. Only the first coord_count
// entries of the owning BlockRegions are meaningful.
struct BlockBounds3d {
  std::array<double, 3> min{};
  std::array<double, 3> max{};
};

// One of the three coordinate blocks of one phase.
//
// The five projection planes of a phase split into three groups that share no
// coordinate, so an 8D area is the Cartesian product of three low-dimensional
// block polytopes and its vertex set is the product of their vertex sets
// (step7_block_reduction.md, Lemma 1). Phase 0 groups the planes as
// {31,36} x {24,27} x {58} over coordinates {0,2,5} x {1,3,6} x {4,7}; phase 1
// as {51,57} x {84,86} x {23} over {0,4,6} x {3,5,7} x {1,2}.
//
// Block C is two-dimensional: coord_count is 2, layer_count is 1, and only the
// first two components of each vertex and of the bounds are meaningful.
//
// REGION ID ORDER. The 8D area ids that the rest of the pipeline uses are
// assigned in the order the assembly loops run, block A outermost:
//   area_id = (j_A * M_B + j_B) * M_C + j_C.
// Anything that maps between an area id and a block triple must use exactly
// this order.
struct BlockRegions {
  // Unused trailing slots are -1, not 0, because 0 is a valid coordinate and a
  // valid layer id: reading past coord_count or layer_count must be obviously
  // wrong rather than quietly plausible.
  std::array<int, 3> coords{};     // state coordinates, ascending
  int coord_count = 0;             // 3 for blocks A and B, 2 for block C
  std::array<int, 2> layer_ids{};  // indices into the phase's five layers
  int layer_count = 0;             // 2 for blocks A and B, 1 for block C
  // Per block region, the simplex id in each of the block's layers. For block
  // C only entry [0] is meaningful.
  std::vector<std::array<size_t, 2>> triangle_ids;
  // Per block region, ALL vertices of the block polytope. This is the data the
  // 8D assembly loops truncate to the first two.
  std::vector<std::vector<hcpwa::Vec<3>>> vertices;
  std::vector<BlockBounds3d> bounds;
};

// Structure to hold phase intersection computation results
struct PhaseIntersectionResult {
  std::vector<std::vector<size_t>> intersection_prism_indices_phase0;
  std::vector<std::vector<hcpwa::Vec<8>>> intersection_points_phase0;
  std::vector<std::vector<size_t>> intersection_prism_indices_phase1;
  std::vector<std::vector<hcpwa::Vec<8>>> intersection_points_phase1;
  // Blocks A, B, C in that order. Always populated, including when the 8D
  // vertex lists above are switched off.
  std::array<BlockRegions, 3> blocks_phase0;
  std::array<BlockRegions, 3> blocks_phase1;
  // Per 8D area, the box composed by concatenating the three block boxes.
  // Exact, because the area is a product and is unconstrained in the other
  // blocks' coordinates.
  std::vector<AreaBounds8d> area_bounds_phase0;
  std::vector<AreaBounds8d> area_bounds_phase1;
};

// Options for the triangle geometry pipeline.
//
// build_8d_regions materialises intersection_points_phase{0,1}, the explicit
// 8D vertex list per area. Those loops take their bounds from the 2-element
// prism-index list while indexing the block vertex lists, so they keep only
// the first two vertices of each 3D block and produce 12 of the 108-192
// vertices an area actually has. See
// docs/barycentric_block_reduction_context.md part I.
//
// The barycentric path switches them off and works from the block data
// instead. The piecewise path still consumes them and keeps the default; the
// truncation is left in place there deliberately, because correcting it would
// cost about 13 GB of vertices for a code path that is slated for deletion
// (same document, part IV.2).
struct TriangleAreasOptions {
  bool verbose = false;
  bool build_8d_regions = true;
};

// One nonempty cell in the common refinement of the phase-0 and phase-1 area
// partitions. The area ids index the already computed per-phase area arrays,
// while the prism index arrays keep the phase-specific tuple order explicit:
//   phase 0: [31, 36, 24, 27, 58]
//   phase 1: [51, 57, 84, 86, 23]
struct CommonRefinementArea {
  size_t phase0_area_id = 0;
  size_t phase1_area_id = 0;
  std::array<size_t, 5> phase0_prism_indices{};
  std::array<size_t, 5> phase1_prism_indices{};
  std::vector<hcpwa::Vec<8>> vertices;
};

struct CommonRefinementResult {
  std::vector<CommonRefinementArea> areas;
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
  CommonRefinementResult common_refinement;
  // Block decomposition of each phase, blocks A, B, C in that order. Populated
  // whether or not the 8D vertex lists above were built.
  std::array<BlockRegions, 3> blocks_phase0;
  std::array<BlockRegions, 3> blocks_phase1;
};

struct PolygonAreasVerticesResult {
  // Phase 0
  std::vector<std::vector<hcpwa::Vec<8>>> intersection_points_phase0;
  std::vector<std::vector<size_t>> intersection_prism_indices_phase0;
  // Phase 1
  std::vector<std::vector<hcpwa::Vec<8>>> intersection_points_phase1;
  std::vector<std::vector<size_t>> intersection_prism_indices_phase1;
  // Block decomposition of each phase, blocks A, B, C in that order. Populated
  // whether or not the 8D vertex lists above were built. Block C here is a
  // general convex polygon cell rather than a triangle, so its vertex count is
  // only bounded below by 3.
  std::array<BlockRegions, 3> blocks_phase0;
  std::array<BlockRegions, 3> blocks_phase1;
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
    TriangleAreasOptions options = {});

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
    TriangleAreasOptions options = {});

CommonRefinementResult compute_common_refinement_area_vertices(
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
    const std::vector<std::vector<size_t>>& phase0_area_prism_indices,
    const std::vector<std::vector<size_t>>& phase1_area_prism_indices,
    // Per-area bounding boxes. These used to be derived here from the 8D area
    // vertex lists, but those lists are truncated to 12 of an area's 108-192
    // vertices, so the boxes came out too small and often degenerate. A
    // degenerate box makes has_positive_width() discard the area outright and
    // an undersized one makes the overlap test miss pairs, which silently
    // shrinks the common refinement and hence the border LP's test set. The
    // caller now composes these boxes by concatenating per-block boxes built
    // from complete block vertex sets: exact, cheaper, and correct. See
    // docs/barycentric_block_reduction_context.md parts I.8 and IV.2.
    const std::vector<AreaBounds8d>& phase0_area_bounds,
    const std::vector<AreaBounds8d>& phase1_area_bounds,
    hcpwa::Float N, bool verbose = false);

// Compute areas vertices - pure C++ computation logic
TriangleAreasVerticesResult compute_triangle_areas_vertices(double N, double F, double v, double w,
                                           double b51, double b57, double b84, double b86,
                                           double b31, double b36, double b24, double b27,
                                           double f2min, double f3min, double f5min, double f8min,
                                           double f2max, double f3max, double f5max, double f8max,
                                           TriangleAreasOptions options = {});

PolygonAreasVerticesResult compute_polygon_areas_vertices(
    double N, double F, double v, double w, double b51, double b57, double b84,
    double b86, double b31, double b36, double b24, double b27, double f2min,
    double f3min, double f5min, double f8min, double f2max, double f3max,
    double f5max, double f8max, TriangleAreasOptions options = {});

}  // namespace hcpwa