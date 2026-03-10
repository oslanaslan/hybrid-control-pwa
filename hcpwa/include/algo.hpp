#pragma once

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

// Structure to hold phase intersection computation results
struct PhaseIntersectionResult {
  std::vector<std::vector<size_t>> intersection_prism_indices_phase0;
  std::vector<std::vector<hcpwa::Vec<8>>> intersection_points_phase0;
  std::vector<std::vector<size_t>> intersection_prism_indices_phase1;
  std::vector<std::vector<hcpwa::Vec<8>>> intersection_points_phase1;
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
                                           bool verbose = false);

PolygonAreasVerticesResult compute_polygon_areas_vertices(
    double N, double F, double v, double w, double b51, double b57, double b84,
    double b86, double b31, double b36, double b24, double b27, double f2min,
    double f3min, double f5min, double f8min, double f2max, double f3max,
    double f5max, double f8max, bool verbose = false);

}  // namespace hcpwa