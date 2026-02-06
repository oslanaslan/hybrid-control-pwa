#pragma once

#include <span>
#include <types.hpp>
#include "uniqie_pool.hpp"
#include <vector>

namespace hcpwa {
// Structure to hold computation results for areas vertices
struct AreasVerticesResult {
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

std::vector<hcpwa::PolygonResolution> SplitAABBWithLines(
    AABB<2> aabb, const LineSet<2>& lines);

void NormalizeVertices(std::vector<hcpwa::PolygonResolution>& data,
                       hcpwa::UniquePool<hcpwa::Vec<2>>& pool);

std::vector<std::pair<hcpwa::Triangle, std::size_t>> Triangulate(
    const std::span<hcpwa::PolygonResolution>& data);

hcpwa::LineSet<8> CalcPrism(const Triangle& triangle,
                            const std::array<int, 2>& dims);

template<int N=8>
std::vector<hcpwa::Vec<N>> LinesToPoints(const hcpwa::LineSet<N>& data);

std::vector<hcpwa::TriangleWithUniqueVertices> GetTrianglesWithUniqueVertices(const AABB<2>& aabb, std::vector<hcpwa::PolygonResolution>& polygons);

// Compute areas vertices - pure C++ computation logic
AreasVerticesResult compute_areas_vertices(double N, double F, double v, double w,
                                           double b51, double b57, double b84, double b86,
                                           double b31, double b36, double b24, double b27,
                                           double f2min, double f3min, double f5min, double f8min,
                                           double f2max, double f3max, double f5max, double f8max,
                                           bool verbose = false);

}  // namespace hcpwa