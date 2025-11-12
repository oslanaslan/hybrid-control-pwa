#pragma once

#include <span>
#include <types.hpp>
#include "uniqie_pool.hpp"
#include <vector>

// Structure to hold computation results for areas vertices
struct AreasVerticesResult {
  std::vector<hcpwa::TriangleWithUniqueVertices> triangles51;
  std::vector<hcpwa::TriangleWithUniqueVertices> triangles57;
  std::vector<hcpwa::TriangleWithUniqueVertices> triangles84;
  std::vector<hcpwa::TriangleWithUniqueVertices> triangles86;
  std::vector<std::vector<hcpwa::Vec<8>>> intersection_points;
  std::vector<std::vector<size_t>> intersection_prism_indices;
};

namespace hcpwa {
std::vector<hcpwa::PolygonResolution> SplitAABBWithLines(
    AABB<2> aabb, const LineSet<2>& lines);

void NormalizeVertices(std::vector<hcpwa::PolygonResolution>& data,
                       hcpwa::UniquePool<hcpwa::Vec<2>>& pool);

std::vector<std::pair<hcpwa::Triangle, std::size_t>> Triangulate(
    const std::span<hcpwa::PolygonResolution>& data);

hcpwa::LineSet<8> CalcPrism(const Triangle& triangle,
                            const std::array<int, 2>& dims);

std::vector<hcpwa::Vec<8>> LinesToPoints(const hcpwa::LineSet<8>& data);

std::vector<hcpwa::TriangleWithUniqueVertices> GetTrianglesWithUniqueVertices(const AABB<2>& aabb, std::vector<hcpwa::PolygonResolution>& polygons);

// Compute areas vertices - pure C++ computation logic
AreasVerticesResult compute_areas_vertices(int phase, float N, float F, float v, float w,
                                           float b51, float b57, float b84, float b86,
                                           float b31, float b36, float b24, float b27,
                                           float f2min, float f3min, float f5min, float f8min,
                                           float f2max, float f3max, float f5max, float f8max);

}  // namespace hcpwa