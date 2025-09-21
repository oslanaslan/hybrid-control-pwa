#pragma once

#include <span>
#include <types.hpp>
#include "uniqie_pool.hpp"

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

std::vector<std::pair<hcpwa::Triangle, std::size_t>> TrianglesWithUniqueVertices(const AABB<2>& aabb, std::vector<hcpwa::PolygonResolution>& polygons);

}  // namespace hcpwa