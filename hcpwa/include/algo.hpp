#pragma once

#include <types.hpp>

namespace hcpwa {
std::vector<PolygonResolution> SplitAABBWithLines(AABB<2> aabb,
                                                  const LineSet<2>& lines);
}