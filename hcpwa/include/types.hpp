#pragma once

#include "float.hpp"
#include <linalg.h>
#include <bitset>

namespace hcpwa {

template <int N>
using Vec = linalg::vec<hcpwa::Float, N>;

constexpr Float modulo(const Vec<2>& v) {  // NOLINT
  return v.x * v.x + v.y * v.y;
}

// Format x*x + y*y + z = 0
template <int N>
using Line = Vec<N + 1>;

template <int N>
using LineSet = std::vector<Line<N>>;

using Polygon = std::vector<Vec<2>>;

struct PolygonResolution {
  Polygon polygon;
  std::bitset<64> is_gt;
};

template <int N>
using AABB = std::pair<Vec<N>, Vec<N>>;

struct Triangle {
  Vec<2> a, b, c;
};

}  // namespace hcpwa