#pragma once

#include "float.hpp"
#include <linalg.h>
#include <bitset>
#include <utility>

namespace linalg {
template <>
struct vec<hcpwa::Float, 8> {
  hcpwa::Float x, y, z, w, x5, x6, x7, x8;
  constexpr hcpwa::Float& operator[](int i) {
    switch (i) {
      case 0:
        return x;
      case 1:
        return y;
      case 2:
        return z;
      case 3:
        return w;
      case 4:
        return x5;
      case 5:
        return x6;
      case 6:
        return x7;
      case 7:
        return x8;
    }
    std::unreachable();
  }

  constexpr const hcpwa::Float& operator[](int i) const {
    switch (i) {
      case 0:
        return x;
      case 1:
        return y;
      case 2:
        return z;
      case 3:
        return w;
      case 4:
        return x5;
      case 5:
        return x6;
      case 6:
        return x7;
      case 7:
        return x8;
    }
    std::unreachable();
  }
};

template <>
struct vec<hcpwa::Float, 9> {
  hcpwa::Float x, y, z, w, x5, x6, x7, x8, x9;
  constexpr hcpwa::Float& operator[](int i) {
    switch (i) {
      case 0:
        return x;
      case 1:
        return y;
      case 2:
        return z;
      case 3:
        return w;
      case 4:
        return x5;
      case 5:
        return x6;
      case 6:
        return x7;
      case 7:
        return x8;
      case 8:
        return x9;
    }
    std::unreachable();
  }

  constexpr const hcpwa::Float& operator[](int i) const {
    switch (i) {
      case 0:
        return x;
      case 1:
        return y;
      case 2:
        return z;
      case 3:
        return w;
      case 4:
        return x5;
      case 5:
        return x6;
      case 6:
        return x7;
      case 7:
        return x8;
      case 8:
        return x9;
    }
    std::unreachable();
  }
};
}  // namespace linalg

namespace hcpwa {

template <int N>
  requires(sizeof(linalg::vec<Float, N>) == sizeof(Float) * N
           && alignof(linalg::vec<Float, N>) == alignof(Float))
using Vec = linalg::vec<hcpwa::Float, N>;

constexpr Float modulo(const Vec<2>& v) {  // NOLINT
  return v.x * v.x + v.y * v.y;
}

// Format x*x + y*y + z = 0
template <int N>
using Line = Vec<N + 1>;

template <int N>
constexpr Float Apply(const Line<N>& line, const Vec<N>& vec) {
  Float result = line[N];
  for (int i = 0; i < N; i++) {
    result += line[i] * vec[i];
  }
  return result;
}

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