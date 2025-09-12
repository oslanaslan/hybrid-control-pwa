#pragma once

#include "float.hpp"
// #include <linalg.h>
#include <bitset>
#include <format>
#include <utility>
#include <vector.hpp>
#include <vector>

// namespace linalg {
// template <>
// struct vec<hcpwa::Float, 8> {
//   hcpwa::Float x, y, z, w, x5, x6, x7, x8;
//   constexpr hcpwa::Float& operator[](int i) {
//     switch (i) {
//       case 0:
//         return x;
//       case 1:
//         return y;
//       case 2:
//         return z;
//       case 3:
//         return w;
//       case 4:
//         return x5;
//       case 5:
//         return x6;
//       case 6:
//         return x7;
//       case 7:
//         return x8;
//     }
//     std::unreachable();
//   }

//   constexpr const hcpwa::Float& operator[](int i) const {
//     switch (i) {
//       case 0:
//         return x;
//       case 1:
//         return y;
//       case 2:
//         return z;
//       case 3:
//         return w;
//       case 4:
//         return x5;
//       case 5:
//         return x6;
//       case 6:
//         return x7;
//       case 7:
//         return x8;
//     }
//     std::unreachable();
//   }
// };

// template <>
// struct vec<hcpwa::Float, 9> {
//   hcpwa::Float x, y, z, w, x5, x6, x7, x8, x9;
//   constexpr hcpwa::Float& operator[](int i) {
//     switch (i) {
//       case 0:
//         return x;
//       case 1:
//         return y;
//       case 2:
//         return z;
//       case 3:
//         return w;
//       case 4:
//         return x5;
//       case 5:
//         return x6;
//       case 6:
//         return x7;
//       case 7:
//         return x8;
//       case 8:
//         return x9;
//     }
//     std::unreachable();
//   }

//   constexpr const hcpwa::Float& operator[](int i) const {
//     switch (i) {
//       case 0:
//         return x;
//       case 1:
//         return y;
//       case 2:
//         return z;
//       case 3:
//         return w;
//       case 4:
//         return x5;
//       case 5:
//         return x6;
//       case 6:
//         return x7;
//       case 7:
//         return x8;
//       case 8:
//         return x9;
//     }
//     std::unreachable();
//   }
// };
// }  // namespace linalg

namespace hcpwa {

template <int N>
  requires(sizeof(linalg::Vector<hcpwa::Float, N>) == sizeof(hcpwa::Float) * N
           && alignof(linalg::Vector<hcpwa::Float, N>) == alignof(hcpwa::Float))
// using Vec = linalg::vec<hcpwa::Float, N>;
using Vec = linalg::Vector<hcpwa::Float, N>;

constexpr Float modulo(const Vec<2>& v) {  // NOLINT
  return v[0] * v[0] + v[1] * v[1];
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

struct Params {
  Float b51, F, v, w, N, b57, b84, b86;
  Float f2min;
  Float f3min;
  Float f5min;
  Float f8min;
  Float f2max;
  Float f3max;
  Float f5max;
  Float f8max;
};

template <typename T>
concept Vector = requires(T v) { []<int N>(Vec<N>) {}(v); };

template <int N>
consteval int VectorSize(Vec<N>) {
  return N;
}

template <Vector T>
consteval int VectorSize() {
  return VectorSize(T{});
}

constexpr struct {
  template <Vector T>
  constexpr operator T() const {  // NOLINT
    T result;
    constexpr int kN = VectorSize<T>();
    for (int i = 0; i < kN; i++) {
      result[i] = 0;
    }
    return result;
  }
} kZeroVec;

}  // namespace hcpwa
namespace std {
template <int N>
struct formatter<hcpwa::Vec<N>, char> {
  // NOLINTNEXTLINE
  constexpr auto parse(std::format_parse_context& ctx) {
    return ctx.begin();
  }

  // NOLINTNEXTLINE
  auto format(const hcpwa::Vec<N>& id, std::format_context& ctx) const {
    std::format_to(ctx.out(), "<");
    for (int i = 0; i < N; i++) {
      if (i + 1 < N) {
        std::format_to(ctx.out(), "{}, ", id[i]);

      } else {
        std::format_to(ctx.out(), "{}", id[i]);
      }
    }
    return std::format_to(ctx.out(), ">");
  }
};

template <int N>
struct formatter<std::vector<hcpwa::Vec<N>>, char> {
  // NOLINTNEXTLINE
  constexpr auto parse(std::format_parse_context& ctx) {
    return ctx.begin();
  }

  // NOLINTNEXTLINE
  auto format(const std::vector<hcpwa::Vec<N>>& ids,
              std::format_context& ctx) const {
    std::format_to(ctx.out(), "[");
    for (const auto& id : ids) {
      std::format_to(ctx.out(), "\n<");
      for (int i = 0; i < N; i++) {
        if (i + 1 < N) {
          std::format_to(ctx.out(), "{}, ", id[i]);

        } else {
          std::format_to(ctx.out(), "{}", id[i]);
        }
      }
      std::format_to(ctx.out(), ">");
    }
    return std::format_to(ctx.out(), "\n]");
  }
};
}  // namespace std