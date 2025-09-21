#pragma once

#include <bitset>
#include <types.hpp>

namespace hcpwa {

template <int N, int M>
constexpr Vec<N> Promote(const Vec<M>& v, const std::array<int, M>& pos) {
  Vec<N> result = kZeroVec;
  for (int i = 0; i < M; i++) {
    result[pos[i]] = v[i];
  }
  return result;
}

template <int N, typename T>
  requires(N >= VectorSize<T>() - 1)
constexpr Line<N> Expand(const T& line) {
  constexpr int kM = VectorSize<T>() - 1;
  Line<N> result = kZeroVec;
  for (int i = 0; i < kM; i++) {
    result[i] = line[i];
  }
  result[N] = line[kM];
  return result;
}

template <int N>
hcpwa::LineSet<N> AABBBounds(const AABB<N>& aabb) {
  hcpwa::LineSet<N> result;
  result.reserve(2 * N);
  for (int i = 0; i < N; i++) {
    hcpwa::Line<N> near(0);
    hcpwa::Line<N> far(0);
    far[i] = 1;
    far[N] = -aabb.second[i];
    near[i] = -1;
    near[N] = aabb.first[i];
    result.push_back(near);
    result.push_back(far);
  }
  return result;
}


template <int N>
constexpr auto ResolutionsToMasks(
    const std::vector<std::pair<LineSet<N>, Line<N>>>& resolutions) {
  LineSet<N> lineset;
  std::vector<std::bitset<64>> masks;

  for (std::size_t i = 0; i < resolutions.size(); i++) {
    std::size_t start = lineset.size();
    auto& ls = resolutions[i].first;
    lineset.append_range(ls);
    std::bitset<64> mask;
    for (int n = 0; n < ls.size(); ++n) {
      mask.set(n + start);
    }
    masks.push_back(mask);
  }
  return std::pair{std::move(lineset), std::move(masks)};
}

template <int N, int M>
constexpr auto DimensionCast(const Line<M>& line, const std::array<int, N>& pos) {
  Line<N> result = kZeroVec;
  for (int i = 0; i < N; i++) {
    result[i] = line[pos[i]];
  }
  result[N] = line[M];
  return result;
}
template <int N, int M>
constexpr auto DimensionCast(const LineSet<M>& lines, const std::array<int, N>& pos) {
  LineSet<N> result;
  for (const auto& line : lines) {
    result.push_back(DimensionCast<N, M>(line, pos));
  }
  return result;
}
}  // namespace hcpwa