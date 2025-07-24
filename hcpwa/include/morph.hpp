#pragma once

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

}  // namespace hcpwa