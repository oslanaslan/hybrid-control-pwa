#pragma once

#include <types.hpp>

namespace hcpwa {

template <int N, int M>
constexpr Vec<N> Promote(const Vec<M>& v, const std::array<int, M>& pos) {
  Vec<N> result(0);
  for (int i = 0; i < M; i++) {
    result[pos[i]] = v[i];
  }
  return result;
}

}  // namespace hcpwa