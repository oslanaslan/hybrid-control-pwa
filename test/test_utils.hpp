#pragma once

#include <ostream>
#include <vector>
#include <types.hpp>

#define GTEST_COUT std::cerr << "[          ] [ INFO ]"

namespace tu {
template <typename T>
struct PrettyPrint {
  const T& ref;
  constexpr explicit PrettyPrint(const T& ref)
      : ref(ref) {
  }
  PrettyPrint(const PrettyPrint&) = delete;
  PrettyPrint(PrettyPrint&&) = delete;
};
template <typename T>
PrettyPrint(const T&) -> PrettyPrint<T>;

template <typename T>
  requires requires { std::declval<std::ostream>() << std::declval<T>(); }
constexpr std::ostream& operator<<(std::ostream& os, const PrettyPrint<T>& pp) {
  return os << pp.ref;
}

template <typename T>
using VecValueType = decltype([]<typename F, int N>(linalg::Vector<F, N>) -> F {
  throw 1;
}(std::declval<T>()));

template <int N>
constexpr std::ostream& operator<<(std::ostream& os,
                                   const PrettyPrint<hcpwa::Vec<N>>& pp) {
  // using ValueType = VecValueType<hcpwa::Vec<N>>;
  using ValueType = hcpwa::Vec<N>::ValueType;
  constexpr std::size_t kCount = sizeof(hcpwa::Vec<N>) / sizeof(ValueType);
  os << "Vec<" << N << ">{";
  for (std::size_t i = 0; i < kCount - 1; i++) {
    if constexpr (requires { os << std::declval<ValueType>(); }) {
      os << PrettyPrint((reinterpret_cast<const ValueType*>(&pp.ref)[i]))
         << ", ";
    }
  }
  os << PrettyPrint(reinterpret_cast<const ValueType*>(&pp.ref)[kCount - 1])
     << "}";
  return os;
}

template <typename T>
constexpr std::ostream& operator<<(std::ostream& os,
                                   const PrettyPrint<std::vector<T>>& pp) {
  os << "[";
  for (auto it = pp.ref.begin(); it != pp.ref.end(); ++it) {
    os << PrettyPrint(*it);
    if (std::next(it) != pp.ref.end()) {
      os << ", ";
    }
  }
  os << "]";
  return os;
}

constexpr std::ostream& operator<<(
    std::ostream& os, const PrettyPrint<hcpwa::PolygonResolution>& pp) {
  return os << "PolygonResolution{\n.is_gt = " << pp.ref.is_gt.to_string()
            << "\n.polygon = " << PrettyPrint(pp.ref.polygon) << "\n}";
}
}  // namespace tu