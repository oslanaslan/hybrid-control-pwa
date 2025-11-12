#pragma once

#include <cstddef>
#include <vector>

namespace hcpwa {

template <typename T, auto P>
using FixedPrecisionComparator = decltype([](T a, T b) {
  if constexpr (requires { modulo(a - b); }) {
    return modulo(a - b) < P;
  } else {
    return abs(a - b) < P;
  }
});
template <typename T>
using DefaultComparator = decltype([](T a, T b) {
  return a == b;
});

template <typename T>
struct UniquePool {
 public:
  using Comp = std::add_pointer_t<bool(T, T)>;
  constexpr UniquePool()
      : comparator_(DefaultComparator<T>{}) {
  }

  constexpr explicit UniquePool(Comp comparator)
      : comparator_(comparator) {
  }

  constexpr T Unique(T value) {
    for (T elem : elements_) {
      if (comparator_(elem, value)) {
        return elem;
      }
    }
    elements_.push_back(value);
    return value;
  }

  constexpr std::size_t Index(T value) {
    for (std::size_t i = 0; i < elements_.size(); i++) {
      if (comparator_(elements_[i], value)) {
        return i;
      }
    }
    return static_cast<std::size_t>(-1);
  }

 private:
  std::vector<T> elements_;
  Comp comparator_;
};
}  // namespace hcpwa