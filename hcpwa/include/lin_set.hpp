#pragma once

#include <array>
#include <vector>
#include <algorithm>
#include <stdexcept>

namespace hcpwa {

template <typename T, size_t N>
struct LinearSet {
  void Insert(T value) {
    if (tail_ >= N) {
      throw std::runtime_error("out of bounds");
    }
    if (!std::ranges::contains(values_, value)) {
      values_[tail_++] = std::move(value);
    }
  }

  void Remove(const T& value) {
    for (std::size_t i = 0; i < tail_; i++) {
      if (values_[i] == value) {
        swap(values_[i], values_[tail_ - 1]);
        tail_--;
        return;
      }
    }
  }

 private:
  std::array<T, N> values_;
  std::size_t tail_ = 0;
};

template <typename T>
struct LinearSet<T, 0> {
  void Insert(T value) {
    if (!std::ranges::contains(values_, value)) {
      values_.push_back(std::move(value));
    }
  }

  void Remove(T value) {
    for (std::size_t i = 0; i < values_.size(); i++) {
      if (values_[i] == value) {
        swap(values_[i], values_.back());
        values_.pop_back();
        return;
      }
    }
  }

 private:
  std::vector<T> values_;
};
}  // namespace hcpwa