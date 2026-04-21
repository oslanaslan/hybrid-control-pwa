#pragma once

#include <chrono>
#include <cstddef>
#include <utility>
#include <functional>
#include <array>
#include <stdexcept>
#include <vector>

template <typename T>
struct defer {
  // NOLINTNEXTLINE
  constexpr defer(T&& callback)
      : callback_(std::forward<T>(callback)) {
  }
  constexpr ~defer() {
    std::invoke(std::forward<T>(callback_));
  }

 private:
  T callback_;
};

struct Timer {
  decltype(std::chrono::system_clock::now()) start_time;

  void Start() {
    start_time = std::chrono::system_clock::now();
  }

  size_t Elapsed() {
    return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time).count();
  }

  size_t Round() {
    auto res = Elapsed();
    Start();
    return res;
  }
};

// vec is [V_0..V_{space_dim-1}, v]; returns {min(V), max(V), v}.
inline std::array<double, 3> MinMaxVAndAffineTerm(
    const std::vector<double>& vec, int space_dim = 8) {
  const std::size_t need = static_cast<std::size_t>(space_dim) + 1;
  if (vec.size() < need) {
    throw std::invalid_argument(
        "MinMaxVAndAffineTerm: vector must have at least space_dim+1 elements");
  }
  double min_val = vec[0];
  double max_val = vec[0];
  for (int i = 1; i < space_dim; ++i) {
    const std::size_t j = static_cast<std::size_t>(i);
    if (vec[j] < min_val) {
      min_val = vec[j];
    }
    if (vec[j] > max_val) {
      max_val = vec[j];
    }
  }
  return {min_val, max_val, vec[static_cast<std::size_t>(space_dim)]};
}

#define temp_variable(S) v_##S##__
#define temp_variable_line(L) temp_variable(L)
#define _ temp_variable_line(__LINE__)