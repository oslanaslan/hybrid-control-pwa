#pragma once

#include <chrono>
#include <cstddef>
#include <utility>
#include <functional>

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

#define temp_variable(S) v_##S##__
#define temp_variable_line(L) temp_variable(L)
#define _ temp_variable_line(__LINE__)