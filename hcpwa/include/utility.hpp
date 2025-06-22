#pragma once

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

#define temp_variable(S) v_##S##__
#define temp_variable_line(L) temp_variable(L)
#define _ temp_variable_line(__LINE__)