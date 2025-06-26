#pragma once

#include <span>
#include <vector>
#include <utility>
#include <ranges>

namespace cddwrap {

template <typename T>
struct matrix {  // NOLINT
  constexpr matrix(std::size_t cols, std::initializer_list<T> elems)
      : values_(elems),
        cols_(cols) {
    std::size_t rows = values_.size() / cols + (values_.size() % cols);
    while (values_.size() < cols * rows) {
      values_.emplace_back();
    }
  }
  constexpr explicit matrix(std::size_t cols, std::ranges::range auto elems)
      : values_(elems),
        cols_(cols) {
    std::size_t rows = values_.size() / cols + (values_.size() % cols);
    while (values_.size() < rows * cols) {
      values_.emplace_back();
    }
  }
  constexpr matrix(std::size_t cols, std::size_t rows) {
    cols_ = cols;
    values_.resize(rows * cols);
  }

  constexpr matrix(const matrix& other)
      : values_(other.values_),
        cols_(other.cols_) {
  }

  constexpr matrix(matrix&& other)
      : values_(std::move(other.values_)),
        cols_(other.cols_) {
  }

  matrix& operator=(const matrix& other) {
    if (&other == this) {
      return *this;
    }
    values_ = other.values_;
    cols_ = other.cols_;
    return *this;
  }

  matrix& operator=(matrix&& other) {
    if (&other == this) {
      return *this;
    }
    swap(other);
    return *this;
  }

  void swap(matrix& other) {
    values_.swap(other.values_);
    cols_ = std::exchange(other.cols_, cols_);
  }

  constexpr T& operator[](std::size_t row, std::size_t col) {
    return values_[row * cols_ + col];
  }

  constexpr const T& operator[](std::size_t row, std::size_t col) const {
    return values_[row * cols_ + col];
  }

  constexpr std::span<T> operator[](std::size_t row) {
    return std::span<T>(values_.begin() + (row * cols_), cols_);
  }

  constexpr std::span<const T> operator[](std::size_t row) const {
    return std::span<const T>(values_.begin() + (row * cols_), cols_);
  }

  std::size_t rows() const {
    return values_.size() / cols_;
  }

  std::size_t cols() const {
    return cols_;
  }

  template <std::ranges::range R>
  void push_back(R&& row) {
    values_.append_range(row | std::ranges::views::take(cols_));
    while (values_.size() % cols_ != 0) {
      values_.emplace_back();
    }
  }

  void emplace_back() {
    values_.resize(values_.size() + cols_);
  }

  std::span<T> back() {
    return std::span<T>(values_.end() - cols_, cols_);
  }

  std::span<const T> back() const {
    return std::span<T>(values_.end() - cols_, cols_);
  }

 private:
  std::vector<T> values_;
  std::size_t cols_;
};
}  // namespace cddwrap