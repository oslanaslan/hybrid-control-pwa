#include <gtest/gtest.h>
#include <algorithm>
#include <cddwrap/lineareq.hpp>
#include <cddwrap/cdd.hpp>
#include "utility.hpp"

static bool CheckPoints(const cddwrap::matrix<double>& m,
                        std::vector<std::vector<double>> points) {
  for (std::size_t i = 0; i < m.rows(); i++) {
    for (std::size_t j = 0; j < points.size(); j++) {
      if (std::ranges::equal(m[i], points[j])) {
        goto next;
      }
    }
    return false;
  next:
  }
  return true;
}

TEST(cdd, hull2d) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;
  // Ax<b <=> -Ax+b >= 0
  // Format: A|b
  {
    cddwrap::matrix<double> values(3, {1, 0, 1, 0, 1, 1, -1, 0, 0, 0, -1, 0});
    ASSERT_EQ(values.rows(), 4);
    ASSERT_EQ(values.cols(), 3);

    auto result = cddwrap::GetHullPoints(values);
    ASSERT_GT(result.cols(), 0) << result.rows() << " x " << result.cols();
    for (std::size_t i = 0; i < result.rows(); i++) {
      for (std::size_t j = 0; j < result.cols(); j++) {
        std::cerr << result[i, j] << ' ';
      }
      std::cerr << '\n';
    }
    ASSERT_EQ(result.cols(), 2);
    ASSERT_EQ(result.rows(), 4);
    ASSERT_TRUE(CheckPoints(result, {{0, 0}, {0, 1}, {1, 0}, {1, 1}}));
  }
  {
    cddwrap::matrix<double> values(3, {1, 1, 2, 0, 1, 1, -1, 0, 0, 0, -1, 0});
    ASSERT_EQ(values.rows(), 4);
    ASSERT_EQ(values.cols(), 3);

    auto result = GetHullPoints(values);
    ASSERT_GT(result.cols(), 0) << result.rows() << " x " << result.cols();
    for (std::size_t i = 0; i < result.rows(); i++) {
      for (std::size_t j = 0; j < result.cols(); j++) {
        std::cerr << result[i, j] << ' ';
      }
      std::cerr << '\n';
    }
    ASSERT_EQ(result.cols(), 2);
    ASSERT_EQ(result.rows(), 4);
    ASSERT_TRUE(CheckPoints(result, {{0, 0}, {0, 1}, {2, 0}, {1, 1}}));
  }
}