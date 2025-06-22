#include <gtest/gtest.h>
#include <algorithm>
#include <types.hpp>
#include "algo.hpp"
#include "cdd.hpp"
#include "utility.hpp"
#include "test_utils.hpp"

TEST(math, vec) {
  hcpwa::Vec<2> a{1, 2};
  hcpwa::Vec<2> b{1, 2};
  hcpwa::Vec<2> c{2, 4};
  ASSERT_EQ(a - b, hcpwa::Vec<2>{});
  ASSERT_EQ(a, b);
  ASSERT_EQ(a + b, c);
}

TEST(math, split) {
  GlobalInit();
  defer _ = &GlobalFree;
  hcpwa::AABB<2> aabb = {{0, 0}, {1, 1}};
  hcpwa::LineSet<2> lines = {
      {1, 1, -1},
      {1, -1, 0},
      {0, 1, -0.2},
  };
  auto polygons = hcpwa::SplitAABBWithLines(aabb, lines);

  ASSERT_EQ(polygons.size(), 7) << tu::PrettyPrint(polygons);

  const std::size_t triangle_count
      = std::ranges::count_if(polygons, [](const hcpwa::PolygonResolution& p) {
          return p.polygon.size() == 3;
        });
  const std::size_t quad_count
      = std::ranges::count_if(polygons, [](const hcpwa::PolygonResolution& p) {
          return p.polygon.size() == 4;
        });

  ASSERT_EQ(triangle_count, 4);
  ASSERT_EQ(quad_count, 3);
}