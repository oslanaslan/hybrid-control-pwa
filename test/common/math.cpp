#include <gtest/gtest.h>
#include <algorithm>
#include <types.hpp>
#include "algo.hpp"
#include "cddwrap/cdd.hpp"
#include "morph.hpp"
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

TEST(math, line) {
  hcpwa::Line<2> a{1, 0, 2};
  hcpwa::Line<8> b = hcpwa::Expand<8>(a);
  ASSERT_EQ(b[0], a[0]);
  ASSERT_EQ(b[8], a[2]);
}

TEST(math, floating) {
  auto test = []<typename T>(T) {
    T a = 1;
    ASSERT_EQ(a, 1);
    T b = 2;
    T c = a + b;
    ASSERT_GE(c, a);
    ASSERT_GE(c, b);
    ASSERT_LT(a, b);
    auto x = c + a;
    auto y = c - a;
    static_assert(std::same_as<decltype(x), T>);
    static_assert(std::same_as<decltype(y), T>);

    float m = 1.5;
    a += m;
    b -= m;
    c = m + a;
    c = m - b;
    c = m + 1;
  };
  test(float{0});
  test(double{0});
  test(hcpwa::CustomFloat<float>(0));
  test(hcpwa::CustomFloat<double>(0));
}

TEST(math, split) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;
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

TEST(math, normalize) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;
  hcpwa::AABB<2> aabb = {{0, 0}, {1, 1}};
  hcpwa::LineSet<2> lines = {
      {1, 1, -1},
      {1, -1, 0},
      {0, 1, -0.2},
      {1, 0, -0.999},
  };
  auto polygons = hcpwa::SplitAABBWithLines(aabb, lines);

  ASSERT_EQ(polygons.size(), 11) << tu::PrettyPrint(polygons);

  using Comp = hcpwa::FixedPrecisionComparator<hcpwa::Vec<2>, 0.01f>;
  hcpwa::UniquePool<hcpwa::Vec<2>> pool(Comp{});
  pool.Unique(aabb.first);
  pool.Unique(aabb.second);
  pool.Unique({aabb.first[0], aabb.second[1]});
  pool.Unique({aabb.second[0], aabb.first[1]});

  hcpwa::NormalizeVertices(polygons, pool);

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

TEST(math, triangulate) {
  std::vector<hcpwa::PolygonResolution> data
      = {hcpwa::PolygonResolution{.polygon = {
                                      {0, 0},
                                      {0, 1},
                                      {1, 1},
                                      {1, 0},
                                  },
                                },
                              hcpwa::PolygonResolution{.polygon = {
                                      {0, 0},
                                      {0, 2},
                                      {2, 2},
                                      {2, 0},
                                  },
                                }};

  auto tri = hcpwa::Triangulate(data);

  ASSERT_EQ(tri.size(), 4);
  ASSERT_EQ(tri[0].second, 0);
  ASSERT_EQ(tri[1].second, 0);
  ASSERT_EQ(tri[2].second, 1);
  ASSERT_EQ(tri[3].second, 1);
}

TEST(math, prism1) {
  hcpwa::Triangle triangle = {{0, 0}, {0, 1}, {1, 0}};
  auto prism = hcpwa::CalcPrism(triangle, {0, 1});

  const hcpwa::Vec<8> inside{0.4, 0.4, 0};
  const hcpwa::Vec<8> outside{-0.5, 0.5, 0};

  ASSERT_EQ(prism.size(), 3);

  ASSERT_LT(hcpwa::Apply(prism[0], inside), 0);
  ASSERT_LT(hcpwa::Apply(prism[1], inside), 0);
  ASSERT_LT(hcpwa::Apply(prism[2], inside), 0);

  ASSERT_TRUE(hcpwa::Apply(prism[0], outside) > 0
              || hcpwa::Apply(prism[1], outside) > 0
              || hcpwa::Apply(prism[2], outside) > 0);
}

TEST(math, prism2) {
  hcpwa::Triangle triangle = {{0, 0}, {0, 1}, {1, 0}};
  auto prism = hcpwa::CalcPrism(triangle, {0, 2});

  const hcpwa::Vec<8> inside{0.4, 1, 0.4, 0};
  const hcpwa::Vec<8> outside{-0.5, 0, 0.5, 0};

  ASSERT_EQ(prism.size(), 3);

  ASSERT_LT(hcpwa::Apply(prism[0], inside), 0);
  ASSERT_LT(hcpwa::Apply(prism[1], inside), 0);
  ASSERT_LT(hcpwa::Apply(prism[2], inside), 0);

  ASSERT_TRUE(hcpwa::Apply(prism[0], outside) > 0
              || hcpwa::Apply(prism[1], outside) > 0
              || hcpwa::Apply(prism[2], outside) > 0);
}

TEST(math, hull8d) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  hcpwa::Triangle triangle = {{0, 0}, {0, 1}, {1, 0}};
  auto prism = hcpwa::CalcPrism(triangle, {0, 2});

  ASSERT_EQ(prism.size(), 3);

  const hcpwa::AABB<8> aabb{{0, 0, 0, 0, 0, 0, 0, 0}, {1, 1, 1, 1, 1, 1, 1, 1}};
  auto aabb_lines = hcpwa::AABBBounds(aabb);

  for (auto& v : prism) {
    aabb_lines.push_back(v);
  }

  auto points = hcpwa::LinesToPoints(aabb_lines);

  ASSERT_EQ(points.size(), 192);
}