#include <gtest/gtest.h>
#include <cstdio>
#include <format>
#include <hcpwa.hpp>
#include <symbolic.hpp>
#include "algo.hpp"
#include "morph.hpp"
#include "test_utils.hpp"
#include "types.hpp"
#include "cddwrap/cdd.hpp"
#include "utility.hpp"
#include <uniqie_pool.hpp>

TEST(complex, symbolic2d) {
  // NOLINTNEXTLINE
  using namespace hcpwa::symbols;

  cddwrap::global_init();
  defer _ = &cddwrap::global_free;
  hcpwa::AABB<2> aabb = {{0, 0}, {2, 2}};

  {
    constexpr auto x = X<0>{};  // NOLINT
    constexpr auto y = X<1>{};  // NOLINT
    auto v = SymMin(3 * x + y, x + y + 1);

    auto resolutions = MinResolutions(v);
    constexpr int kDim
        = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);

    ASSERT_EQ(masks.size(), 2);
    ASSERT_EQ(masks[0], 0b01);
    ASSERT_EQ(masks[1], 0b10);

    auto polygons = hcpwa::SplitAABBWithLines(aabb, lines);

    using Comp = hcpwa::FixedPrecisionComparator<hcpwa::Vec<2>, 0.01f>;
    hcpwa::UniquePool<hcpwa::Vec<2>> pool(Comp{});
    pool.Unique(aabb.first);
    pool.Unique(aabb.second);
    pool.Unique({aabb.first[0], aabb.second[1]});
    pool.Unique({aabb.second[0], aabb.first[1]});

    hcpwa::NormalizeVertices(polygons, pool);

    std::vector<hcpwa::Line<2>> polygon_line;

    for (std::size_t i = 0; i < polygons.size(); i++) {
      for (std::size_t mask_index = 0; mask_index < masks.size();
           mask_index++) {
        if (((~polygons[i].is_gt) & masks[mask_index]) == masks[mask_index]) {
          polygon_line.push_back(resolutions[mask_index].second);
          break;
        }
      }
    }

    for (int i = 0; i < polygons.size(); i++) {
      GTEST_COUT << std::format("poly {}: {}\n", i, polygons[i].polygon);
    }

    GTEST_COUT << std::format("{}\n", polygon_line[0]);
    GTEST_COUT << std::format("{}\n", polygon_line[1]);
  }
}

TEST(complex, symbolic3d) {
  // NOLINTNEXTLINE
  using namespace hcpwa::symbols;

  cddwrap::global_init();
  defer _ = &cddwrap::global_free;
  hcpwa::AABB<3> aabb = {{0, 0, 0}, {1, 1, 1}};

  {
    constexpr auto x = X<0>{};  // NOLINT
    constexpr auto y = X<1>{};  // NOLINT
    constexpr auto z = X<2>{};  // NOLINT
    auto vxy = SymMin(3 * x + y, x + y + 1);

    auto resolutions = MinResolutions(vxy);
    constexpr int kDim
        = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    auto polygons
        = hcpwa::SplitAABBWithLines({aabb.first[0, 1], aabb.second[0, 1]}, lines);

    using Comp = hcpwa::FixedPrecisionComparator<hcpwa::Vec<2>, 0.01f>;
    hcpwa::UniquePool<hcpwa::Vec<2>> pool(Comp{});
    pool.Unique(aabb.first[0, 1]);
    pool.Unique(aabb.second[0, 1]);
    pool.Unique({aabb.first[0], aabb.second[1]});
    pool.Unique({aabb.second[0], aabb.first[1]});

    hcpwa::NormalizeVertices(polygons, pool);

    std::vector<hcpwa::Line<2>> polygon_line;

    for (std::size_t i = 0; i < polygons.size(); i++) {
      for (std::size_t mask_index = 0; mask_index < masks.size();
           mask_index++) {
        if (((~polygons[i].is_gt) & masks[mask_index]) == masks[mask_index]) {
          polygon_line.push_back(resolutions[mask_index].second);
          break;
        }
      }
    }

    for (int i = 0; i < polygons.size(); i++) {
      GTEST_COUT << std::format("poly {}: {}\n", i, polygons[i].polygon);
    }

    GTEST_COUT << std::format("{}\n", polygon_line[0]);
    GTEST_COUT << std::format("{}\n", polygon_line[1]);
  }
}

TEST(complex, intersect) {
  // NOLINTNEXTLINE
  using namespace hcpwa::symbols;

  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  constexpr hcpwa::Float N = 100;
  constexpr hcpwa::Float F = 20;
  constexpr hcpwa::Float v = 10;
  constexpr hcpwa::Float w = 10;
  constexpr hcpwa::Float b51 = 10;
  constexpr hcpwa::Float b57 = 10;
  constexpr hcpwa::Float b84 = 10;
  constexpr hcpwa::Float b86 = 10;
  constexpr hcpwa::Float f2min = 10;
  constexpr hcpwa::Float f3min = 10;
  constexpr hcpwa::Float f5min = 10;
  constexpr hcpwa::Float f8min = 10;
  constexpr hcpwa::Float f2max = 10;
  constexpr hcpwa::Float f3max = 10;
  constexpr hcpwa::Float f5max = 10;
  constexpr hcpwa::Float f8max = 10;

  hcpwa::AABB<8> aabb = {{0, 0, 0, 0, 0, 0, 0, 0}, {N, N, N, N, N, N, N, N}};

  constexpr auto n1 = X<0>{};
  constexpr auto n2 = X<1>{};
  constexpr auto n3 = X<2>{};
  constexpr auto n4 = X<3>{};
  constexpr auto n5 = X<4>{};
  constexpr auto n6 = X<5>{};
  constexpr auto n7 = X<6>{};
  constexpr auto n8 = X<7>{};

  // f51(t) = min{β51F,β51vn5(t),w(N−n1(t))}
  auto f51 = SymMin(b51 * F, (v * b51) * n5, w * (N - n1));
  // f57(t) = min{β57F,β57vn5(t),w(N−n7(t))}
  auto f57 = SymMin(b57 * F, v * b57 * n5, w * (N - n7));
  // f84(t) = min{β84F,β84vn8(t),w(N−n4(t))}
  auto f84 = SymMin(b84 * F, v * b84 * n8, w * (N - n4));
  // f86(t) = min{β86F,β86vn8(t),w(N−n6(t))}
  auto f86 = SymMin(b86 * F, v * b86 * n8, w * (N - n6));

}