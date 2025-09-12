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
    pool.Unique({aabb.first.x, aabb.second.y});
    pool.Unique({aabb.second.x, aabb.first.y});

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
        = hcpwa::SplitAABBWithLines({aabb.first.xy(), aabb.second.xy()}, lines);

    using Comp = hcpwa::FixedPrecisionComparator<hcpwa::Vec<2>, 0.01f>;
    hcpwa::UniquePool<hcpwa::Vec<2>> pool(Comp{});
    pool.Unique(aabb.first.xy());
    pool.Unique(aabb.second.xy());
    pool.Unique({aabb.first.x, aabb.second.y});
    pool.Unique({aabb.second.x, aabb.first.y});

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
  
}