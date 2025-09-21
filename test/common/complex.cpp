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
  hcpwa::AABB<2> aabb2d = {{0, 0}, {N, N}};

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
  // f1out(t) = min{vn1(τ),F}
  auto f1out = SymMin(v * n1, F);
  // f4out(t) = min{vn4(τ),F}
  auto f4out = SymMin(v * n4, F);
  // f6out(t) = min{vn6(τ),F}
  auto f6out = SymMin(v * n6, F);
  // f7out(t) = min{vn7(τ),F}
  auto f7out = SymMin(v * n7, F);
  // f_i_in(t) = min{f_i_min,w(N−n_i(τ))} for i = 2, 3, 5, 8
  auto f2in = SymMin(f2min, w * (N - n2));
  auto f3in = SymMin(f3min, w * (N - n3));
  auto f5in = SymMin(f5min, w * (N - n5));
  auto f8in = SymMin(f8min, w * (N - n8));

  // Plane (5, 1)
  hcpwa::LineSet<8> lines51;
  {
    auto resolutions = MinResolutions(f51);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines51.append_range(lines);
  }
  {
    auto resolutions = MinResolutions(f1out);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines51.append_range(lines);
  }
  {
    auto resolutions = MinResolutions(f5in);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines51.append_range(lines);
  }
  // Plane (5, 7)
  hcpwa::LineSet<8> lines57;
  {
    auto resolutions = MinResolutions(f57);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines57.append_range(lines);
  }
  {
    auto resolutions = MinResolutions(f5in);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines57.append_range(lines);
  }
  {
    auto resolutions = MinResolutions(f7out);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines57.append_range(lines);
  }
  // Plane (8, 4)
  hcpwa::LineSet<8> lines84;
  {
    auto resolutions = MinResolutions(f84);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines84.append_range(lines);
  }
  {
    auto resolutions = MinResolutions(f4out);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines84.append_range(lines);
  }
  {
    auto resolutions = MinResolutions(f8in);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines84.append_range(lines);
  }
  // Plane (8, 6)
  hcpwa::LineSet<8> lines86;
  {
    auto resolutions = MinResolutions(f86);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines86.append_range(lines);
  }
  {
    auto resolutions = MinResolutions(f6out);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines86.append_range(lines);
  }
  {
    auto resolutions = MinResolutions(f8in);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines86.append_range(lines);
  }

  auto polygon51 = hcpwa::SplitAABBWithLines(aabb2d, hcpwa::DimensionCast<2, 8>(lines51, {1 - 1, 5 - 1}));
  auto polygon57 = hcpwa::SplitAABBWithLines(aabb2d, hcpwa::DimensionCast<2, 8>(lines57, {5 - 1, 7 - 1}));
  auto polygon84 = hcpwa::SplitAABBWithLines(aabb2d, hcpwa::DimensionCast<2, 8>(lines84, {4 - 1, 8 - 1}));
  auto polygon86 = hcpwa::SplitAABBWithLines(aabb2d, hcpwa::DimensionCast<2, 8>(lines86, {6 - 1, 8 - 1}));


  GTEST_COUT << "Original polygons:\n";
  for (int i = 0; i < polygon51.size(); i++) {
    GTEST_COUT << std::format("poly {}: {}\n", i, polygon51[i].polygon);
  }
  for (int i = 0; i < polygon57.size(); i++) {
    GTEST_COUT << std::format("poly {}: {}\n", i, polygon57[i].polygon);
  }
  for (int i = 0; i < polygon84.size(); i++) {
    GTEST_COUT << std::format("poly {}: {}\n", i, polygon84[i].polygon);
  }
  for (int i = 0; i < polygon86.size(); i++) {
    GTEST_COUT << std::format("poly {}: {}\n", i, polygon86[i].polygon);
  }

  // Triangulated polygons with unique vertices
  auto triangles51 = hcpwa::TrianglesWithUniqueVertices(aabb2d, polygon51);
  auto triangles57 = hcpwa::TrianglesWithUniqueVertices(aabb2d, polygon57);
  auto triangles84 = hcpwa::TrianglesWithUniqueVertices(aabb2d, polygon84);
  auto triangles86 = hcpwa::TrianglesWithUniqueVertices(aabb2d, polygon86);

  GTEST_COUT << "Unique vertices polygons:\n";
  for (int i = 0; i < polygon51.size(); i++) {
    GTEST_COUT << std::format("poly {}: {}\n", i, polygon51[i].polygon);
  }
  for (int i = 0; i < polygon57.size(); i++) {
    GTEST_COUT << std::format("poly {}: {}\n", i, polygon57[i].polygon);
  }
  for (int i = 0; i < polygon84.size(); i++) {
    GTEST_COUT << std::format("poly {}: {}\n", i, polygon84[i].polygon);
  }
  for (int i = 0; i < polygon86.size(); i++) {
    GTEST_COUT << std::format("poly {}: {}\n", i, polygon86[i].polygon);
  }

  GTEST_COUT << "Triangulated polygons with unique vertices:\n";
  for (int i = 0; i < triangles51.size(); i++) {
    GTEST_COUT << std::format("poly {}: {}\n", i, triangles51[i].first.a);
    GTEST_COUT << std::format("poly {}: {}\n", i, triangles51[i].first.b);
    GTEST_COUT << std::format("poly {}: {}\n", i, triangles51[i].first.c);
  }
  for (int i = 0; i < triangles57.size(); i++) {
    GTEST_COUT << std::format("poly {}: {}\n", i, triangles57[i].first.a);
    GTEST_COUT << std::format("poly {}: {}\n", i, triangles57[i].first.b);
    GTEST_COUT << std::format("poly {}: {}\n", i, triangles57[i].first.c);
  }
  for (int i = 0; i < triangles84.size(); i++) {
    GTEST_COUT << std::format("poly {}: {}\n", i, triangles84[i].first.a);
    GTEST_COUT << std::format("poly {}: {}\n", i, triangles84[i].first.b);
    GTEST_COUT << std::format("poly {}: {}\n", i, triangles84[i].first.c);
  }
  for (int i = 0; i < triangles86.size(); i++) {
    GTEST_COUT << std::format("poly {}: {}\n", i, triangles86[i].first.a);
    GTEST_COUT << std::format("poly {}: {}\n", i, triangles86[i].first.b);
    GTEST_COUT << std::format("poly {}: {}\n", i, triangles86[i].first.c);
  }

  // Calculate prisms
  std::vector<hcpwa::LineSet<8>> prisms51;
  std::vector<hcpwa::LineSet<8>> prisms57;
  std::vector<hcpwa::LineSet<8>> prisms84;
  std::vector<hcpwa::LineSet<8>> prisms86;
  for (auto& triangle : triangles51) {
    auto prism = hcpwa::CalcPrism(triangle.first, {1 - 1, 5 - 1});
    prisms51.push_back(prism);
  }
  for (auto& triangle : triangles57) {
    auto prism = hcpwa::CalcPrism(triangle.first, {5 - 1, 7 - 1});
    prisms57.push_back(prism);
  }
  for (auto& triangle : triangles84) {
    auto prism = hcpwa::CalcPrism(triangle.first, {4 - 1, 8 - 1});
    prisms84.push_back(prism);
  }
  for (auto& triangle : triangles86) {
    auto prism = hcpwa::CalcPrism(triangle.first, {6 - 1, 8 - 1});
    prisms86.push_back(prism);
  }

  GTEST_COUT << "Prisms count: " << prisms51.size() << " " << prisms57.size() << " " << prisms84.size() << " " << prisms86.size() << "\n";

  // Intersect all prisms
  std::vector<hcpwa::LineSet<8>> intersections;
  std::vector<std::vector<hcpwa::Vec<8>>> intersection_points;
  for (auto& prism51 : prisms51) {
    for (auto& prism57 : prisms57) {
      for (auto& prism84 : prisms84) {
        for (auto& prism86 : prisms86) {
          hcpwa::LineSet<8> concatenated_prisms;
          concatenated_prisms.append_range(prism51);
          concatenated_prisms.append_range(prism57);
          concatenated_prisms.append_range(prism84);
          concatenated_prisms.append_range(prism86);
          auto intersection = hcpwa::LinesToPoints(concatenated_prisms);
          if (intersection.size() > 0) {
            intersections.push_back(concatenated_prisms);
            intersection_points.push_back(intersection);
          }
        }
      }
    }
  }
  GTEST_COUT << "Intersections count: " << intersections.size() << "\n";
}