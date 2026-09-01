#include <gtest/gtest.h>
#include <cstddef>
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

namespace {

hcpwa::Vec<8> MakeUniformVec8(double value) {
  hcpwa::Vec<8> v = hcpwa::kZeroVec;
  for (int d = 0; d < 8; ++d) {
    v[d] = value;
  }
  return v;
}

}  // namespace

TEST(user_algo, compute_triangle_areas_vertices) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  constexpr hcpwa::Float N = 100;
  constexpr hcpwa::Float F = 15;
  constexpr hcpwa::Float v = 0.2;
  constexpr hcpwa::Float w = 0.5;
  constexpr hcpwa::Float b51 = 0.5;
  constexpr hcpwa::Float b57 = 0.5;
  constexpr hcpwa::Float b84 = 0.5;
  constexpr hcpwa::Float b86 = 0.5;
  constexpr hcpwa::Float b31 = 0.5;
  constexpr hcpwa::Float b36 = 0.5;
  constexpr hcpwa::Float b24 = 0.5;
  constexpr hcpwa::Float b27 = 0.5;
  constexpr hcpwa::Float f2min = 5;
  constexpr hcpwa::Float f3min = 5;
  constexpr hcpwa::Float f5min = 5;
  constexpr hcpwa::Float f8min = 5;
  constexpr hcpwa::Float f2max = 10;
  constexpr hcpwa::Float f3max = 10;
  constexpr hcpwa::Float f5max = 10;
  constexpr hcpwa::Float f8max = 10;

  hcpwa::TriangleAreasVerticesResult result
      = hcpwa::compute_triangle_areas_vertices(
          N, F, v, w, b51, b57, b84, b86, b31, b36, b24, b27, f2min, f3min,
          f5min, f8min, f2max, f3max, f5max, f8max);

  // Print triangle sizes for phase 0
  GTEST_COUT << "result.triangles31.size(): " << result.triangles31.size()
             << '\n';
  GTEST_COUT << "result.triangles36.size(): " << result.triangles36.size()
             << '\n';
  GTEST_COUT << "result.triangles24.size(): " << result.triangles24.size()
             << '\n';
  GTEST_COUT << "result.triangles27.size(): " << result.triangles27.size()
             << '\n';
  // Print triangle sizes for phase 1
  GTEST_COUT << "result.triangles51.size(): " << result.triangles51.size()
             << '\n';
  GTEST_COUT << "result.triangles57.size(): " << result.triangles57.size()
             << '\n';
  GTEST_COUT << "result.triangles84.size(): " << result.triangles84.size()
             << '\n';
  GTEST_COUT << "result.triangles86.size(): " << result.triangles86.size()
             << '\n';
  // Print intersection points sizes for phase 0
  GTEST_COUT << "result.intersection_points_phase0.size(): "
             << result.intersection_points_phase0.size() << '\n';
  GTEST_COUT << "result.intersection_prism_indices_phase0.size(): "
             << result.intersection_prism_indices_phase0.size() << '\n';
  // Print intersection points sizes for phase 1
  GTEST_COUT << "result.intersection_points_phase1.size(): "
             << result.intersection_points_phase1.size() << '\n';
  GTEST_COUT << "result.intersection_prism_indices_phase1.size(): "
             << result.intersection_prism_indices_phase1.size() << '\n';

  // Asserts for phase 0 triangles
  ASSERT_GT(result.triangles31.size(), 0);
  ASSERT_GT(result.triangles36.size(), 0);
  ASSERT_GT(result.triangles24.size(), 0);
  ASSERT_GT(result.triangles27.size(), 0);
  // Asserts for phase 1 triangles
  ASSERT_GT(result.triangles51.size(), 0);
  ASSERT_GT(result.triangles57.size(), 0);
  ASSERT_GT(result.triangles84.size(), 0);
  ASSERT_GT(result.triangles86.size(), 0);

  // Verify that intersection points and indices have the same size for phase
  // 1
  ASSERT_EQ(result.intersection_points_phase0.size(),
            result.intersection_prism_indices_phase0.size());
  ASSERT_EQ(result.intersection_points_phase1.size(),
            result.intersection_prism_indices_phase1.size());
  // TODO Fix intersection points computing so that areas and points counts
  // matches for pahse 0 and 1
  // ASSERT_EQ(result.intersection_prism_indices_phase0.size(),
  // result.intersection_prism_indices_phase1.size());
  // ASSERT_EQ(result.intersection_points_phase0.size(),
  // result.intersection_points_phase1.size());

  // Verify that each intersection has exactly 5 prism indices for phase 1
  for (const auto& prism_indices : result.intersection_prism_indices_phase1) {
    ASSERT_EQ(prism_indices.size(), 5);
    // Verify indices are within valid ranges (phase 1 uses triangles51, 57,
    // 84, 86, 23)
    ASSERT_LT(prism_indices[0], result.triangles51.size());
    ASSERT_LT(prism_indices[1], result.triangles57.size());
    ASSERT_LT(prism_indices[2], result.triangles84.size());
    ASSERT_LT(prism_indices[3], result.triangles86.size());
    ASSERT_LT(prism_indices[4], result.triangles23.size());
  }
  // Verify that each intersection has exactly 5 prism indices for phase 0
  for (const auto& prism_indices : result.intersection_prism_indices_phase0) {
    ASSERT_EQ(prism_indices.size(), 5);
    // Verify indices are within valid ranges (phase 0 uses triangles31, 36,
    // 24, 27, 58)
    ASSERT_LT(prism_indices[0], result.triangles31.size());
    ASSERT_LT(prism_indices[1], result.triangles36.size());
    ASSERT_LT(prism_indices[2], result.triangles24.size());
    ASSERT_LT(prism_indices[3], result.triangles27.size());
    ASSERT_LT(prism_indices[4], result.triangles58.size());
  }

  // Verify that intersection points are valid (Vec<8> is a fixed-size type
  // with 8 elements)
  for (const auto& point_group : result.intersection_points_phase1) {
    ASSERT_GT(point_group.size(), 0);
  }
  for (const auto& point_group : result.intersection_points_phase0) {
    ASSERT_GT(point_group.size(), 0);
  }

  size_t areas_vertices_count_phase0 = 0;
  size_t areas_vertices_count_phase1 = 0;
  for (auto& area : result.intersection_points_phase0) {
    areas_vertices_count_phase0 += area.size();
  }
  for (auto& area : result.intersection_points_phase1) {
    areas_vertices_count_phase1 += area.size();
  }

  GTEST_COUT << "Phase 0 results:\n";
  GTEST_COUT << std::format("  triangles31: {}\n", result.triangles31.size());
  GTEST_COUT << std::format("  triangles36: {}\n", result.triangles36.size());
  GTEST_COUT << std::format("  triangles24: {}\n", result.triangles24.size());
  GTEST_COUT << std::format("  triangles27: {}\n", result.triangles27.size());
  GTEST_COUT << std::format("  intersection areas: {}\n",
                            result.intersection_points_phase0.size());
  GTEST_COUT << std::format("  intersection areas vertices: {}\n",
                            areas_vertices_count_phase0);
  GTEST_COUT << "Phase 1 results:\n";
  GTEST_COUT << std::format("  triangles51: {}\n", result.triangles51.size());
  GTEST_COUT << std::format("  triangles57: {}\n", result.triangles57.size());
  GTEST_COUT << std::format("  triangles84: {}\n", result.triangles84.size());
  GTEST_COUT << std::format("  triangles86: {}\n", result.triangles86.size());
  GTEST_COUT << std::format("  intersection areas: {}\n",
                            result.intersection_points_phase1.size());
  GTEST_COUT << std::format("  intersection areas vertices: {}\n",
                            areas_vertices_count_phase1);
}

TEST(user_algo, compute_polygon_areas_vertices) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  constexpr hcpwa::Float N = 100;
  constexpr hcpwa::Float F = 15;
  constexpr hcpwa::Float v = 0.2;
  constexpr hcpwa::Float w = 0.5;
  constexpr hcpwa::Float b51 = 0.5;
  constexpr hcpwa::Float b57 = 0.5;
  constexpr hcpwa::Float b84 = 0.5;
  constexpr hcpwa::Float b86 = 0.5;
  constexpr hcpwa::Float b31 = 0.5;
  constexpr hcpwa::Float b36 = 0.5;
  constexpr hcpwa::Float b24 = 0.5;
  constexpr hcpwa::Float b27 = 0.5;
  constexpr hcpwa::Float f2min = 5;
  constexpr hcpwa::Float f3min = 5;
  constexpr hcpwa::Float f5min = 5;
  constexpr hcpwa::Float f8min = 5;
  constexpr hcpwa::Float f2max = 10;
  constexpr hcpwa::Float f3max = 10;
  constexpr hcpwa::Float f5max = 10;
  constexpr hcpwa::Float f8max = 10;

  hcpwa::PolygonAreasVerticesResult result = hcpwa::compute_polygon_areas_vertices(
      N, F, v, w, b51, b57, b84, b86, b31, b36, b24, b27, f2min, f3min, f5min,
      f8min, f2max, f3max, f5max, f8max);

  // Print intersection points sizes for phase 0
  GTEST_COUT << "result.intersection_points_phase0.size(): "
             << result.intersection_points_phase0.size() << '\n';
  GTEST_COUT << "result.intersection_prism_indices_phase0.size(): "
             << result.intersection_prism_indices_phase0.size() << '\n';
  // Print intersection points sizes for phase 1
  GTEST_COUT << "result.intersection_points_phase1.size(): "
             << result.intersection_points_phase1.size() << '\n';
  GTEST_COUT << "result.intersection_prism_indices_phase1.size(): "
             << result.intersection_prism_indices_phase1.size() << '\n';

  // Verify that intersection points and indices have the same size
  ASSERT_EQ(result.intersection_points_phase0.size(),
            result.intersection_prism_indices_phase0.size());
  ASSERT_EQ(result.intersection_points_phase1.size(),
            result.intersection_prism_indices_phase1.size());

  // Verify that each intersection has exactly 5 prism indices (phase 0: 31, 36,
  // 24, 27, 58; phase 1: 51, 57, 84, 86, 23)
  for (const auto& prism_indices : result.intersection_prism_indices_phase0) {
    ASSERT_EQ(prism_indices.size(), 5);
  }
  for (const auto& prism_indices : result.intersection_prism_indices_phase1) {
    ASSERT_EQ(prism_indices.size(), 5);
  }

  // Verify that intersection points are valid (each point group has at least one
  // Vec<8> vertex)
  for (const auto& point_group : result.intersection_points_phase0) {
    ASSERT_GT(point_group.size(), 0);
  }
  for (const auto& point_group : result.intersection_points_phase1) {
    ASSERT_GT(point_group.size(), 0);
  }

  size_t areas_vertices_count_phase0 = 0;
  size_t areas_vertices_count_phase1 = 0;
  for (const auto& area : result.intersection_points_phase0) {
    areas_vertices_count_phase0 += area.size();
  }
  for (const auto& area : result.intersection_points_phase1) {
    areas_vertices_count_phase1 += area.size();
  }

  GTEST_COUT << "Phase 0 results:\n";
  GTEST_COUT << std::format("  intersection areas: {}\n",
                            result.intersection_points_phase0.size());
  GTEST_COUT << std::format("  intersection areas vertices: {}\n",
                            areas_vertices_count_phase0);
  GTEST_COUT << "Phase 1 results:\n";
  GTEST_COUT << std::format("  intersection areas: {}\n",
                            result.intersection_points_phase1.size());
  GTEST_COUT << std::format("  intersection areas vertices: {}\n",
                            areas_vertices_count_phase1);
}

TEST(user_algo, compute_areas_vertices_exhaustive) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  constexpr hcpwa::Float N = 100;
  constexpr hcpwa::Float F = 15;
  constexpr hcpwa::Float v = 0.2;
  constexpr hcpwa::Float w = 0.5;
  constexpr hcpwa::Float b51 = 0.5;
  constexpr hcpwa::Float b57 = 0.5;
  constexpr hcpwa::Float b84 = 0.5;
  constexpr hcpwa::Float b86 = 0.5;
  constexpr hcpwa::Float b31 = 0.5;
  constexpr hcpwa::Float b36 = 0.5;
  constexpr hcpwa::Float b24 = 0.5;
  constexpr hcpwa::Float b27 = 0.5;
  constexpr hcpwa::Float f2min = 5;
  constexpr hcpwa::Float f3min = 5;
  constexpr hcpwa::Float f5min = 5;
  constexpr hcpwa::Float f8min = 5;
  constexpr hcpwa::Float f2max = 10;
  constexpr hcpwa::Float f3max = 10;
  constexpr hcpwa::Float f5max = 10;
  constexpr hcpwa::Float f8max = 10;

  hcpwa::TriangleAreasVerticesResult result
      = hcpwa::compute_triangle_areas_vertices(
          N, F, v, w, b51, b57, b84, b86, b31, b36, b24, b27, f2min, f3min,
          f5min, f8min, f2max, f3max, f5max, f8max);
}
// Regression test for the 8D area assembly in compute_intersection_points().
//
// Each 8D area is the product of two 3D group cells and one 2D simplex, so its
// vertex set must be the full product of the three factors. The assembly loops
// used to be bounded by intersection_prism_indices_*, which always holds
// exactly the two prism ids that formed the cell, so every area came out with
// 2 * 2 * 3 = 12 vertices no matter how many vertices the factors really had.
//
// Deliberately built on a hand-made one-triangle-per-plane geometry: calling
// the full pipeline here would run the whole arrangement and is far too slow
// for a unit test.
TEST(user_algo, area_vertices_are_the_full_product_of_group_cells) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  constexpr hcpwa::Float N = 10;

  hcpwa::TriangleWithUniqueVertices tri;
  tri.a = {0, 0};
  tri.b = {N, 0};
  tri.c = {0, N};
  tri.a_index = 0;
  tri.b_index = 1;
  tri.c_index = 2;
  tri.polygon_index = 0;

  // Phase 0 tuple order [31, 36, 24, 27, 58], phase 1 [51, 57, 84, 86, 23].
  auto prism = [&tri](std::array<int, 2> dims) {
    return std::vector<hcpwa::LineSet<8>>{hcpwa::CalcPrism(tri, dims)};
  };
  const std::vector<hcpwa::TriangleWithUniqueVertices> tris = {tri};

  const hcpwa::PhaseIntersectionResult result = hcpwa::compute_intersection_points(
      prism({0, 2}), prism({2, 5}), prism({1, 3}), prism({1, 6}), prism({4, 7}),
      prism({0, 4}), prism({4, 6}), prism({3, 7}), prism({5, 7}), prism({1, 2}),
      tris, tris, N, /*verbose=*/false);

  // Independently rebuild each 3D group cell the way computend() does, so the
  // expected factor sizes come from the geometry rather than from the code
  // under test.
  const hcpwa::AABB<3> aabb3d = {{0, 0, 0}, {N, N, N}};
  auto group_cell_vertices = [&aabb3d](const hcpwa::LineSet<8>& p0,
                                       const hcpwa::LineSet<8>& p1,
                                       std::array<int, 3> dims) {
    hcpwa::LineSet<3> lines = hcpwa::AABBBounds(aabb3d);
    for (const auto& l : hcpwa::DimensionCast<3, 8>(p0, dims)) {
      lines.push_back(l);
    }
    for (const auto& l : hcpwa::DimensionCast<3, 8>(p1, dims)) {
      lines.push_back(l);
    }
    return hcpwa::LinesToPoints<3>(lines).size();
  };

  const std::size_t n136 = group_cell_vertices(
      hcpwa::CalcPrism(tri, {0, 2}), hcpwa::CalcPrism(tri, {2, 5}), {0, 2, 5});
  const std::size_t n247 = group_cell_vertices(
      hcpwa::CalcPrism(tri, {1, 3}), hcpwa::CalcPrism(tri, {1, 6}), {1, 3, 6});

  ASSERT_EQ(result.intersection_points_phase0.size(), 1U);
  ASSERT_GT(n136, 2U) << "degenerate fixture: the group cell must have more "
                         "than the 2 vertices the old bound would emit";
  ASSERT_GT(n247, 2U);

  EXPECT_EQ(result.intersection_points_phase0[0].size(), n136 * n247 * 3)
      << "phase-0 area is not the full product of its two group cells and the "
         "2D simplex";

  const std::size_t n157 = group_cell_vertices(
      hcpwa::CalcPrism(tri, {0, 4}), hcpwa::CalcPrism(tri, {4, 6}), {0, 4, 6});
  const std::size_t n468 = group_cell_vertices(
      hcpwa::CalcPrism(tri, {3, 7}), hcpwa::CalcPrism(tri, {5, 7}), {3, 5, 7});

  ASSERT_EQ(result.intersection_points_phase1.size(), 1U);
  EXPECT_EQ(result.intersection_points_phase1[0].size(), n157 * n468 * 3)
      << "phase-1 area is not the full product of its two group cells and the "
         "2D simplex";

  GTEST_COUT << " group cell sizes: 136=" << n136 << " 247=" << n247
             << " 157=" << n157 << " 468=" << n468 << ", phase-0 area has "
             << result.intersection_points_phase0[0].size() << " vertices\n";
}
