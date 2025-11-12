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

TEST(user_algo, compute_areas_vertices_phase0) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  constexpr hcpwa::Float N = 100;
  constexpr hcpwa::Float F = 15;
  constexpr hcpwa::Float v = 0.3;
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

  int phase = 0;

  AreasVerticesResult result = hcpwa::compute_areas_vertices(
      phase, N, F, v, w,
      b51, b57, b84, b86,
      b31, b36, b24, b27,
      f2min, f3min, f5min, f8min,
      f2max, f3max, f5max, f8max);

  // Verify that triangles are computed for all four planes
  ASSERT_GT(result.triangles51.size(), 0);
  ASSERT_GT(result.triangles57.size(), 0);
  ASSERT_GT(result.triangles84.size(), 0);
  ASSERT_GT(result.triangles86.size(), 0);

  // Verify that intersection points and indices have the same size
  ASSERT_EQ(result.intersection_points.size(), result.intersection_prism_indices.size());

  // Verify that each intersection has exactly 4 prism indices
  for (const auto& prism_indices : result.intersection_prism_indices) {
    ASSERT_EQ(prism_indices.size(), 4);
    // Verify indices are within valid ranges
    ASSERT_LT(prism_indices[0], result.triangles51.size());
    ASSERT_LT(prism_indices[1], result.triangles57.size());
    ASSERT_LT(prism_indices[2], result.triangles84.size());
    ASSERT_LT(prism_indices[3], result.triangles86.size());
  }

  // Verify that intersection points are valid (Vec<8> is a fixed-size type with 8 elements)
  for (const auto& point_group : result.intersection_points) {
    ASSERT_GT(point_group.size(), 0);
  }

  GTEST_COUT << "Phase 0 results:\n";
  GTEST_COUT << std::format("  triangles51: {}\n", result.triangles51.size());
  GTEST_COUT << std::format("  triangles57: {}\n", result.triangles57.size());
  GTEST_COUT << std::format("  triangles84: {}\n", result.triangles84.size());
  GTEST_COUT << std::format("  triangles86: {}\n", result.triangles86.size());
  GTEST_COUT << std::format("  intersections: {}\n", result.intersection_points.size());
}

TEST(user_algo, compute_areas_vertices_phase1) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  constexpr hcpwa::Float N = 100;
  constexpr hcpwa::Float F = 15;
  constexpr hcpwa::Float v = 0.3;
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

  int phase = 1;

  AreasVerticesResult result = hcpwa::compute_areas_vertices(
      phase, N, F, v, w,
      b51, b57, b84, b86,
      b31, b36, b24, b27,
      f2min, f3min, f5min, f8min,
      f2max, f3max, f5max, f8max);

  // Verify that triangles are computed for all four planes
  ASSERT_GT(result.triangles51.size(), 0);
  ASSERT_GT(result.triangles57.size(), 0);
  ASSERT_GT(result.triangles84.size(), 0);
  ASSERT_GT(result.triangles86.size(), 0);

  // Verify that intersection points and indices have the same size
  ASSERT_EQ(result.intersection_points.size(), result.intersection_prism_indices.size());

  // Verify that each intersection has exactly 4 prism indices
  for (const auto& prism_indices : result.intersection_prism_indices) {
    ASSERT_EQ(prism_indices.size(), 4);
    // Verify indices are within valid ranges
    ASSERT_LT(prism_indices[0], result.triangles51.size());
    ASSERT_LT(prism_indices[1], result.triangles57.size());
    ASSERT_LT(prism_indices[2], result.triangles84.size());
    ASSERT_LT(prism_indices[3], result.triangles86.size());
  }

  // Verify that intersection points are valid (Vec<8> is a fixed-size type with 8 elements)
  for (const auto& point_group : result.intersection_points) {
    ASSERT_GT(point_group.size(), 0);
  }

  GTEST_COUT << "Phase 1 results:\n";
  GTEST_COUT << std::format("  triangles51: {}\n", result.triangles51.size());
  GTEST_COUT << std::format("  triangles57: {}\n", result.triangles57.size());
  GTEST_COUT << std::format("  triangles84: {}\n", result.triangles84.size());
  GTEST_COUT << std::format("  triangles86: {}\n", result.triangles86.size());
  GTEST_COUT << std::format("  intersections: {}\n", result.intersection_points.size());
}

TEST(user_algo, compute_areas_vertices_phase_comparison) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  constexpr hcpwa::Float N = 100;
  constexpr hcpwa::Float F = 15;
  constexpr hcpwa::Float v = 0.3;
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

  // Test both phases and verify they produce valid results
  AreasVerticesResult result0 = hcpwa::compute_areas_vertices(
      0, N, F, v, w,
      b51, b57, b84, b86,
      b31, b36, b24, b27,
      f2min, f3min, f5min, f8min,
      f2max, f3max, f5max, f8max);

  AreasVerticesResult result1 = hcpwa::compute_areas_vertices(
      1, N, F, v, w,
      b51, b57, b84, b86,
      b31, b36, b24, b27,
      f2min, f3min, f5min, f8min,
      f2max, f3max, f5max, f8max);

  // Both phases should produce valid results
  ASSERT_GT(result0.triangles51.size(), 0);
  ASSERT_GT(result1.triangles51.size(), 0);

  // Verify that results are consistent (same structure)
  ASSERT_EQ(result0.intersection_points.size(), result0.intersection_prism_indices.size());
  ASSERT_EQ(result1.intersection_points.size(), result1.intersection_prism_indices.size());

  GTEST_COUT << "Phase comparison:\n";
  GTEST_COUT << std::format("  Phase 0 - triangles51: {}, intersections: {}\n",
                            result0.triangles51.size(), result0.intersection_points.size());
  GTEST_COUT << std::format("  Phase 1 - triangles51: {}, intersections: {}\n",
                            result1.triangles51.size(), result1.intersection_points.size());
}

