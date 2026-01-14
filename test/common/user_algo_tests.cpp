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

TEST(user_algo, compute_areas_vertices) {
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

  AreasVerticesResult result = hcpwa::compute_areas_vertices(
      N, F, v, w,
      b51, b57, b84, b86,
      b31, b36, b24, b27,
      f2min, f3min, f5min, f8min,
      f2max, f3max, f5max, f8max
  );

  // Print triangle sizes for phase 0
  GTEST_COUT << "result.triangles31.size(): " << result.triangles31.size() << '\n';
  GTEST_COUT << "result.triangles36.size(): " << result.triangles36.size() << '\n';
  GTEST_COUT << "result.triangles24.size(): " << result.triangles24.size() << '\n';
  GTEST_COUT << "result.triangles27.size(): " << result.triangles27.size() << '\n';
  // Print triangle sizes for phase 1
  GTEST_COUT << "result.triangles51.size(): " << result.triangles51.size() << '\n';
  GTEST_COUT << "result.triangles57.size(): " << result.triangles57.size() << '\n';
  GTEST_COUT << "result.triangles84.size(): " << result.triangles84.size() << '\n';
  GTEST_COUT << "result.triangles86.size(): " << result.triangles86.size() << '\n';
  // Print intersection points sizes for phase 0
  GTEST_COUT << "result.intersection_points_phase0.size(): " << result.intersection_points_phase0.size() << '\n';
  GTEST_COUT << "result.intersection_prism_indices_phase0.size(): " << result.intersection_prism_indices_phase0.size() << '\n';
  // Print intersection points sizes for phase 1
  GTEST_COUT << "result.intersection_points_phase1.size(): " << result.intersection_points_phase1.size() << '\n';
  GTEST_COUT << "result.intersection_prism_indices_phase1.size(): " << result.intersection_prism_indices_phase1.size() << '\n';

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

  // Verify that intersection points and indices have the same size for phase 1
  ASSERT_EQ(result.intersection_points_phase0.size(), result.intersection_prism_indices_phase0.size());
  ASSERT_EQ(result.intersection_points_phase1.size(), result.intersection_prism_indices_phase1.size());
  ASSERT_EQ(result.intersection_prism_indices_phase0.size(), result.intersection_prism_indices_phase1.size());
  ASSERT_EQ(result.intersection_points_phase0.size(), result.intersection_points_phase1.size());

  // Verify that each intersection has exactly 4 prism indices for phase 1
  for (const auto& prism_indices : result.intersection_prism_indices_phase1) {
    ASSERT_EQ(prism_indices.size(), 4);
    // Verify indices are within valid ranges (phase 1 uses triangles51, 57, 84, 86)
    ASSERT_LT(prism_indices[0], result.triangles51.size());
    ASSERT_LT(prism_indices[1], result.triangles57.size());
    ASSERT_LT(prism_indices[2], result.triangles84.size());
    ASSERT_LT(prism_indices[3], result.triangles86.size());
  }
  // Verify that each intersection has exactly 4 prism indices for phase 0
  for (const auto& prism_indices : result.intersection_prism_indices_phase0) {
    ASSERT_EQ(prism_indices.size(), 4);
    // Verify indices are within valid ranges (phase 0 uses triangles31, 36, 24, 27)
    ASSERT_LT(prism_indices[0], result.triangles31.size());
    ASSERT_LT(prism_indices[1], result.triangles36.size());
    ASSERT_LT(prism_indices[2], result.triangles24.size());
    ASSERT_LT(prism_indices[3], result.triangles27.size());
  }

  // Verify that intersection points are valid (Vec<8> is a fixed-size type with 8 elements)
  for (const auto& point_group : result.intersection_points_phase1) {
    ASSERT_GT(point_group.size(), 0);
  }
  for (const auto& point_group : result.intersection_points_phase0) {
    ASSERT_GT(point_group.size(), 0);
  }

  GTEST_COUT << "Phase 0 results:\n";
  GTEST_COUT << std::format("  triangles31: {}\n", result.triangles31.size());
  GTEST_COUT << std::format("  triangles36: {}\n", result.triangles36.size());
  GTEST_COUT << std::format("  triangles24: {}\n", result.triangles24.size());
  GTEST_COUT << std::format("  triangles27: {}\n", result.triangles27.size());
  GTEST_COUT << std::format("  intersections: {}\n", result.intersection_points_phase0.size());
  GTEST_COUT << "Phase 1 results:\n";
  GTEST_COUT << std::format("  triangles51: {}\n", result.triangles51.size());
  GTEST_COUT << std::format("  triangles57: {}\n", result.triangles57.size());
  GTEST_COUT << std::format("  triangles84: {}\n", result.triangles84.size());
  GTEST_COUT << std::format("  triangles86: {}\n", result.triangles86.size());
  GTEST_COUT << std::format("  intersections: {}\n", result.intersection_points_phase1.size());
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

  AreasVerticesResult result = hcpwa::compute_areas_vertices(
      N, F, v, w,
      b51, b57, b84, b86,
      b31, b36, b24, b27,
      f2min, f3min, f5min, f8min,
      f2max, f3max, f5max, f8max
  );

  
  
}
