#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

#include "algo.hpp"
#include "cddwrap/cdd.hpp"
#include "types.hpp"
#include "utility.hpp"

// Tier 2 of the block reduction: the invariants the block factorisation of the
// geometry has to satisfy, checked without solving anything.
//
// The fixture is two triangles per projection plane, which is the smallest
// geometry that still has more than one cell per block and therefore exercises
// the region encoding. It runs in milliseconds; the same invariants on real
// arrangement geometry are checked by the pipeline test.
namespace {

constexpr hcpwa::Float kN = 10;

std::vector<hcpwa::TriangleWithUniqueVertices> squareTriangles() {
  hcpwa::TriangleWithUniqueVertices lower;
  lower.a = {0, 0};
  lower.b = {kN, 0};
  lower.c = {0, kN};
  lower.a_index = 0;
  lower.b_index = 1;
  lower.c_index = 2;
  lower.polygon_index = 0;

  hcpwa::TriangleWithUniqueVertices upper;
  upper.a = {kN, 0};
  upper.b = {kN, kN};
  upper.c = {0, kN};
  upper.a_index = 1;
  upper.b_index = 3;
  upper.c_index = 2;
  upper.polygon_index = 0;

  return {lower, upper};
}

std::vector<hcpwa::LineSet<8>> prismsFor(
    const std::vector<hcpwa::TriangleWithUniqueVertices>& triangles,
    std::array<int, 2> dims) {
  std::vector<hcpwa::LineSet<8>> out;
  out.reserve(triangles.size());
  for (const auto& triangle : triangles) {
    out.push_back(hcpwa::CalcPrism(triangle, dims));
  }
  return out;
}

hcpwa::PhaseIntersectionResult buildFixture(
    const hcpwa::TriangleGeometryOptions& options) {
  const std::vector<hcpwa::TriangleWithUniqueVertices> tris = squareTriangles();
  // Phase 0 tuple order [31, 36, 24, 27, 58], phase 1 [51, 57, 84, 86, 23].
  return hcpwa::compute_intersection_points(
      prismsFor(tris, {0, 2}), prismsFor(tris, {2, 5}), prismsFor(tris, {1, 3}),
      prismsFor(tris, {1, 6}), prismsFor(tris, {4, 7}), prismsFor(tris, {0, 4}),
      prismsFor(tris, {4, 6}), prismsFor(tris, {3, 7}), prismsFor(tris, {5, 7}),
      prismsFor(tris, {1, 2}), tris, tris, kN, /*verbose=*/false, options);
}

// j = (jA * MB + jB) * MC + jC, the order the phase loops nest in.
std::size_t encodeRegion(const std::array<hcpwa::BlockRegions, 3>& blocks,
                         std::size_t ja, std::size_t jb, std::size_t jc) {
  return (ja * blocks[1].num_regions() + jb) * blocks[2].num_regions() + jc;
}

}  // namespace

TEST(barycentric_block_geometry, blocks_carry_the_expected_axes) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;
  const hcpwa::PhaseIntersectionResult result = buildFixture({});

  // Phase 0 splits as {0,2,5} / {1,3,6} / {4,7} with planes (31,36) / (24,27) /
  // (58); phase 1 as {0,4,6} / {3,5,7} / {1,2} with (51,57) / (84,86) / (23).
  const std::array<std::vector<int>, 3> coords0
      = {std::vector<int>{0, 2, 5}, {1, 3, 6}, {4, 7}};
  const std::array<std::vector<int>, 3> coords1
      = {std::vector<int>{0, 4, 6}, {3, 5, 7}, {1, 2}};
  const std::array<std::vector<int>, 3> layers
      = {std::vector<int>{0, 1}, {2, 3}, {4}};

  for (int phase = 0; phase < 2; ++phase) {
    const auto& blocks
        = phase == 0 ? result.blocks_phase0 : result.blocks_phase1;
    const auto& expected_coords = phase == 0 ? coords0 : coords1;

    std::vector<int> seen;
    for (std::size_t b = 0; b < blocks.size(); ++b) {
      EXPECT_EQ(blocks[b].coord_count,
                static_cast<int>(expected_coords[b].size()))
          << "phase " << phase << " block " << b;
      for (int c = 0; c < blocks[b].coord_count; ++c) {
        EXPECT_EQ(blocks[b].coords[static_cast<std::size_t>(c)],
                  expected_coords[b][static_cast<std::size_t>(c)])
            << "phase " << phase << " block " << b << " coordinate " << c;
        seen.push_back(blocks[b].coords[static_cast<std::size_t>(c)]);
      }
      EXPECT_EQ(blocks[b].layer_count, static_cast<int>(layers[b].size()));
      for (int l = 0; l < blocks[b].layer_count; ++l) {
        EXPECT_EQ(blocks[b].layer_ids[static_cast<std::size_t>(l)],
                  layers[b][static_cast<std::size_t>(l)]);
      }
    }
    // The three blocks must partition the eight state coordinates.
    std::sort(seen.begin(), seen.end());
    EXPECT_EQ(seen, (std::vector<int>{0, 1, 2, 3, 4, 5, 6, 7}))
        << "phase " << phase;
  }
}

TEST(barycentric_block_geometry, region_count_is_the_product_of_block_counts) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;
  const hcpwa::PhaseIntersectionResult result = buildFixture({});

  for (int phase = 0; phase < 2; ++phase) {
    const auto& blocks
        = phase == 0 ? result.blocks_phase0 : result.blocks_phase1;
    const auto& indices = phase == 0 ? result.intersection_prism_indices_phase0
                                     : result.intersection_prism_indices_phase1;
    for (const auto& block : blocks) {
      ASSERT_GT(block.num_regions(), 0U) << "phase " << phase;
    }
    EXPECT_EQ(blocks[0].num_regions() * blocks[1].num_regions()
                  * blocks[2].num_regions(),
              indices.size())
        << "phase " << phase;
  }
}

TEST(barycentric_block_geometry, prism_ids_are_the_block_triangle_ids) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;
  const hcpwa::PhaseIntersectionResult result = buildFixture({});

  // region_triangle_ids is read straight out of these five ids, so the order
  // they are concatenated in is a contract between the geometry and the LP.
  for (int phase = 0; phase < 2; ++phase) {
    const auto& blocks
        = phase == 0 ? result.blocks_phase0 : result.blocks_phase1;
    const auto& indices = phase == 0 ? result.intersection_prism_indices_phase0
                                     : result.intersection_prism_indices_phase1;
    for (std::size_t ja = 0; ja < blocks[0].num_regions(); ++ja) {
      for (std::size_t jb = 0; jb < blocks[1].num_regions(); ++jb) {
        for (std::size_t jc = 0; jc < blocks[2].num_regions(); ++jc) {
          const std::size_t j = encodeRegion(blocks, ja, jb, jc);
          ASSERT_LT(j, indices.size());
          const std::vector<std::size_t> expected = {
              static_cast<std::size_t>(blocks[0].triangle_ids[ja][0]),
              static_cast<std::size_t>(blocks[0].triangle_ids[ja][1]),
              static_cast<std::size_t>(blocks[1].triangle_ids[jb][0]),
              static_cast<std::size_t>(blocks[1].triangle_ids[jb][1]),
              static_cast<std::size_t>(blocks[2].triangle_ids[jc][0])};
          EXPECT_EQ(indices[j], expected)
              << "phase " << phase << " region " << j;
        }
      }
    }
  }
}

TEST(barycentric_block_geometry, product_of_block_vertices_is_the_8d_region) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;
  // The 8D product is off in production; this test is the reason it still
  // exists, so it asks for it explicitly.
  const hcpwa::PhaseIntersectionResult result = buildFixture({true});

  // The one place the block path and the 8D path are compared directly. After
  // the 8D product is removed this test is what still pins the scatter rule
  //   v[coords[c]] = block_vertex[c]
  // that every phi row is built on.
  for (int phase = 0; phase < 2; ++phase) {
    const auto& blocks
        = phase == 0 ? result.blocks_phase0 : result.blocks_phase1;
    const auto& points = phase == 0 ? result.intersection_points_phase0
                                    : result.intersection_points_phase1;
    for (std::size_t ja = 0; ja < blocks[0].num_regions(); ++ja) {
      for (std::size_t jb = 0; jb < blocks[1].num_regions(); ++jb) {
        for (std::size_t jc = 0; jc < blocks[2].num_regions(); ++jc) {
          const std::size_t j = encodeRegion(blocks, ja, jb, jc);
          const auto& va = blocks[0].vertices[ja];
          const auto& vb = blocks[1].vertices[jb];
          const auto& vc = blocks[2].vertices[jc];
          ASSERT_EQ(points[j].size(), va.size() * vb.size() * vc.size())
              << "phase " << phase << " region " << j;

          std::size_t at = 0;
          for (std::size_t ka = 0; ka < va.size(); ++ka) {
            for (std::size_t kb = 0; kb < vb.size(); ++kb) {
              for (std::size_t kc = 0; kc < vc.size(); ++kc, ++at) {
                hcpwa::Vec<8> expected = hcpwa::kZeroVec;
                for (std::size_t b = 0; b < blocks.size(); ++b) {
                  const auto& source = b == 0 ? va[ka] : (b == 1 ? vb[kb]
                                                                 : vc[kc]);
                  for (int c = 0; c < blocks[b].coord_count; ++c) {
                    expected[blocks[b].coords[static_cast<std::size_t>(c)]]
                        = static_cast<hcpwa::Float>(
                            source[static_cast<std::size_t>(c)]);
                  }
                }
                for (int d = 0; d < 8; ++d) {
                  EXPECT_NEAR(static_cast<double>(points[j][at][d]),
                              static_cast<double>(expected[d]), 1e-9)
                      << "phase " << phase << " region " << j << " vertex "
                      << at << " coordinate " << d;
                }
              }
            }
          }
        }
      }
    }
  }
}

TEST(barycentric_block_geometry, block_aabb_bounds_its_own_vertices) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;
  const hcpwa::PhaseIntersectionResult result = buildFixture({});

  for (int phase = 0; phase < 2; ++phase) {
    const auto& blocks
        = phase == 0 ? result.blocks_phase0 : result.blocks_phase1;
    for (std::size_t b = 0; b < blocks.size(); ++b) {
      ASSERT_EQ(blocks[b].aabb.size(), blocks[b].vertices.size());
      for (std::size_t j = 0; j < blocks[b].vertices.size(); ++j) {
        for (const auto& vertex : blocks[b].vertices[j]) {
          for (int c = 0; c < blocks[b].coord_count; ++c) {
            const std::size_t cc = static_cast<std::size_t>(c);
            EXPECT_GE(vertex[cc], blocks[b].aabb[j].lower[cc] - 1e-12);
            EXPECT_LE(vertex[cc], blocks[b].aabb[j].upper[cc] + 1e-12);
          }
        }
      }
    }
  }
}

TEST(barycentric_block_geometry, skipping_the_8d_product_keeps_the_blocks) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;
  const hcpwa::PhaseIntersectionResult full = buildFixture({true});
  hcpwa::TriangleGeometryOptions options;
  options.build_8d_vertices = false;
  const hcpwa::PhaseIntersectionResult lean = buildFixture(options);

  // The region list and its prism ids stay; only the vertex product is gone.
  EXPECT_EQ(lean.intersection_prism_indices_phase0,
            full.intersection_prism_indices_phase0);
  EXPECT_EQ(lean.intersection_prism_indices_phase1,
            full.intersection_prism_indices_phase1);
  for (const auto& region : lean.intersection_points_phase0) {
    EXPECT_TRUE(region.empty());
  }
  for (const auto& region : lean.intersection_points_phase1) {
    EXPECT_TRUE(region.empty());
  }

  for (int phase = 0; phase < 2; ++phase) {
    const auto& a = phase == 0 ? full.blocks_phase0 : full.blocks_phase1;
    const auto& b = phase == 0 ? lean.blocks_phase0 : lean.blocks_phase1;
    for (std::size_t i = 0; i < a.size(); ++i) {
      ASSERT_EQ(a[i].num_regions(), b[i].num_regions());
      EXPECT_EQ(a[i].triangle_ids, b[i].triangle_ids);
      EXPECT_EQ(a[i].vertices, b[i].vertices);
    }
  }
}
