#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "algo.hpp"
#include "barycentric_geometry_types.hpp"
#include "cddwrap/cdd.hpp"
#include "types.hpp"
#include "utility.hpp"

// Tier 2, approximator side: the barycentric charts of the block factorisation.
//
// local_axis -- where a plane's two axes sit inside its block's coordinate list
// -- is the one table whose mistakes are silent: a wrong entry projects the
// vertex somewhere else on the plane and produces a wrong phi row rather than a
// crash. It is derived by search in production and pinned here.
namespace {

namespace baa = barycentric_affine_approximator;

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

hcpwa::PhaseIntersectionResult buildFixture() {
  const std::vector<hcpwa::TriangleWithUniqueVertices> tris = squareTriangles();
  const auto ax0 = baa::projectionAxesForPhase(0);
  const auto ax1 = baa::projectionAxesForPhase(1);
  return hcpwa::compute_intersection_points(
      prismsFor(tris, ax0[0]), prismsFor(tris, ax0[1]), prismsFor(tris, ax0[2]),
      prismsFor(tris, ax0[3]), prismsFor(tris, ax0[4]), prismsFor(tris, ax1[0]),
      prismsFor(tris, ax1[1]), prismsFor(tris, ax1[2]), prismsFor(tris, ax1[3]),
      prismsFor(tris, ax1[4]), tris, tris, kN, /*verbose=*/false,
      // The product path is what this file compares the block path against.
      {/*build_8d_vertices=*/true});
}

// Mirrors the basis construction in getIntersectionPoints(): invert
// [[ax,bx,cx],[ay,by,cy],[1,1,1]] so that alpha(z) = H z + h.
baa::TriangleBasis makeBasis(const hcpwa::TriangleWithUniqueVertices& tri,
                             const std::vector<Eigen::Vector2d>& unique) {
  baa::TriangleBasis basis;
  const std::array<Eigen::Vector2d, 3> v = {
      Eigen::Vector2d(static_cast<double>(tri.a[0]),
                      static_cast<double>(tri.a[1])),
      Eigen::Vector2d(static_cast<double>(tri.b[0]),
                      static_cast<double>(tri.b[1])),
      Eigen::Vector2d(static_cast<double>(tri.c[0]),
                      static_cast<double>(tri.c[1]))};
  for (int k = 0; k < 3; ++k) {
    basis.vertex_ids[static_cast<std::size_t>(k)] = -1;
    for (int i = 0; i < static_cast<int>(unique.size()); ++i) {
      if ((unique[static_cast<std::size_t>(i)]
           - v[static_cast<std::size_t>(k)]).norm() < 1e-9) {
        basis.vertex_ids[static_cast<std::size_t>(k)] = i;
        break;
      }
    }
  }
  Eigen::Matrix3d g;
  g << v[0](0), v[1](0), v[2](0), v[0](1), v[1](1), v[2](1), 1.0, 1.0, 1.0;
  const Eigen::Matrix3d m = g.inverse();
  basis.H = m.block<3, 2>(0, 0);
  basis.h = m.col(2);
  return basis;
}

struct LpFixture {
  std::array<baa::PhaseGeometry, baa::kPhases> geometries;
  std::array<baa::BarycentricVarLayout, baa::kPhases> layouts;
};

LpFixture makeLpFixture() {
  const auto tris = squareTriangles();
  const std::vector<Eigen::Vector2d> unique = {{0, 0}, {kN, 0}, {0, kN},
                                               {kN, kN}};
  const hcpwa::PhaseIntersectionResult inter = buildFixture();

  LpFixture f;
  for (int phase = 0; phase < baa::kPhases; ++phase) {
    const auto axes = baa::projectionAxesForPhase(phase);
    auto& geom = f.geometries[static_cast<std::size_t>(phase)];
    for (int s = 0; s < baa::kSubsystemCount; ++s) {
      auto& layer = geom.layers[static_cast<std::size_t>(s)];
      layer.axes = axes[static_cast<std::size_t>(s)];
      layer.triangles = tris;
      layer.unique_vertices = unique;
      for (const auto& t : tris) {
        layer.bases.push_back(makeBasis(t, unique));
      }
    }

    const auto& pts = phase == 0 ? inter.intersection_points_phase0
                                 : inter.intersection_points_phase1;
    const auto& idx = phase == 0 ? inter.intersection_prism_indices_phase0
                                 : inter.intersection_prism_indices_phase1;
    for (std::size_t j = 0; j < pts.size(); ++j) {
      std::vector<Eigen::VectorXd> verts;
      for (const auto& p : pts[j]) {
        Eigen::VectorXd v(baa::kSpaceDim);
        for (int d = 0; d < baa::kSpaceDim; ++d) {
          v(d) = static_cast<double>(p[d]);
        }
        verts.push_back(std::move(v));
      }
      geom.region_vertices.push_back(std::move(verts));
      std::array<int, baa::kSubsystemCount> tuple{};
      for (int s = 0; s < baa::kSubsystemCount; ++s) {
        tuple[static_cast<std::size_t>(s)]
            = static_cast<int>(idx[j][static_cast<std::size_t>(s)]);
      }
      geom.region_triangle_ids.push_back(tuple);
    }
    geom.blocks = baa::blockGeometryFromRegions(
        phase, phase == 0 ? inter.blocks_phase0 : inter.blocks_phase1);

    baa::BarycentricVarLayout layout;
    for (int s = 0; s < baa::kSubsystemCount; ++s) {
      layout.offset_s[static_cast<std::size_t>(s)] = layout.num_x;
      layout.eta_s[static_cast<std::size_t>(s)]
          = static_cast<int>(unique.size());
      layout.num_x += layout.eta_s[static_cast<std::size_t>(s)];
    }
    f.layouts[static_cast<std::size_t>(phase)] = layout;
  }
  return f;
}

}  // namespace

TEST(barycentric_block_charts, local_axes_match_the_pinned_table) {
  // The six tuples, written down once. Two of them are non-adjacent -- phase 0
  // layer 3 is (0,2) and phase 1 layer 2 is (0,2) -- and a wrong entry produces
  // a silently wrong phi row rather than a crash, so they are pinned here as
  // well as derived by search in production.
  const std::array<std::vector<int>, 3> coords0
      = {std::vector<int>{0, 2, 5}, {1, 3, 6}, {4, 7}};
  const std::array<std::vector<int>, 3> coords1
      = {std::vector<int>{0, 4, 6}, {3, 5, 7}, {1, 2}};
  const std::array<std::vector<int>, 3> layers
      = {std::vector<int>{0, 1}, {2, 3}, {4}};

  const std::array<std::array<std::array<int, 2>, 2>, 3> expected0
      = {{{{{0, 1}, {1, 2}}}, {{{0, 1}, {0, 2}}}, {{{0, 1}, {0, 0}}}}};
  const std::array<std::array<std::array<int, 2>, 2>, 3> expected1
      = {{{{{0, 1}, {1, 2}}}, {{{0, 2}, {1, 2}}}, {{{0, 1}, {0, 0}}}}};

  for (int phase = 0; phase < 2; ++phase) {
    const auto& coords = phase == 0 ? coords0 : coords1;
    const auto& expected = phase == 0 ? expected0 : expected1;
    for (std::size_t b = 0; b < 3; ++b) {
      std::array<int, baa::kMaxBlockCoords> c{};
      for (std::size_t i = 0; i < coords[b].size(); ++i) {
        c[i] = coords[b][i];
      }
      std::array<int, baa::kMaxBlockLayers> l{};
      for (std::size_t i = 0; i < layers[b].size(); ++i) {
        l[i] = layers[b][i];
      }
      const auto local = baa::localAxesForBlock(
          phase, c, static_cast<int>(coords[b].size()), l,
          static_cast<int>(layers[b].size()));
      for (std::size_t layer = 0; layer < layers[b].size(); ++layer) {
        EXPECT_EQ(local[layer][0], expected[b][layer][0])
            << "phase " << phase << " block " << b << " layer " << layer;
        EXPECT_EQ(local[layer][1], expected[b][layer][1])
            << "phase " << phase << " block " << b << " layer " << layer;
      }
    }
  }
}

TEST(barycentric_block_charts, coordinate_blocks_are_the_graph_components) {
  const std::array<std::vector<std::vector<int>>, 2> expected = {
      std::vector<std::vector<int>>{{0, 2, 5}, {1, 3, 6}, {4, 7}},
      std::vector<std::vector<int>>{{0, 4, 6}, {1, 2}, {3, 5, 7}}};
  for (int phase = 0; phase < 2; ++phase) {
    const auto blocks = baa::coordinateBlocksForPhase(phase);
    std::vector<std::vector<int>> got(blocks.begin(), blocks.end());
    std::sort(got.begin(), got.end());
    std::vector<std::vector<int>> want
        = expected[static_cast<std::size_t>(phase)];
    std::sort(want.begin(), want.end());
    EXPECT_EQ(got, want) << "phase " << phase;
  }
}

TEST(barycentric_block_charts, block_phi_rows_sum_to_the_product_phi_row) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;
  const LpFixture f = makeLpFixture();

  for (int phase = 0; phase < baa::kPhases; ++phase) {
    const auto& geom = f.geometries[static_cast<std::size_t>(phase)];
    const auto& layout = f.layouts[static_cast<std::size_t>(phase)];
    const auto& blocks = geom.blocks;

    for (std::size_t ja = 0;
         ja < static_cast<std::size_t>(blocks[0].numRegions()); ++ja) {
      for (std::size_t jb = 0;
           jb < static_cast<std::size_t>(blocks[1].numRegions()); ++jb) {
        for (std::size_t jc = 0;
             jc < static_cast<std::size_t>(blocks[2].numRegions()); ++jc) {
          const std::size_t j
              = (ja * static_cast<std::size_t>(blocks[1].numRegions()) + jb)
                    * static_cast<std::size_t>(blocks[2].numRegions())
                + jc;
          const std::array<std::size_t, 3> jbs = {ja, jb, jc};

          // Take a vertex of each block and both the product point it builds
          // and the three block rows it produces.
          for (std::size_t ka = 0; ka < blocks[0].vertices[ja].size(); ++ka) {
            const std::array<std::size_t, 3> ks = {ka, 0, 0};
            Eigen::VectorXd point = Eigen::VectorXd::Zero(baa::kSpaceDim);
            std::vector<double> accumulated(
                static_cast<std::size_t>(layout.num_x), 0.0);
            double mass_total = 0.0;

            for (std::size_t b = 0; b < 3; ++b) {
              const auto& vertex
                  = blocks[b].vertices[jbs[b]][ks[b] % blocks[b]
                                                   .vertices[jbs[b]].size()];
              for (int c = 0; c < blocks[b].coord_count; ++c) {
                point(blocks[b].coords[static_cast<std::size_t>(c)])
                    = vertex(c);
              }
              const baa::SparseVec row = baa::buildPhiRowBlock(
                  geom, layout, static_cast<int>(b),
                  static_cast<int>(jbs[b]), vertex);
              double mass = 0.0;
              for (std::size_t i = 0; i < row.cols.size(); ++i) {
                accumulated[static_cast<std::size_t>(row.cols[i])]
                    += row.vals[i];
                mass += row.vals[i];
              }
              // One unit of barycentric mass per plane of the block.
              EXPECT_NEAR(mass, static_cast<double>(blocks[b].layer_count),
                          1e-9)
                  << "phase " << phase << " block " << b;
              mass_total += mass;
            }
            EXPECT_NEAR(mass_total, 5.0, 1e-9) << "phase " << phase;

            const baa::SparseVec product
                = baa::buildPhiRow(geom, layout, static_cast<int>(j), point);
            std::vector<double> expected(
                static_cast<std::size_t>(layout.num_x), 0.0);
            for (std::size_t i = 0; i < product.cols.size(); ++i) {
              expected[static_cast<std::size_t>(product.cols[i])]
                  = product.vals[i];
            }
            for (int c = 0; c < layout.num_x; ++c) {
              EXPECT_NEAR(accumulated[static_cast<std::size_t>(c)],
                          expected[static_cast<std::size_t>(c)], 1e-9)
                  << "phase " << phase << " region " << j << " column " << c;
            }
          }
        }
      }
    }
  }
}
