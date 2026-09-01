#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <random>
#include <vector>

#include <algo.hpp>
#include <barycentric_geometry_types.hpp>
#include <courier_border_solver.hpp>

#include "cddwrap/cdd.hpp"
#include "morph.hpp"
#include "types.hpp"
#include "utility.hpp"

namespace {

using barycentric_affine_approximator::kSubsystemCount;
using barycentric_affine_approximator::projectionAxesForPhase;
using barycentric_affine_approximator::detail::clipTriangleToRect;
using barycentric_affine_approximator::detail::Point2;

// Coordinate groups of a phase are the connected components of the graph whose
// edges are that phase's projection planes. Derived here rather than hard-coded
// so the test tracks projectionAxesForPhase() instead of duplicating it.
std::array<int, 8> groupOf(int phase) {
  const auto axes = projectionAxesForPhase(phase);
  std::array<int, 8> parent{};
  for (int i = 0; i < 8; ++i) {
    parent[static_cast<std::size_t>(i)] = i;
  }
  std::function<int(int)> find = [&](int x) {
    while (parent[static_cast<std::size_t>(x)] != x) {
      parent[static_cast<std::size_t>(x)]
          = parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(x)])];
      x = parent[static_cast<std::size_t>(x)];
    }
    return x;
  };
  for (const auto& e : axes) {
    parent[static_cast<std::size_t>(find(e[0]))] = find(e[1]);
  }
  std::array<int, 8> group{};
  for (int i = 0; i < 8; ++i) {
    group[static_cast<std::size_t>(i)] = find(i);
  }
  return group;
}

bool sameSet(std::vector<Point2> a, std::vector<Point2> b, double tol) {
  if (a.size() != b.size()) {
    return false;
  }
  for (const auto& p : a) {
    const bool found = std::any_of(b.begin(), b.end(), [&](const Point2& q) {
      return std::abs(p.x - q.x) <= tol && std::abs(p.y - q.y) <= tol;
    });
    if (!found) {
      return false;
    }
  }
  return true;
}

// Reference clip through the production vertex enumerator, so the hand-written
// Sutherland-Hodgman is checked against cdd rather than against itself.
// Deliberately not used in the solver: it throws on an empty intersection and
// needs cdd's global state, and prepare() performs hundreds of thousands of
// clips.
std::vector<Point2> clipViaCdd(const std::array<Point2, 3>& tri, double lo_x,
                               double hi_x, double lo_y, double hi_y) {
  hcpwa::Triangle t;
  t.a = {static_cast<hcpwa::Float>(tri[0].x), static_cast<hcpwa::Float>(tri[0].y)};
  t.b = {static_cast<hcpwa::Float>(tri[1].x), static_cast<hcpwa::Float>(tri[1].y)};
  t.c = {static_cast<hcpwa::Float>(tri[2].x), static_cast<hcpwa::Float>(tri[2].y)};

  const hcpwa::AABB<2> box = {{static_cast<hcpwa::Float>(lo_x),
                               static_cast<hcpwa::Float>(lo_y)},
                              {static_cast<hcpwa::Float>(hi_x),
                               static_cast<hcpwa::Float>(hi_y)}};
  hcpwa::LineSet<2> lines = hcpwa::AABBBounds(box);
  for (const auto& l :
       hcpwa::DimensionCast<2, 8>(hcpwa::CalcPrism(t, {0, 1}), {0, 1})) {
    lines.push_back(l);
  }
  std::vector<Point2> out;
  try {
    for (const auto& p : hcpwa::LinesToPoints<2>(lines)) {
      out.push_back(Point2{static_cast<double>(p[0]), static_cast<double>(p[1])});
    }
  } catch (const std::exception&) {
    return {};  // cdd reports an empty intersection by throwing
  }
  return out;
}

}  // namespace

// Every plane of one phase must draw its two axes from two DIFFERENT coordinate
// groups of the other phase. This is what makes the projection of a region onto
// a source plane an axis-aligned rectangle, which is the whole reason the
// courier certificate can be checked plane by plane on a finite vertex set.
TEST(courier_border_geometry, source_planes_split_target_groups) {
  for (int target = 0; target < 2; ++target) {
    const std::array<int, 8> group = groupOf(target);
    const auto source_axes = projectionAxesForPhase(1 - target);
    for (int s = 0; s < kSubsystemCount; ++s) {
      const int ga = group[static_cast<std::size_t>(source_axes[s][0])];
      const int gb = group[static_cast<std::size_t>(source_axes[s][1])];
      EXPECT_NE(ga, gb)
          << "target phase " << target << ", source plane " << s << " (axes "
          << source_axes[s][0] << "," << source_axes[s][1]
          << ") lies inside a single target group, so its projection is not a "
             "rectangle and the courier certificate is unsound";
    }
  }
}

// The mirror property: every plane of a phase lies wholly inside one of that
// phase's own groups. True by construction -- the groups are the components of
// exactly this graph -- but it is what makes the value function separate over
// the three groups, so a future edit that broke it should fail loudly here.
TEST(courier_border_geometry, target_planes_lie_inside_one_group) {
  for (int phase = 0; phase < 2; ++phase) {
    const std::array<int, 8> group = groupOf(phase);
    const auto axes = projectionAxesForPhase(phase);
    for (int s = 0; s < kSubsystemCount; ++s) {
      EXPECT_EQ(group[static_cast<std::size_t>(axes[s][0])],
                group[static_cast<std::size_t>(axes[s][1])])
          << "phase " << phase << " plane " << s;
    }
  }
  // Both phases must split the 8 coordinates into exactly three groups.
  for (int phase = 0; phase < 2; ++phase) {
    const std::array<int, 8> group = groupOf(phase);
    std::vector<int> roots(group.begin(), group.end());
    std::sort(roots.begin(), roots.end());
    roots.erase(std::unique(roots.begin(), roots.end()), roots.end());
    EXPECT_EQ(roots.size(), 3U) << "phase " << phase;
  }
}

TEST(courier_border_geometry, clip_triangle_fully_inside_is_unchanged) {
  const std::array<Point2, 3> tri = {Point2{2, 2}, Point2{6, 3}, Point2{3, 7}};
  const auto out = clipTriangleToRect(tri, 0, 10, 0, 10);
  ASSERT_EQ(out.size(), 3U);
  EXPECT_TRUE(sameSet(out, {tri[0], tri[1], tri[2]}, 1e-9));
}

TEST(courier_border_geometry, clip_triangle_fully_outside_is_empty) {
  const std::array<Point2, 3> tri
      = {Point2{20, 20}, Point2{25, 20}, Point2{20, 25}};
  EXPECT_TRUE(clipTriangleToRect(tri, 0, 10, 0, 10).empty());
}

TEST(courier_border_geometry, clip_triangle_straddling_an_edge) {
  // Right triangle with legs 8, clipped at x <= 4: the part kept is the
  // trapezoid (0,0), (4,0), (4,4), (0,8).
  const std::array<Point2, 3> tri = {Point2{0, 0}, Point2{8, 0}, Point2{0, 8}};
  const auto out = clipTriangleToRect(tri, 0, 4, 0, 8);
  EXPECT_TRUE(sameSet(out,
                      {Point2{0, 0}, Point2{4, 0}, Point2{4, 4}, Point2{0, 8}},
                      1e-9))
      << "got " << out.size() << " vertices";
}

TEST(courier_border_geometry, clip_keeps_a_vertex_lying_on_the_edge) {
  // The apex sits exactly on x = 4. Dropping it would drop an LP row.
  const std::array<Point2, 3> tri = {Point2{0, 0}, Point2{4, 2}, Point2{0, 4}};
  const auto out = clipTriangleToRect(tri, 0, 4, 0, 4);
  const bool has_apex = std::any_of(out.begin(), out.end(), [](const Point2& p) {
    return std::abs(p.x - 4.0) <= 1e-9 && std::abs(p.y - 2.0) <= 1e-9;
  });
  EXPECT_TRUE(has_apex);
}

TEST(courier_border_geometry, clip_keeps_a_degenerate_triangle) {
  // Zero-area triangle: it must collapse to points, never vanish, because a
  // dropped constraint can break soundness.
  const std::array<Point2, 3> tri = {Point2{1, 1}, Point2{3, 3}, Point2{5, 5}};
  const auto out = clipTriangleToRect(tri, 0, 10, 0, 10);
  EXPECT_FALSE(out.empty());
}

// Cross-validation against the production vertex enumerator on random cases.
TEST(courier_border_geometry, clip_agrees_with_cdd_vertex_enumeration) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  std::mt19937 rng(12345);
  std::uniform_real_distribution<double> coord(0.0, 10.0);
  int compared = 0;
  for (int trial = 0; trial < 200; ++trial) {
    const std::array<Point2, 3> tri = {Point2{coord(rng), coord(rng)},
                                       Point2{coord(rng), coord(rng)},
                                       Point2{coord(rng), coord(rng)}};
    double x0 = coord(rng);
    double x1 = coord(rng);
    double y0 = coord(rng);
    double y1 = coord(rng);
    if (x0 > x1) {
      std::swap(x0, x1);
    }
    if (y0 > y1) {
      std::swap(y0, y1);
    }
    if (x1 - x0 < 0.5 || y1 - y0 < 0.5) {
      continue;
    }
    // Skip near-degenerate triangles: cdd is entitled to differ there, and the
    // solver keeps slivers on purpose while cdd may reject them.
    const double area = std::abs((tri[1].x - tri[0].x) * (tri[2].y - tri[0].y)
                                 - (tri[2].x - tri[0].x) * (tri[1].y - tri[0].y));
    if (area < 1.0) {
      continue;
    }

    const auto mine = clipTriangleToRect(tri, x0, x1, y0, y1);
    const auto reference = clipViaCdd(tri, x0, x1, y0, y1);
    if (reference.empty() && mine.empty()) {
      ++compared;
      continue;
    }
    // A sliver of near-zero area is where the two legitimately disagree.
    if (mine.size() < 3 || reference.size() < 3) {
      continue;
    }
    EXPECT_TRUE(sameSet(mine, reference, 1e-6))
        << "trial " << trial << ": clip produced " << mine.size()
        << " vertices, cdd produced " << reference.size();
    ++compared;
  }
  EXPECT_GT(compared, 50) << "too few usable random cases to trust the check";
}
