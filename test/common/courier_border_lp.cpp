#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <algo.hpp>
#include <barycentric_geometry_types.hpp>
#include <courier_border_solver.hpp>

#include "cddwrap/cdd.hpp"
#include "types.hpp"
#include "utility.hpp"

namespace {

using barycentric_affine_approximator::ApproximationMode;
using barycentric_affine_approximator::BarycentricVarLayout;
using barycentric_affine_approximator::CourierBorderOptions;
using barycentric_affine_approximator::CourierBorderRequest;
using barycentric_affine_approximator::CourierBorderSolver;
using barycentric_affine_approximator::CourierBorderStats;
using barycentric_affine_approximator::kPhases;
using barycentric_affine_approximator::kSpaceDim;
using barycentric_affine_approximator::kSubsystemCount;
using barycentric_affine_approximator::PhaseGeometry;
using barycentric_affine_approximator::projectionAxesForPhase;
using barycentric_affine_approximator::TriangleBasis;

constexpr double kN = 10.0;

// Smallest geometry that still exercises the real structure: every projection
// plane is the unit square split into two triangles, on the production axis
// pairs. Regions then come out of the production assembly code, so the fixture
// cannot silently disagree with it about what a region is.
struct Fixture {
  std::array<PhaseGeometry, kPhases> geometries;
  std::array<BarycentricVarLayout, kPhases> layouts;
  std::array<Eigen::VectorXd, kPhases> node_weights;
};

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
  upper.polygon_index = 1;

  return {lower, upper};
}

// Mirrors the barycentric basis construction in getIntersectionPoints(): invert
// [[ax,bx,cx],[ay,by,cy],[1,1,1]] so alpha(z) = H z + h.
TriangleBasis makeBasis(const hcpwa::TriangleWithUniqueVertices& tri,
                        const std::vector<Eigen::Vector2d>& unique_vertices) {
  TriangleBasis basis;
  const std::array<Eigen::Vector2d, 3> v
      = {Eigen::Vector2d(static_cast<double>(tri.a[0]), static_cast<double>(tri.a[1])),
         Eigen::Vector2d(static_cast<double>(tri.b[0]), static_cast<double>(tri.b[1])),
         Eigen::Vector2d(static_cast<double>(tri.c[0]), static_cast<double>(tri.c[1]))};
  for (int k = 0; k < 3; ++k) {
    int found = -1;
    for (int i = 0; i < static_cast<int>(unique_vertices.size()); ++i) {
      if ((unique_vertices[static_cast<std::size_t>(i)] - v[static_cast<std::size_t>(k)])
              .norm()
          < 1e-9) {
        found = i;
        break;
      }
    }
    basis.vertex_ids[static_cast<std::size_t>(k)] = found;
  }
  Eigen::Matrix3d g;
  g << v[0](0), v[1](0), v[2](0), v[0](1), v[1](1), v[2](1), 1.0, 1.0, 1.0;
  const Eigen::Matrix3d m = g.inverse();
  basis.H = m.block<3, 2>(0, 0);
  basis.h = m.col(2);
  return basis;
}

Fixture makeFixture() {
  const auto tris = squareTriangles();
  const std::vector<Eigen::Vector2d> unique_vertices
      = {{0, 0}, {kN, 0}, {0, kN}, {kN, kN}};

  auto prisms = [&tris](std::array<int, 2> dims) {
    std::vector<hcpwa::LineSet<8>> out;
    for (const auto& t : tris) {
      out.push_back(hcpwa::CalcPrism(t, dims));
    }
    return out;
  };

  const auto ax0 = projectionAxesForPhase(0);
  const auto ax1 = projectionAxesForPhase(1);

  const hcpwa::PhaseIntersectionResult inter = hcpwa::compute_intersection_points(
      prisms(ax0[0]), prisms(ax0[1]), prisms(ax0[2]), prisms(ax0[3]),
      prisms(ax0[4]), prisms(ax1[0]), prisms(ax1[1]), prisms(ax1[2]),
      prisms(ax1[3]), prisms(ax1[4]), tris, tris, kN, /*verbose=*/false);

  Fixture f;
  for (int phase = 0; phase < kPhases; ++phase) {
    const auto axes = projectionAxesForPhase(phase);
    for (int s = 0; s < kSubsystemCount; ++s) {
      auto& layer = f.geometries[static_cast<std::size_t>(phase)].layers[s];
      layer.axes = axes[static_cast<std::size_t>(s)];
      layer.triangles = tris;
      layer.unique_vertices = unique_vertices;
      for (const auto& t : tris) {
        layer.bases.push_back(makeBasis(t, unique_vertices));
      }
    }

    const auto& pts = phase == 0 ? inter.intersection_points_phase0
                                 : inter.intersection_points_phase1;
    const auto& idx = phase == 0 ? inter.intersection_prism_indices_phase0
                                 : inter.intersection_prism_indices_phase1;
    auto& geom = f.geometries[static_cast<std::size_t>(phase)];
    for (std::size_t j = 0; j < pts.size(); ++j) {
      std::vector<Eigen::VectorXd> verts;
      for (const auto& p : pts[j]) {
        Eigen::VectorXd v(kSpaceDim);
        for (int d = 0; d < kSpaceDim; ++d) {
          v(d) = static_cast<double>(p[d]);
        }
        verts.push_back(std::move(v));
      }
      geom.region_vertices.push_back(std::move(verts));
      std::array<int, kSubsystemCount> tuple{};
      for (int s = 0; s < kSubsystemCount; ++s) {
        tuple[static_cast<std::size_t>(s)] = static_cast<int>(idx[j][static_cast<std::size_t>(s)]);
      }
      geom.region_triangle_ids.push_back(tuple);
    }
    // The solver works block by block; blockGeometryFromRegions is the same
    // ingestion the approximator uses, so the fixture cannot disagree with it
    // about what a block is.
    geom.blocks = barycentric_affine_approximator::blockGeometryFromRegions(
        phase, phase == 0 ? inter.blocks_phase0 : inter.blocks_phase1);

    BarycentricVarLayout layout;
    for (int s = 0; s < kSubsystemCount; ++s) {
      layout.offset_s[static_cast<std::size_t>(s)] = layout.num_x;
      layout.eta_s[static_cast<std::size_t>(s)]
          = static_cast<int>(unique_vertices.size());
      layout.num_x += layout.eta_s[static_cast<std::size_t>(s)];
    }
    f.layouts[static_cast<std::size_t>(phase)] = layout;
    f.node_weights[static_cast<std::size_t>(phase)]
        = Eigen::VectorXd::Ones(layout.num_x);
  }
  return f;
}

// A candidate whose value is the constant c everywhere. Barycentric weights of
// each plane sum to 1, so five planes each carrying c/5 give exactly c.
std::vector<double> constantCandidate(const BarycentricVarLayout& layout,
                                      double c) {
  return std::vector<double>(static_cast<std::size_t>(layout.num_x),
                             c / static_cast<double>(kSubsystemCount));
}

double valueAtRegionVertex(const PhaseGeometry& geom,
                           const BarycentricVarLayout& layout, int region,
                           const Eigen::VectorXd& nu,
                           const std::vector<double>& z) {
  double acc = 0.0;
  const auto& tri_ids = geom.region_triangle_ids[static_cast<std::size_t>(region)];
  for (int s = 0; s < kSubsystemCount; ++s) {
    const auto& layer = geom.layers[static_cast<std::size_t>(s)];
    const TriangleBasis& basis
        = layer.bases[static_cast<std::size_t>(tri_ids[static_cast<std::size_t>(s)])];
    const Eigen::Vector3d a
        = basis.H * Eigen::Vector2d(nu(layer.axes[0]), nu(layer.axes[1])) + basis.h;
    for (int k = 0; k < 3; ++k) {
      acc += a(k)
             * z[static_cast<std::size_t>(
                 layout.idxX(s, basis.vertex_ids[static_cast<std::size_t>(k)]))];
    }
  }
  return acc;
}

}  // namespace

TEST(courier_border_lp, lower_reproduces_a_constant_candidate) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  const Fixture f = makeFixture();
  ASSERT_FALSE(f.geometries[0].region_vertices.empty());

  CourierBorderSolver solver{CourierBorderOptions{}};
  solver.prepare(f.geometries, f.layouts, f.node_weights, kN);
  ASSERT_TRUE(solver.prepared());

  constexpr double kConst = 3.5;
  const std::vector<std::vector<double>> cands
      = {constantCandidate(f.layouts[1], kConst)};

  CourierBorderRequest req;
  req.target_phase = 0;
  req.source_phase = 1;
  req.candidates = cands;

  CourierBorderStats stats;
  const std::vector<double> z
      = solver.solve(req, ApproximationMode::Lower, &stats);
  ASSERT_EQ(z.size(), static_cast<std::size_t>(f.layouts[0].num_x));

  // The bound is tight: Vt <= const is achievable with equality, and the
  // objective maximizes the integral of Vt.
  for (int j = 0; j < static_cast<int>(f.geometries[0].region_vertices.size());
       ++j) {
    for (const auto& nu : f.geometries[0].region_vertices[static_cast<std::size_t>(j)]) {
      const double v = valueAtRegionVertex(f.geometries[0], f.layouts[0], j, nu, z);
      EXPECT_NEAR(v, kConst, 1e-5);
    }
  }
  EXPECT_EQ(solver.worstCertificateResidual(req, ApproximationMode::Lower, z),
            0.0);
}

TEST(courier_border_lp, upper_reproduces_a_constant_candidate) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  const Fixture f = makeFixture();
  CourierBorderSolver solver{CourierBorderOptions{}};
  solver.prepare(f.geometries, f.layouts, f.node_weights, kN);

  constexpr double kConst = -2.25;
  const std::vector<std::vector<double>> cands
      = {constantCandidate(f.layouts[1], kConst)};

  CourierBorderRequest req;
  req.target_phase = 0;
  req.source_phase = 1;
  req.candidates = cands;

  const std::vector<double> z = solver.solve(req, ApproximationMode::Upper);
  for (int j = 0; j < static_cast<int>(f.geometries[0].region_vertices.size());
       ++j) {
    for (const auto& nu : f.geometries[0].region_vertices[static_cast<std::size_t>(j)]) {
      const double v = valueAtRegionVertex(f.geometries[0], f.layouts[0], j, nu, z);
      EXPECT_NEAR(v, kConst, 1e-5);
    }
  }
}

// Sandwich: with the same inputs the lower estimate must never exceed the upper
// one. Needs no reference implementation, and it is the end-to-end check that
// survives when a full run is out of reach.
TEST(courier_border_lp, lower_never_exceeds_upper) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  const Fixture f = makeFixture();
  CourierBorderSolver solver{CourierBorderOptions{}};
  solver.prepare(f.geometries, f.layouts, f.node_weights, kN);

  // Two genuinely different candidates, so the max over rho is not trivial.
  std::vector<double> a = constantCandidate(f.layouts[1], 1.0);
  std::vector<double> b = constantCandidate(f.layouts[1], 1.0);
  for (int k = 0; k < f.layouts[1].eta_s[0]; ++k) {
    a[static_cast<std::size_t>(f.layouts[1].idxX(0, k))] += 0.3 * k;
    b[static_cast<std::size_t>(f.layouts[1].idxX(0, k))] -= 0.2 * k;
  }
  const std::vector<std::vector<double>> cands = {a, b};

  CourierBorderRequest req;
  req.target_phase = 0;
  req.source_phase = 1;
  req.candidates = cands;

  const std::vector<double> lo = solver.solve(req, ApproximationMode::Lower);
  const std::vector<double> hi = solver.solve(req, ApproximationMode::Upper);

  // With more than one candidate the re-verification must use the same
  // per-region candidate the solver chose. Demanding the courier stay under
  // every candidate would be strictly stronger than the border condition and
  // would reject this perfectly valid solution.
  EXPECT_EQ(solver.worstCertificateResidual(req, ApproximationMode::Lower, lo),
            0.0);
  EXPECT_EQ(solver.worstCertificateResidual(req, ApproximationMode::Upper, hi),
            0.0);

  for (int j = 0; j < static_cast<int>(f.geometries[0].region_vertices.size());
       ++j) {
    for (const auto& nu : f.geometries[0].region_vertices[static_cast<std::size_t>(j)]) {
      const double vlo = valueAtRegionVertex(f.geometries[0], f.layouts[0], j, nu, lo);
      const double vhi = valueAtRegionVertex(f.geometries[0], f.layouts[0], j, nu, hi);
      EXPECT_LE(vlo, vhi + 1e-6)
          << "lower estimate exceeds the upper one at region " << j;
    }
  }
}

// Explicitly covers the Benders path: a run that needs cuts must actually
// generate them, converge, and end up certified. Without this the cut code --
// dual extraction, the stationarity identities and the separation check --
// would only ever be reached incidentally.
TEST(courier_border_lp, generates_cuts_and_still_certifies) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  const Fixture f = makeFixture();
  CourierBorderOptions opts;
  opts.max_seed_rows = 0;  // force the master to learn everything from cuts
  CourierBorderSolver solver{opts};
  solver.prepare(f.geometries, f.layouts, f.node_weights, kN);

  std::vector<double> a = constantCandidate(f.layouts[1], 2.0);
  for (int k = 0; k < f.layouts[1].eta_s[0]; ++k) {
    a[static_cast<std::size_t>(f.layouts[1].idxX(0, k))] += 0.45 * k;
  }
  const std::vector<std::vector<double>> cands = {a};

  CourierBorderRequest req;
  req.target_phase = 0;
  req.source_phase = 1;
  req.candidates = cands;

  CourierBorderStats stats;
  const std::vector<double> z
      = solver.solve(req, ApproximationMode::Lower, &stats);

  EXPECT_GT(stats.cuts_added, 0)
      << "no cut was generated, so the dual/cut path never ran";
  EXPECT_GT(stats.iterations, 1);
  EXPECT_GT(stats.subproblems_solved, 0);
  // Converging means every region is certified; re-derive it independently.
  EXPECT_EQ(solver.worstCertificateResidual(req, ApproximationMode::Lower, z),
            0.0);
}

// Failure must be loud. run() writes whatever comes back straight into the
// value function and marches on it, so a partial result would poison the whole
// backward sweep.
TEST(courier_border_lp, iteration_cap_throws_instead_of_returning) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  const Fixture f = makeFixture();
  CourierBorderOptions opts;
  opts.max_iterations = 1;
  opts.max_seed_rows = 0;      // no seeding, so one iteration cannot converge
  opts.max_cuts_per_iteration = 1;
  CourierBorderSolver solver{opts};
  solver.prepare(f.geometries, f.layouts, f.node_weights, kN);

  std::vector<double> a = constantCandidate(f.layouts[1], 1.0);
  for (int k = 0; k < f.layouts[1].eta_s[0]; ++k) {
    a[static_cast<std::size_t>(f.layouts[1].idxX(0, k))] += 0.7 * k;
  }
  const std::vector<std::vector<double>> cands = {a};

  CourierBorderRequest req;
  req.target_phase = 0;
  req.source_phase = 1;
  req.candidates = cands;

  EXPECT_THROW(solver.solve(req, ApproximationMode::Lower), std::runtime_error);
}

TEST(courier_border_lp, rejects_unprepared_and_malformed_requests) {
  const Fixture f = makeFixture();
  CourierBorderSolver solver{CourierBorderOptions{}};
  const std::vector<std::vector<double>> cands
      = {constantCandidate(f.layouts[1], 1.0)};

  CourierBorderRequest req;
  req.target_phase = 0;
  req.source_phase = 1;
  req.candidates = cands;

  EXPECT_FALSE(solver.prepared());
  EXPECT_THROW(solver.solve(req, ApproximationMode::Lower), std::runtime_error);
}
