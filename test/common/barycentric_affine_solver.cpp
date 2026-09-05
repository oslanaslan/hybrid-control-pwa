#include <gtest/gtest.h>

#include <cstdlib>
#include <hcpwa.hpp>
#include <algo.hpp>
#include "cddwrap/cdd.hpp"
#include <uniqie_pool.hpp>
#include <symbolic.hpp>
#include <Eigen/Dense>
#include <barycentric_affine_approximator.hpp>
#include "utility.hpp"

#include "morph.hpp"
#include "types.hpp"

// The caller owns the cdd global constants: they have to outlive the whole
// run, not just this factory. They used to be initialised and freed here, so
// the arrangement -- which goes through cdd for every vertex enumeration --
// ran after dd_free_global_constants(). Harmless in the double build, where
// clearing a constant is a no-op, and silent in any build until it is not.
static barycentric_affine_approximator::BarycentricAffineApproximator
create_barycentric_approximator() {
  // Same parameters as user_algo_tests (compute_areas_vertices)
  constexpr double N = 160;
  constexpr double F = 0.5;
  constexpr double v = 0.017;
  constexpr double w = 0.0055;
  constexpr double b51 = 0.6;
  constexpr double b57 = 0.4;
  constexpr double b84 = 0.8;
  constexpr double b86 = 0.2;
  constexpr double b31 = 0.6;
  constexpr double b36 = 0.4;
  constexpr double b24 = 0.7;
  constexpr double b27 = 0.3;
  // constexpr double f2min = 0.44;
  // constexpr double f3min = 0.23;
  // constexpr double f5min = 0.09;
  // constexpr double f8min = 0.23;
  // constexpr double f2max = 0.44;
  // constexpr double f3max = 0.23;
  // constexpr double f5max = 0.09;
  // constexpr double f8max = 0.23;
  constexpr double f2min = 0.44;
  constexpr double f3min = 0.23;
  constexpr double f5min = 0.09;
  constexpr double f8min = 0.23;
  constexpr double f2max = 0.46;
  constexpr double f3max = 0.27;
  constexpr double f5max = 0.11;
  constexpr double f8max = 0.27;

  barycentric_affine_approximator::SystemParams system_params{
      N,   F,   v,     w,     b51,   b57,   b84,   b86,   b31,   b36,
      b24, b27, f2min, f3min, f5min, f8min, f2max, f3max, f5max, f8max};

  // Same timing / mode params as affine_solver.cpp
  constexpr double t_max = 1200.0;
  // constexpr double tau_min = 0.06;
  // constexpr double tau_max = 0.12;
  constexpr double tau_min = 10;
  constexpr double tau_max = 50;
  // constexpr double t_max = 1.0;
  constexpr int t_split_count = 240;
  // constexpr int max_switches = 5;

  barycentric_affine_approximator::ApproximationMode approximation_mode
      = barycentric_affine_approximator::ApproximationMode::Lower;
  bool highs_verbose = true;
  return barycentric_affine_approximator::BarycentricAffineApproximator(
      t_max, t_split_count, tau_min, tau_max, system_params, highs_verbose,
      approximation_mode);
}

// Deliberately coarse, so that a complete solution is reachable at N=160.
//
// certificate_tol decides both whether a region counts as violated and whether
// a cached courier may certify it, so it drives the cost twice over. At the
// library default the screen covered a thousand regions out of the million-odd
// this geometry has, and every sweep was one subproblem LP per region, with a
// dozen sweeps per border condition. Loosening it lets one courier stand in for
// many neighbouring regions, which is what made the N=100 run finish in four
// iterations.
//
// The value function here is of order 1e2 to 1e3, so 1e-2 of absolute slack on
// the border condition is roughly 1e-5 relative. Coarse, and deliberately so;
// move it back towards the library default to tighten.
static barycentric_affine_approximator::CourierBorderOptions
coarse_courier_options() {
  barycentric_affine_approximator::CourierBorderOptions options;
  options.certificate_tol = 1e-2;
  options.max_certificate_cache = 512;
  return options;
}

// One band LP per (level, theta) node covers all the time points of that node
// at once (problem lp_horizon), so this run solves about 12 400 LPs per phase
// rather than the 107 000 one-step LPs the greedy march needed. The grid below
// gives 93% of those nodes nine stages, which on this arrangement is roughly
// 59 000 rows by 12 900 columns and 7 s of interior point.
static barycentric_affine_approximator::BandLpOptions band_lp_options() {
  barycentric_affine_approximator::BandLpOptions options;
  // A node that stops converging must not hold a worker for the rest of the
  // run: at this size the honest solve is 7 s, so this is forty times out of
  // the way of one. What comes back at the limit is kept only if it is primal
  // feasible -- the rows are the statement s * F <= 0, so such a point is a
  // coarser bound rather than a wrong one -- and the exact worst residual is
  // re-derived from the formulas afterwards either way. The other solver of
  // the ladder gets its own attempt first.
  options.time_limit = 300.0;
  return options;
}

// Pre-flight for the run below: the whole pipeline up to and including one
// band LP of the widest kind, on the production parameters, without the
// twenty-odd thousand nodes. Everything that fails in the first minute of a
// real run fails here in a few seconds -- the arrangement, the block
// factorisation, the gauge normalisation checked against the kernel of the
// assembled rows, the band assembler, HiGHS, and the residual check that says
// the result is a bound.
//
// Nine stages is what 93% of the nodes of this grid carry: T = 1200 on 240
// points is dt = 5.02, and tau_max / dt = 9.
TEST(common, DISABLED_barycentric_preflight) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  auto approximator = create_barycentric_approximator();
  approximator.setCourierOptions(coarse_courier_options());
  approximator.setBandLpOptions(band_lp_options());

  approximator.getIntersectionPoints();
  // Builds the block system matrices, cross-checks them against the eight-cell
  // assembly, and verifies the gauge pins against the numerical kernel.
  approximator.precomputeMatrices();

  for (int phase = 0; phase < barycentric_affine_approximator::kPhases;
       ++phase) {
    const int num_x = approximator.layouts()[phase].num_x;
    const std::vector<double> z_terminal(static_cast<std::size_t>(num_x), 0.0);
    barycentric_affine_approximator::BandSolveStats stats;
    // Level 0, the terminal family: z(T) = 0 and the third block's plane
    // identically zero.
    const auto z = approximator.solveBandLp(phase, /*switch_cnt=*/0,
                                            /*num_stages=*/9, z_terminal,
                                            &stats);
    ASSERT_EQ(z.size(), 10u);
    EXPECT_LE(stats.worst_residual, stats.residual_limit);
    GTEST_LOG_(INFO) << "phase " << phase << ": " << stats.rows << " x "
                     << stats.cols << ", " << stats.nnz << " nnz, "
                     << stats.solver << " status " << stats.model_status
                     << " in " << stats.seconds << " s, " << stats.iterations
                     << " iterations, integral " << stats.integral
                     << ", worst residual " << stats.worst_residual
                     << " against a limit of " << stats.residual_limit;
  }
}

TEST(common, barycentric_affine_solver) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  auto barycentric_approximator = create_barycentric_approximator();
  barycentric_approximator.setCourierOptions(coarse_courier_options());
  barycentric_approximator.setBandLpOptions(band_lp_options());

  // The default is the production machine's results tree. Overridable so the
  // same test can be run anywhere without editing it; run() throws if the
  // directory does not exist, which is the whole failure on any other host.
  //
  // The thread count is now just the size of the pool over (level, theta,
  // phase) nodes -- it no longer has to be even, because there is no per-phase
  // solver to split it between. Each worker builds its own band LP and its own
  // HiGHS; measured peak resident set for one nine-stage solve, geometry
  // included, is 330 MB, so sixteen of them is an upper bound of about 5 GB.
  const char* out = std::getenv("HCPWA_RESULTS_DIR");
  barycentric_approximator.run(
      out != nullptr
          ? out
          : "/root/gitlab/hybrid-control-pwa/results/barycentric/lower/"
            "first_run/",
      16);
}
