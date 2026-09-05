#include <gtest/gtest.h>

#include <chrono>
#include <memory>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <Highs.h>

#include <barycentric_affine_approximator.hpp>
#include <courier_border_solver.hpp>
#include <util/band_lp.hpp>
#include <util/gauge_fix.hpp>

#include "cddwrap/cdd.hpp"
#include "utility.hpp"

// Runs the real arrangement once and lets precomputeMatrices() cross-check the
// per-block CTM and disturbance-box data against the same assemblers run on all
// eight cells at once.
//
// getFIJMinResolution returns rows that depend on the representative point only
// through which of the three branches is the minimum, so agreement between the
// two paths is exact; any difference means they resolved different branches.
namespace {

namespace baa = barycentric_affine_approximator;
namespace br = barycentric_affine_approximator::block_reduction;

baa::SystemParams makeParams() {
  // The parameters of the user_algo pipeline tests: a real arrangement that is
  // still small enough to build here.
  return baa::SystemParams{
      /*N=*/100.0,  /*F=*/15.0,   /*v=*/0.2,    /*w=*/0.5,
      /*b51=*/0.5,  /*b57=*/0.5,  /*b84=*/0.5,  /*b86=*/0.5,
      /*b31=*/0.5,  /*b36=*/0.5,  /*b24=*/0.5,  /*b27=*/0.5,
      /*f2min=*/5.0, /*f3min=*/5.0, /*f5min=*/5.0, /*f8min=*/5.0,
      /*f2max=*/10.0, /*f3max=*/10.0, /*f5max=*/10.0, /*f8max=*/10.0};
}

// Geometry and block data, ready for band solves. The 8D product of the block
// cells is 90 million vertices on this geometry and nothing in the LP path
// reads it any more.
// Behind a unique_ptr: the class owns a mutex through its ValueFunction and is
// neither copyable nor movable.
std::unique_ptr<baa::BarycentricAffineApproximator> makePreparedApproximator(
    baa::ApproximationMode mode = baa::ApproximationMode::Upper) {
  auto approximator = std::make_unique<baa::BarycentricAffineApproximator>(
      /*t_max=*/300.0, /*t_split_count=*/10, /*tau_min=*/60.0,
      /*tau_max=*/120.0, makeParams(), /*highs_verbose=*/false, mode);
  hcpwa::TriangleGeometryOptions options;
  options.build_8d_vertices = false;
  approximator->setGeometryOptions(options);
  approximator->getIntersectionPoints();
  approximator->precomputeMatrices();
  return approximator;
}

}  // namespace

TEST(barycentric_block_system, block_matrices_agree_with_the_full_assembly) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;
  // The cross-check lives inside precomputeMatrices() so that it also runs on
  // production geometry, where no test could afford to rebuild it.
  makePreparedApproximator();
}

// Same check with the 8D product materialised, which additionally lets
// validateBlockDecomposition() compare the concatenated block centroids against
// the product centroid. It costs about six minutes and 90 million vertices on
// this geometry, so it is opt-in: run it with
//   --gtest_also_run_disabled_tests --gtest_filter='*centroid*'
TEST(barycentric_block_system,
     DISABLED_block_centroids_agree_with_the_product_centroid) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  baa::BarycentricAffineApproximator approximator(
      /*t_max=*/300.0, /*t_split_count=*/10, /*tau_min=*/60.0,
      /*tau_max=*/120.0, makeParams(), /*highs_verbose=*/false);
  hcpwa::TriangleGeometryOptions options;
  options.build_8d_vertices = true;
  approximator.setGeometryOptions(options);
  approximator.getIntersectionPoints();
  approximator.precomputeMatrices();
}

// The gauge pins derived from the geometry meet the numerically computed
// kernel transversally -- precomputeMatrices() already throws otherwise, so
// what this test adds is the shape of the answer and that the check is not
// vacuous: dropping one line pin has to be caught.
TEST(barycentric_block_system, gauge_pins_span_the_kernel) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;
  const auto approximator_ptr = makePreparedApproximator();
  auto& approximator = *approximator_ptr;

  for (int phase = 0; phase < baa::kPhases; ++phase) {
    const baa::GaugeFix& gauge
        = approximator.gaugeFixes()[static_cast<std::size_t>(phase)];
    const auto& layout = approximator.layouts()[static_cast<std::size_t>(phase)];
    ASSERT_EQ(gauge.pairs.size(), 2u) << "phase " << phase;
    EXPECT_EQ(gauge.constant_pins.size(), 2u) << "phase " << phase;
    EXPECT_GE(gauge.line_pins.size(), 4u) << "phase " << phase;
    EXPECT_EQ(static_cast<int>(gauge.group_c_columns.size()),
              layout.eta_s[static_cast<std::size_t>(gauge.group_c_layer)])
        << "phase " << phase;
    for (const baa::GaugeFix::Pair& pair : gauge.pairs) {
      ASSERT_GE(pair.xi.size(), 2u);
      EXPECT_NEAR(pair.xi.front(), 0.0, 1e-9);
      EXPECT_NEAR(pair.xi.back(), 100.0, 1e-9);
    }

    const baa::GaugeFixReport report = baa::verifyGaugeFix(
        gauge, approximator.reducedLpInput(phase), layout);
    EXPECT_EQ(report.line_rank, report.kernel_dim) << "phase " << phase;
    EXPECT_EQ(report.full_rank, report.kernel_dim + 2) << "phase " << phase;
    int predicted = 0;
    for (const baa::GaugeFix::Pair& pair : gauge.pairs) {
      predicted += static_cast<int>(pair.xi.size());
    }
    EXPECT_EQ(report.kernel_dim, predicted) << "phase " << phase;

    // The two levels pin different sets, and the constant pins never appear
    // at r = 0.
    const std::vector<int> pins0 = gauge.pinsForLevel(0);
    const std::vector<int> pins1 = gauge.pinsForLevel(1);
    EXPECT_EQ(pins0.size(),
              gauge.line_pins.size() + gauge.group_c_columns.size());
    EXPECT_EQ(pins1.size(),
              gauge.line_pins.size() + gauge.constant_pins.size());
    for (const int col : gauge.constant_pins) {
      EXPECT_EQ(std::find(pins0.begin(), pins0.end(), col), pins0.end());
    }

    // Not vacuous: without one line pin a kernel direction survives.
    baa::GaugeFix crippled = gauge;
    crippled.line_pins.pop_back();
    EXPECT_THROW(baa::verifyGaugeFix(crippled, approximator.reducedLpInput(phase),
                                     layout),
                 std::runtime_error)
        << "phase " << phase;
  }
}

// The band LP on production-shaped geometry, driven from outside run(): the
// assembler, HiGHS, the gauge pins and the unconditional residual check that
// says the result is a bound on every stage. Level 0: z_L = 0 and the block-C
// plane pinned to zero.
TEST(barycentric_block_system, band_lp_solves_and_is_a_bound) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  for (const auto mode :
       {baa::ApproximationMode::Upper, baa::ApproximationMode::Lower}) {
    const auto approximator_ptr = makePreparedApproximator(mode);
  auto& approximator = *approximator_ptr;
    for (int phase = 0; phase < baa::kPhases; ++phase) {
      const int num_x
          = approximator.layouts()[static_cast<std::size_t>(phase)].num_x;
      const std::vector<double> z_terminal(static_cast<std::size_t>(num_x),
                                           0.0);
      for (const int stages : {1, 3}) {
        baa::BandSolveStats stats;
        // Throws unless s * F <= 0 everywhere on every stage.
        const std::vector<std::vector<double>> z = approximator.solveBandLp(
            phase, /*switch_cnt=*/0, stages, z_terminal, &stats);
        ASSERT_EQ(static_cast<int>(z.size()), stages + 1);
        EXPECT_EQ(stats.rows, stages * approximator.reducedLpRows(phase).num_rows);
        EXPECT_EQ(stats.cols, stages * approximator.reducedLpCols(phase).num_cols);
        EXPECT_EQ(stats.model_status, static_cast<int>(HighsModelStatus::kOptimal))
            << "phase " << phase << " stages " << stages;
        EXPECT_FALSE(stats.accepted_non_optimal);
        EXPECT_TRUE(std::isfinite(stats.integral));
        // The integral of a lower bound of a nonnegative cost is nonnegative
        // and of an upper one too; both must be finite and of the right sign
        // relative to zero terminal data: V(t) >= 0 for t < T in both modes,
        // because the running cost g is nonnegative.
        EXPECT_GE(stats.integral, -1e-6) << "phase " << phase;
        // The block-C plane is identically zero at level 0.
        const baa::GaugeFix& gauge
            = approximator.gaugeFixes()[static_cast<std::size_t>(phase)];
        for (int k = 0; k < stages; ++k) {
          for (const int col : gauge.group_c_columns) {
            EXPECT_EQ(z[static_cast<std::size_t>(k)][static_cast<std::size_t>(col)],
                      0.0)
                << "phase " << phase << " stage " << k << " column " << col;
          }
          for (const int col : gauge.line_pins) {
            EXPECT_EQ(z[static_cast<std::size_t>(k)][static_cast<std::size_t>(col)],
                      0.0)
                << "phase " << phase << " stage " << k << " column " << col;
          }
        }
        GTEST_LOG_(INFO) << "mode " << (mode == baa::ApproximationMode::Upper
                                            ? "upper" : "lower")
                         << " phase " << phase << " stages " << stages
                         << ": " << stats.rows << "x" << stats.cols << ", "
                         << stats.nnz << " nnz, " << stats.solver << " in "
                         << stats.seconds << " s, " << stats.iterations
                         << " iterations, integral " << stats.integral
                         << ", worst residual " << stats.worst_residual;
      }
    }
  }
}

// The remark on the greedy scheme, on real geometry: the greedy march (the
// band LP with one stage, repeated backwards) is feasible for the band LP and
// therefore no better by the integral criterion.
TEST(barycentric_block_system, band_beats_greedy_on_real_geometry) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  for (const auto mode :
       {baa::ApproximationMode::Upper, baa::ApproximationMode::Lower}) {
    const auto approximator_ptr = makePreparedApproximator(mode);
  auto& approximator = *approximator_ptr;
    const int phase = 0;
    const int stages = 3;
    const int num_x
        = approximator.layouts()[static_cast<std::size_t>(phase)].num_x;
    const std::vector<double> z_terminal(static_cast<std::size_t>(num_x), 0.0);

    std::vector<std::vector<double>> greedy(static_cast<std::size_t>(stages) + 1);
    greedy.back() = z_terminal;
    for (int k = stages - 1; k >= 0; --k) {
      greedy[static_cast<std::size_t>(k)] = approximator.solveBandLp(
          phase, 0, 1, greedy[static_cast<std::size_t>(k + 1)])[0];
    }
    baa::BandSolveStats stats;
    const std::vector<std::vector<double>> band
        = approximator.solveBandLp(phase, 0, stages, z_terminal, &stats);

    const Eigen::VectorXd& q
        = approximator.nodeWeights()[static_cast<std::size_t>(phase)];
    const double integral_greedy
        = br::bandIntegral(q, approximator.tDelta(), greedy);
    const double integral_band = br::bandIntegral(q, approximator.tDelta(), band);
    EXPECT_NEAR(integral_band, stats.integral,
                1e-9 * std::max(1.0, std::abs(integral_band)));
    // The greedy point is certified on every stage ...
    EXPECT_GE(approximator.validateBandResiduals(phase, greedy), -1e300);
    // ... and hence no better.
    const double tol = 1e-7 * std::max(1.0, std::abs(integral_greedy));
    if (mode == baa::ApproximationMode::Lower) {
      EXPECT_GE(integral_band, integral_greedy - tol);
    } else {
      EXPECT_LE(integral_band, integral_greedy + tol);
    }
    GTEST_LOG_(INFO) << (mode == baa::ApproximationMode::Upper ? "upper"
                                                                 : "lower")
                     << ": greedy integral " << integral_greedy
                     << ", band integral " << integral_band;
  }
}

// Two identical band solves return the same node values, bit for bit.
TEST(barycentric_block_system, band_lp_is_reproducible) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  const auto approximator_ptr = makePreparedApproximator();
  auto& approximator = *approximator_ptr;
  const int num_x = approximator.layouts()[0].num_x;
  const std::vector<double> z_terminal(static_cast<std::size_t>(num_x), 0.0);
  const auto first = approximator.solveBandLp(0, 0, 2, z_terminal);
  const auto second = approximator.solveBandLp(0, 0, 2, z_terminal);
  ASSERT_EQ(first.size(), second.size());
  for (std::size_t k = 0; k < first.size(); ++k) {
    ASSERT_EQ(first[k].size(), second[k].size());
    for (std::size_t i = 0; i < first[k].size(); ++i) {
      EXPECT_EQ(first[k][i], second[k][i]) << "stage " << k << " column " << i;
    }
  }
}

// A terminal condition that violates the gauge pins is refused: it would be
// infeasible by the band LP's own column bounds, and it means the border
// problem was normalised differently.
TEST(barycentric_block_system, band_lp_refuses_an_unpinned_terminal_value) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  const auto approximator_ptr = makePreparedApproximator();
  auto& approximator = *approximator_ptr;
  const int num_x = approximator.layouts()[0].num_x;
  std::vector<double> z_terminal(static_cast<std::size_t>(num_x), 0.0);
  const int pinned = approximator.gaugeFixes()[0].pinsForLevel(1).front();
  z_terminal[static_cast<std::size_t>(pinned)] = 1e-3;
  EXPECT_THROW(approximator.solveBandLp(0, 1, 2, z_terminal),
               std::runtime_error);
}

// Timing of the solver ladder's members on one band LP, for choosing the
// default. Informational: run with --gtest_also_run_disabled_tests.
TEST(barycentric_block_system, DISABLED_band_lp_solver_timing) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  const auto approximator_ptr = makePreparedApproximator(baa::ApproximationMode::Lower);
  auto& approximator = *approximator_ptr;
  const int num_x = approximator.layouts()[0].num_x;
  const std::vector<double> z_terminal(static_cast<std::size_t>(num_x), 0.0);
  struct Config {
    const char* name;
    const char* solver;
    bool crossover;
  };
  for (const Config config : {Config{"ipm+crossover", "ipm", true},
                              Config{"ipm", "ipm", false},
                              Config{"simplex", "simplex", false}}) {
    baa::BandLpOptions options;
    options.solver = config.solver;
    options.run_crossover = config.crossover;
    approximator.setBandLpOptions(options);
    for (const int stages : {1, 3, 9}) {
      baa::BandSolveStats stats;
      approximator.solveBandLp(0, 0, stages, z_terminal, &stats);
      GTEST_LOG_(INFO) << config.name << " stages " << stages << ": "
                       << stats.rows << "x" << stats.cols << ", "
                       << stats.seconds << " s, " << stats.iterations
                       << " iterations, status " << stats.model_status
                       << (stats.accepted_non_optimal ? " (non-optimal)" : "")
                       << ", integral " << stats.integral
                       << ", worst residual " << stats.worst_residual;
    }
  }
}

// Runs the border problem to convergence on production-shaped geometry: 1.58
// million target regions, certified region by region.
//
// This is the gate that says the courier method is usable at this scale. It
// converges in four Benders iterations, and its last sweep -- the complete pass
// over every region that certification requires -- takes a third of a second,
// because one cached courier covers all of them.
TEST(barycentric_block_system, courier_converges_on_real_geometry) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  const auto approximator_ptr = makePreparedApproximator();
  auto& approximator = *approximator_ptr;

  baa::CourierBorderOptions courier_options;
  baa::CourierBorderSolver solver(courier_options);
  solver.prepare(approximator.phaseGeometries(), approximator.layouts(),
                 approximator.nodeWeights(), 100.0);

  const int num_x = approximator.layouts()[1].num_x;
  const std::vector<std::vector<double>> candidates
      = {std::vector<double>(static_cast<std::size_t>(num_x), 0.2),
         std::vector<double>(static_cast<std::size_t>(num_x), -0.1)};
  // Level 1 pins of the target phase, as run() would hand them over.
  const std::vector<int> pins = approximator.gaugeFixes()[0].pinsForLevel(1);
  baa::CourierBorderRequest request;
  request.target_phase = 0;
  request.source_phase = 1;
  request.candidates = candidates;
  request.pinned_columns = pins;
  request.switch_cnt = 1;

  baa::CourierBorderStats stats;
  const std::vector<double> z
      = solver.solve(request, baa::ApproximationMode::Lower, &stats);
  EXPECT_EQ(static_cast<int>(z.size()), approximator.layouts()[0].num_x);
  for (const int col : pins) {
    EXPECT_EQ(z[static_cast<std::size_t>(col)], 0.0) << "column " << col;
  }
  // Deliberately not EXPECT_LE(stats.worst_zeta, certificate_tol): that
  // cannot fail. worst is only raised past the certificate_tol
  // early-continue, and any region that raises it also pushes a
  // violation, which prevents the full-pass return. The certification is
  // the return itself -- solve() throws unless one complete pass over
  // every region found nothing. What is worth pinning is that the run did
  // real work rather than certifying vacuously.
  EXPECT_GT(stats.subproblems_solved, 0);
  EXPECT_GT(stats.cuts_added, 0);
  EXPECT_GE(stats.iterations, 1);
  GTEST_LOG_(INFO) << "converged in " << stats.iterations << " iterations, "
                   << stats.cuts_added << " cuts, "
                   << stats.subproblems_solved << " subproblems";

  // And the node values it returned are a usable terminal condition of the
  // level-1 band LP: same pins, so the band accepts them.
  baa::BandSolveStats band_stats;
  const auto band = approximator.solveBandLp(0, 1, 2, z, &band_stats);
  EXPECT_EQ(static_cast<int>(band.size()), 3);
  GTEST_LOG_(INFO) << "level-1 band from the courier's z: integral "
                   << band_stats.integral << ", worst residual "
                   << band_stats.worst_residual;

  // worstCertificateResidual() would re-derive every region's courier from
  // scratch, which is the point of it and also why it is not called here: it
  // has no screening and no early exit, so on this geometry it is an hour of
  // subproblems. solve() already returned only after a complete pass over all
  // 1.58 million regions found nothing, which is the certification itself.
}

// Two structural invariants of the assembled block LP that nothing else pins.
//
// Psi_b 1 = 0 says a uniform shift of the barycentric node values leaves the
// gradient alone -- it is what makes the gauge freedom a nullspace of the
// dynamics and, through that, what the boundedness argument for the objective
// weights rests on. And Psi_b must touch only the x columns of its own block's
// projection planes, which is the disjointness the three block rows are summed
// under.
TEST(barycentric_block_system, block_psi_rows_are_gauge_free_and_local) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  const auto approximator_ptr = makePreparedApproximator();
  auto& approximator = *approximator_ptr;

  for (int phase = 0; phase < baa::kPhases; ++phase) {
    const auto& input = approximator.reducedLpInput(phase);
    const auto& geometry
        = approximator.phaseGeometries()[static_cast<std::size_t>(phase)];
    const auto& layout
        = approximator.layouts()[static_cast<std::size_t>(phase)];

    for (int b = 0; b < baa::kBlockCount; ++b) {
      const auto& block_geometry
          = geometry.blocks[static_cast<std::size_t>(b)];
      // The x columns this block is allowed to touch.
      std::vector<bool> owned(static_cast<std::size_t>(layout.num_x), false);
      for (int l = 0; l < block_geometry.layer_count; ++l) {
        const int s = block_geometry.layer_ids[static_cast<std::size_t>(l)];
        for (int k = 0; k < layout.eta_s[static_cast<std::size_t>(s)]; ++k) {
          owned[static_cast<std::size_t>(layout.idxX(s, k))] = true;
        }
      }

      const auto& regions = input.blocks[static_cast<std::size_t>(b)].regions;
      for (std::size_t j = 0; j < regions.size(); ++j) {
        for (const auto& row : regions[j].psi_rows) {
          double mass = 0.0;
          for (std::size_t i = 0; i < row.cols.size(); ++i) {
            EXPECT_TRUE(owned[static_cast<std::size_t>(row.cols[i])])
                << "phase " << phase << " block " << b << " cell " << j
                << " touches column " << row.cols[i]
                << " outside its own projection planes";
            mass += row.vals[i];
          }
          EXPECT_NEAR(mass, 0.0, 1e-9)
              << "phase " << phase << " block " << b << " cell " << j;
        }
      }
    }
  }
}

