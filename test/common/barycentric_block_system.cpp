#include <gtest/gtest.h>

#include <numeric>
#include <stdexcept>
#include <vector>

#include <Highs.h>

#include <barycentric_affine_approximator.hpp>
#include <courier_border_solver.hpp>

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

barycentric_affine_approximator::SystemParams makeParams() {
  // The parameters of the user_algo pipeline tests: a real arrangement that is
  // still small enough to build here.
  return barycentric_affine_approximator::SystemParams{
      /*N=*/100.0,  /*F=*/15.0,   /*v=*/0.2,    /*w=*/0.5,
      /*b51=*/0.5,  /*b57=*/0.5,  /*b84=*/0.5,  /*b86=*/0.5,
      /*b31=*/0.5,  /*b36=*/0.5,  /*b24=*/0.5,  /*b27=*/0.5,
      /*f2min=*/5.0, /*f3min=*/5.0, /*f5min=*/5.0, /*f8min=*/5.0,
      /*f2max=*/10.0, /*f3max=*/10.0, /*f5max=*/10.0, /*f8max=*/10.0};
}

}  // namespace

TEST(barycentric_block_system, block_matrices_agree_with_the_full_assembly) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  barycentric_affine_approximator::BarycentricAffineApproximator approximator(
      /*t_max=*/300.0, /*t_split_count=*/10, /*tau_min=*/60.0,
      /*tau_max=*/120.0, makeParams(), /*highs_verbose=*/false);

  // The 8D product of the block cells is 90 million vertices on this geometry
  // and nothing in the LP path reads it any more.
  hcpwa::TriangleGeometryOptions options;
  options.build_8d_vertices = false;
  approximator.setGeometryOptions(options);

  approximator.getIntersectionPoints();
  // The cross-check lives inside precomputeMatrices() so that it also runs on
  // production geometry, where no test could afford to rebuild it.
  approximator.precomputeMatrices();
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

  barycentric_affine_approximator::BarycentricAffineApproximator approximator(
      /*t_max=*/300.0, /*t_split_count=*/10, /*tau_min=*/60.0,
      /*tau_max=*/120.0, makeParams(), /*highs_verbose=*/false);
  hcpwa::TriangleGeometryOptions options;
  options.build_8d_vertices = true;
  approximator.setGeometryOptions(options);
  approximator.getIntersectionPoints();
  approximator.precomputeMatrices();
}

// One backward step of the real LP, driven from outside run(). This is the
// first point where the whole chain is exercised on production-shaped
// geometry: the reduced assembler, HiGHS, the per-step RHS update, and the
// unconditional residual check that says the result is a bound.
TEST(barycentric_block_system, one_reduced_step_solves_and_is_a_bound) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  for (const auto mode :
       {barycentric_affine_approximator::ApproximationMode::Upper,
        barycentric_affine_approximator::ApproximationMode::Lower}) {
    barycentric_affine_approximator::BarycentricAffineApproximator approximator(
        /*t_max=*/300.0, /*t_split_count=*/10, /*tau_min=*/60.0,
        /*tau_max=*/120.0, makeParams(), /*highs_verbose=*/false, mode);
    hcpwa::TriangleGeometryOptions options;
    options.build_8d_vertices = false;
    approximator.setGeometryOptions(options);
    approximator.getIntersectionPoints();
    approximator.precomputeMatrices();

    for (int phase = 0; phase < barycentric_affine_approximator::kPhases;
         ++phase) {
      auto [highs, row_lower, base_upper] = approximator.initializeHighs(phase);
      const auto& cols = approximator.reducedLpCols(phase);
      ASSERT_LT(cols.num_cols, 10000);
      ASSERT_LT(approximator.reducedLpRows(phase).num_rows, 100000);

      // Terminal condition of the backward recursion: the value function is
      // zero at the last time layer.
      const std::vector<double> x_next(
          static_cast<std::size_t>(cols.num_x), 0.0);
      const std::vector<double> upper
          = barycentric_affine_approximator::block_reduction::
              updateReducedLpRowUpper(approximator.reducedLpInput(phase),
                                      approximator.reducedLpRows(phase),
                                      base_upper, x_next);
      std::vector<int> row_ids(upper.size());
      std::iota(row_ids.begin(), row_ids.end(), 0);
      ASSERT_EQ(highs->changeRowsBounds(static_cast<int>(row_ids.size()),
                                        row_ids.data(), row_lower.data(),
                                        upper.data()),
                HighsStatus::kOk);
      ASSERT_EQ(highs->run(), HighsStatus::kOk);
      ASSERT_EQ(highs->getModelStatus(), HighsModelStatus::kOptimal)
          << "phase " << phase << " model status "
          << static_cast<int>(highs->getModelStatus());

      const std::vector<double> z(
          highs->getSolution().col_value.begin(),
          highs->getSolution().col_value.begin() + cols.num_x);
      // Throws unless s * F <= 0 everywhere on the product of regions and
      // vertices, which is the whole point of the construction.
      approximator.validateStepResiduals(phase, x_next, z);
    }
  }
}

// The courier border solver on the same production-shaped geometry. Its
// subproblem is reduced over the same blocks: condition (a) splits because
// every target plane lies inside one block and every courier term
// c_s[k] * nu_axis is a single coordinate, so the maximum over the product of
// vertices is the sum of the per-block maxima. The same vertices are checked,
// just counted differently.
TEST(barycentric_block_system, courier_certifies_on_real_geometry) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  barycentric_affine_approximator::BarycentricAffineApproximator approximator(
      /*t_max=*/300.0, /*t_split_count=*/10, /*tau_min=*/60.0,
      /*tau_max=*/120.0, makeParams(), /*highs_verbose=*/false);
  hcpwa::TriangleGeometryOptions options;
  options.build_8d_vertices = false;
  approximator.setGeometryOptions(options);
  approximator.getIntersectionPoints();

  barycentric_affine_approximator::CourierBorderOptions courier_options;
  courier_options.max_iterations = 8;
  barycentric_affine_approximator::CourierBorderSolver solver(courier_options);
  solver.prepare(approximator.phaseGeometries(), approximator.layouts(),
                 approximator.nodeWeights(), 100.0);
  ASSERT_TRUE(solver.prepared());
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

  barycentric_affine_approximator::BarycentricAffineApproximator approximator(
      /*t_max=*/300.0, /*t_split_count=*/10, /*tau_min=*/60.0,
      /*tau_max=*/120.0, makeParams(), /*highs_verbose=*/false);
  hcpwa::TriangleGeometryOptions options;
  options.build_8d_vertices = false;
  approximator.setGeometryOptions(options);
  approximator.getIntersectionPoints();

  barycentric_affine_approximator::CourierBorderOptions courier_options;
  barycentric_affine_approximator::CourierBorderSolver solver(courier_options);
  solver.prepare(approximator.phaseGeometries(), approximator.layouts(),
                 approximator.nodeWeights(), 100.0);

  const int num_x = approximator.layouts()[1].num_x;
  const std::vector<std::vector<double>> candidates
      = {std::vector<double>(static_cast<std::size_t>(num_x), 0.2),
         std::vector<double>(static_cast<std::size_t>(num_x), -0.1)};
  barycentric_affine_approximator::CourierBorderRequest request;
  request.target_phase = 0;
  request.source_phase = 1;
  request.candidates = candidates;

  barycentric_affine_approximator::CourierBorderStats stats;
  const std::vector<double> z = solver.solve(
      request, barycentric_affine_approximator::ApproximationMode::Lower,
      &stats);
  EXPECT_EQ(static_cast<int>(z.size()), approximator.layouts()[0].num_x);
  EXPECT_LE(stats.worst_zeta, courier_options.certificate_tol);
  GTEST_LOG_(INFO) << "converged in " << stats.iterations << " iterations, "
                   << stats.cuts_added << " cuts, "
                   << stats.subproblems_solved << " subproblems";

  // worstCertificateResidual() would re-derive every region's courier from
  // scratch, which is the point of it and also why it is not called here: it
  // has no screening and no early exit, so on this geometry it is an hour of
  // subproblems. solve() already returned only after a complete pass over all
  // 1.58 million regions found nothing, which is the certification itself.
}

// The tie-break is a weak eps * ||z||_1 term on auxiliary columns. It must not
// cost the bound property, and two identical runs must return the same vertex
// of what would otherwise be a degenerate optimal face.
TEST(barycentric_block_system, tie_break_keeps_the_bound_and_repeats) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  auto solveOnce = [](double eps) {
    barycentric_affine_approximator::BarycentricAffineApproximator approximator(
        /*t_max=*/300.0, /*t_split_count=*/10, /*tau_min=*/60.0,
        /*tau_max=*/120.0, makeParams(), /*highs_verbose=*/false);
    approximator.setTieBreakEps(eps);
    approximator.getIntersectionPoints();
    approximator.precomputeMatrices();

    auto [highs, row_lower, base_upper] = approximator.initializeHighs(0);
    const auto& cols = approximator.reducedLpCols(0);
    const std::vector<double> x_next(
        static_cast<std::size_t>(cols.num_x), 0.0);
    const std::vector<double> upper
        = barycentric_affine_approximator::block_reduction::
            updateReducedLpRowUpper(approximator.reducedLpInput(0),
                                    approximator.reducedLpRows(0), base_upper,
                                    x_next);
    std::vector<int> row_ids(upper.size());
    std::iota(row_ids.begin(), row_ids.end(), 0);
    highs->changeRowsBounds(static_cast<int>(row_ids.size()), row_ids.data(),
                            row_lower.data(), upper.data());
    highs->run();
    EXPECT_EQ(highs->getModelStatus(), HighsModelStatus::kOptimal);
    std::vector<double> z(highs->getSolution().col_value.begin(),
                          highs->getSolution().col_value.begin() + cols.num_x);
    approximator.validateStepResiduals(0, x_next, z);
    return z;
  };

  const std::vector<double> first = solveOnce(1e-8);
  const std::vector<double> second = solveOnce(1e-8);
  ASSERT_EQ(first.size(), second.size());
  for (std::size_t k = 0; k < first.size(); ++k) {
    EXPECT_EQ(first[k], second[k]) << "column " << k;
  }
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

  barycentric_affine_approximator::BarycentricAffineApproximator approximator(
      /*t_max=*/300.0, /*t_split_count=*/10, /*tau_min=*/60.0,
      /*tau_max=*/120.0, makeParams(), /*highs_verbose=*/false);
  approximator.getIntersectionPoints();
  approximator.precomputeMatrices();

  for (int phase = 0; phase < barycentric_affine_approximator::kPhases;
       ++phase) {
    const auto& input = approximator.reducedLpInput(phase);
    const auto& geometry
        = approximator.phaseGeometries()[static_cast<std::size_t>(phase)];
    const auto& layout
        = approximator.layouts()[static_cast<std::size_t>(phase)];

    for (int b = 0; b < barycentric_affine_approximator::kBlockCount; ++b) {
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
