#include <gtest/gtest.h>

#include <numeric>
#include <vector>

#include <Highs.h>

#include <barycentric_affine_approximator.hpp>

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
