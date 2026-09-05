#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <random>
#include <stdexcept>
#include <vector>

#include <Eigen/Core>
#include <Highs.h>
#include <lp_data/HighsModelUtils.h>
#include <string>

#include <util/band_lp.hpp>
#include <util/block_reduction_lp.hpp>

#include "barycentric_fixture.hpp"

// Tier 1 of the band LP: the L-stage LP against the one-step LP it generalises,
// and the two statements of the paper it exists for -- every feasible point is
// a certified bound on every stage (theorem on the optimality of the band LP,
// part 1), and the band is never worse than the greedy march by the integral
// criterion (remark on the greedy scheme).
namespace {

using namespace barycentric_test_fixture;

// Quadrature weights for the fixture: positive, like the integrals of hat
// functions they stand for.
Eigen::VectorXd makeWeights() {
  std::mt19937 rng(kSeed + 202);
  std::uniform_real_distribution<double> pos(0.5, 2.0);
  Eigen::VectorXd q(totalNumX());
  for (int i = 0; i < q.size(); ++i) {
    q(i) = pos(rng);
  }
  return q;
}

br::BandLpInput makeBandInput(const br::ReducedLpInput& stage, int num_stages,
                              const std::vector<double>& z_terminal,
                              std::vector<int> pins = {}) {
  br::BandLpInput input;
  input.stage = &stage;
  input.num_stages = num_stages;
  input.z_terminal = z_terminal;
  input.node_weights = makeWeights();
  input.pinned_columns = std::move(pins);
  return input;
}

// Solves the band LP and unpacks z_0 .. z_{L-1}, z_L.
//
// The x columns get the same +-kBox box the one-step tests use. On real
// geometry the band LP is bounded by the theorem (every feasible point is a
// certified bound, so its integral is bounded by the value function's), but
// the fixture's m, rho, g and q are random and certify nothing: without a box
// it is unbounded along a transfer of constants between blocks through the
// coupling rows, the very direction the paper's gauge pins remove on the
// real problem.
std::vector<std::vector<double>> solveBand(const br::BandLpMatrices& lp,
                                           const br::BandLpInput& input) {
  const Eigen::RowVectorXd cost = lp.cost / lp.cost_scale;
  std::vector<double> col_lower = lp.col_lower;
  std::vector<double> col_upper = lp.col_upper;
  for (int k = 0; k < input.num_stages; ++k) {
    for (int i = 0; i < lp.cols.stage.num_x; ++i) {
      const std::size_t c = static_cast<std::size_t>(lp.cols.idxX(k, i));
      if (col_lower[c] == 0.0 && col_upper[c] == 0.0) {
        continue;  // pinned
      }
      col_lower[c] = -kBox;
      col_upper[c] = kBox;
    }
  }
  const LpResult result
      = solveLp(lp.starts, lp.col_index, lp.value, lp.row_lower, lp.row_upper,
                col_lower, col_upper, cost);
  if (!result.optimal) {
    throw std::runtime_error(
        "band LP not optimal, model status "
        + std::to_string(static_cast<int>(result.status)) + " ("
        + utilModelStatusToString(result.status) + ")");
  }
  std::vector<std::vector<double>> z;
  for (int k = 0; k < input.num_stages; ++k) {
    std::vector<double> z_k(static_cast<std::size_t>(lp.cols.stage.num_x));
    for (int i = 0; i < lp.cols.stage.num_x; ++i) {
      z_k[static_cast<std::size_t>(i)]
          = result.col_value[static_cast<std::size_t>(lp.cols.idxX(k, i))];
    }
    z.push_back(std::move(z_k));
  }
  z.push_back(input.z_terminal);
  return z;
}

// The greedy march of the remark: L one-step solves backwards in time, each
// one the band LP with a single stage and the previous answer as z_L.
std::vector<std::vector<double>> solveGreedy(const br::ReducedLpInput& stage,
                                             int num_stages,
                                             const std::vector<double>& z_L) {
  std::vector<std::vector<double>> z(static_cast<std::size_t>(num_stages + 1));
  z.back() = z_L;
  for (int k = num_stages - 1; k >= 0; --k) {
    const br::BandLpInput input
        = makeBandInput(stage, 1, z[static_cast<std::size_t>(k + 1)]);
    const br::BandLpMatrices lp = br::assembleBandLp(input);
    z[static_cast<std::size_t>(k)] = solveBand(lp, input)[0];
  }
  return z;
}

double worstStage(const std::vector<br::WorstResidual>& residuals) {
  double worst = -1e300;
  for (const auto& r : residuals) {
    worst = std::max(worst, r.worst());
  }
  return worst;
}

}  // namespace

TEST(barycentric_band_lp, one_stage_matches_the_one_step_lp) {
  const std::vector<double> x_next = makeXNext();
  for (const double s : {1.0, -1.0}) {
    br::ReducedLpInput stage
        = makeInput(s, br::ObjectiveWeights::ProductCount);
    br::clampBlockRho(stage);

    const br::ReducedLpMatrices one_step = br::assembleReducedLp(stage);
    const std::vector<double> one_step_upper = br::updateReducedLpRowUpper(
        stage, one_step.rows, one_step.row_upper, x_next);

    const br::BandLpInput input = makeBandInput(stage, 1, x_next);
    const br::BandLpMatrices band = br::assembleBandLp(input);

    ASSERT_EQ(band.rows.num_rows, one_step.rows.num_rows) << "s=" << s;
    ASSERT_EQ(band.cols.num_cols, one_step.cols.num_cols) << "s=" << s;
    ASSERT_EQ(band.starts, one_step.starts) << "s=" << s;
    ASSERT_EQ(band.col_index, one_step.col_index) << "s=" << s;
    ASSERT_EQ(band.value.size(), one_step.value.size());
    for (std::size_t i = 0; i < band.value.size(); ++i) {
      EXPECT_NEAR(band.value[i], one_step.value[i], 1e-12)
          << "s=" << s << " entry " << i;
    }
    for (int r = 0; r < band.rows.num_rows; ++r) {
      EXPECT_EQ(band.row_lower[static_cast<std::size_t>(r)],
                one_step.row_lower[static_cast<std::size_t>(r)])
          << "s=" << s << " row " << r;
      EXPECT_NEAR(band.row_upper[static_cast<std::size_t>(r)],
                  one_step_upper[static_cast<std::size_t>(r)], 1e-12)
          << "s=" << s << " row " << r;
    }
    for (int c = 0; c < band.cols.num_cols; ++c) {
      EXPECT_EQ(band.col_lower[static_cast<std::size_t>(c)],
                one_step.col_lower[static_cast<std::size_t>(c)])
          << "s=" << s << " column " << c;
      EXPECT_EQ(band.col_upper[static_cast<std::size_t>(c)],
                one_step.col_upper[static_cast<std::size_t>(c)])
          << "s=" << s << " column " << c;
    }
  }
}

TEST(barycentric_band_lp, layout_sizes_scale_with_stages) {
  br::ReducedLpInput stage
      = makeInput(1.0, br::ObjectiveWeights::ProductCount);
  br::clampBlockRho(stage);
  const int L = 3;
  const br::BandLpInput input = makeBandInput(stage, L, makeXNext());
  const br::BandLpMatrices band = br::assembleBandLp(input);

  // 90 rows / 43 columns per stage: the hand-checked one-step sizes.
  EXPECT_EQ(band.rows.stage_rows, 90);
  EXPECT_EQ(band.cols.stage_cols, 43);
  EXPECT_EQ(band.rows.num_rows, 3 * 90);
  EXPECT_EQ(band.cols.num_cols, 3 * 43);
  EXPECT_EQ(static_cast<int>(band.starts.size()), band.rows.num_rows + 1);

  for (int k = 0; k < L; ++k) {
    EXPECT_EQ(band.cols.idxX(k, 5), k * 43 + 5);
    EXPECT_EQ(band.cols.idxMuL(k, 2), k * 43 + band.cols.stage.idxMuL(2));
    EXPECT_EQ(band.rows.rowSumL(k), k * 90 + 89);
    EXPECT_EQ(band.rows.rowSumR(k), k * 90 + 88);
    EXPECT_EQ(band.rows.rowLeft(k, 0, 0, 0), k * 90);
  }
  EXPECT_THROW(band.cols.idxX(L, 0), std::invalid_argument);
  EXPECT_THROW(band.rows.rowSumL(-1), std::invalid_argument);

  // Every row of a non-terminal stage that involves z_{k+1} touches columns of
  // the next stage; the terminal stage touches only its own.
  for (int k = 0; k < L; ++k) {
    int max_col = -1;
    for (int r = k * 90; r < (k + 1) * 90; ++r) {
      for (int e = band.starts[static_cast<std::size_t>(r)];
           e < band.starts[static_cast<std::size_t>(r) + 1]; ++e) {
        max_col = std::max(max_col, band.col_index[static_cast<std::size_t>(e)]);
      }
    }
    if (k + 1 < L) {
      EXPECT_GE(max_col, (k + 1) * 43) << "stage " << k;
      EXPECT_LT(max_col, (k + 2) * 43) << "stage " << k;
    } else {
      EXPECT_LT(max_col, (k + 1) * 43) << "stage " << k;
    }
  }
}

TEST(barycentric_band_lp, objective_is_the_trapezoid_rule) {
  const int L = 3;
  for (const double s : {1.0, -1.0}) {
    br::ReducedLpInput stage
        = makeInput(s, br::ObjectiveWeights::ProductCount);
    br::clampBlockRho(stage);
    const br::BandLpInput input = makeBandInput(stage, L, makeXNext());
    const br::BandLpMatrices band = br::assembleBandLp(input);

    double max_abs = 0.0;
    for (int k = 0; k < L; ++k) {
      const double omega = k == 0 ? 0.5 : 1.0;
      for (int i = 0; i < stage.num_x; ++i) {
        EXPECT_DOUBLE_EQ(band.cost(band.cols.idxX(k, i)),
                         s * omega * input.node_weights(i))
            << "s=" << s << " stage " << k << " node " << i;
        max_abs = std::max(max_abs, std::abs(band.cost(band.cols.idxX(k, i))));
      }
      for (int c = k * band.cols.stage_cols + stage.num_x;
           c < (k + 1) * band.cols.stage_cols; ++c) {
        EXPECT_EQ(band.cost(c), 0.0) << "s=" << s << " column " << c;
      }
    }
    EXPECT_DOUBLE_EQ(band.cost_scale, max_abs);
    EXPECT_DOUBLE_EQ(band.cost_scale, input.node_weights.maxCoeff());
  }

  // bandIntegral is the trapezoid rule with both ends halved.
  Eigen::VectorXd q(2);
  q << 1.0, 2.0;
  const std::vector<std::vector<double>> z = {{1.0, 1.0}, {2.0, 0.0}, {0.0, 1.0}};
  // q^T z_k = 3, 2, 2 ; weights 1/2, 1, 1/2 ; dt = 4.
  EXPECT_DOUBLE_EQ(br::bandIntegral(q, 4.0, z), 4.0 * (1.5 + 2.0 + 1.0));
  EXPECT_DOUBLE_EQ(br::bandIntegral(q, 4.0, {{1.0, 1.0}}), 0.0);
}

TEST(barycentric_band_lp, solution_is_a_bound_on_every_stage) {
  const int L = 3;
  const std::vector<int> pins = {0, 7};
  for (const double s : {1.0, -1.0}) {
    br::ReducedLpInput stage
        = makeInput(s, br::ObjectiveWeights::ProductCount);
    br::clampBlockRho(stage);
    std::vector<double> z_L = makeXNext();
    for (const int c : pins) {
      z_L[static_cast<std::size_t>(c)] = 0.0;
    }
    const br::BandLpInput input = makeBandInput(stage, L, z_L, pins);
    const br::BandLpMatrices band = br::assembleBandLp(input);
    const std::vector<std::vector<double>> z = solveBand(band, input);
    ASSERT_EQ(static_cast<int>(z.size()), L + 1);

    // s * F <= 0 at both ends of every segment: the certificate of part 1 of
    // the theorem, recomputed from the formulas rather than read off HiGHS.
    const std::vector<br::WorstResidual> residuals
        = br::bandStageResiduals(stage, z);
    ASSERT_EQ(static_cast<int>(residuals.size()), L);
    for (int k = 0; k < L; ++k) {
      EXPECT_LE(residuals[static_cast<std::size_t>(k)].worst(), 1e-7)
          << "s=" << s << " stage " << k;
    }
    for (int k = 0; k < L; ++k) {
      for (const int c : pins) {
        EXPECT_EQ(z[static_cast<std::size_t>(k)][static_cast<std::size_t>(c)],
                  0.0)
            << "s=" << s << " stage " << k << " pin " << c;
      }
    }
  }
}

TEST(barycentric_band_lp, band_is_no_worse_than_greedy) {
  const int L = 3;
  for (const double s : {1.0, -1.0}) {
    br::ReducedLpInput stage
        = makeInput(s, br::ObjectiveWeights::ProductCount);
    br::clampBlockRho(stage);
    const std::vector<double> z_L = makeXNext();
    const Eigen::VectorXd q = makeWeights();

    const std::vector<std::vector<double>> greedy = solveGreedy(stage, L, z_L);
    const br::BandLpInput input = makeBandInput(stage, L, z_L);
    const br::BandLpMatrices lp = br::assembleBandLp(input);
    const std::vector<std::vector<double>> band = solveBand(lp, input);

    // The greedy point is feasible for the band LP ...
    EXPECT_LE(worstStage(br::bandStageResiduals(stage, greedy)), 1e-7)
        << "s=" << s;
    EXPECT_LE(worstStage(br::bandStageResiduals(stage, band)), 1e-7)
        << "s=" << s;

    // ... and therefore no better by the integral criterion: the lower bound
    // maximises the integral, the upper bound minimises it.
    const double integral_greedy = br::bandIntegral(q, stage.t_delta, greedy);
    const double integral_band = br::bandIntegral(q, stage.t_delta, band);
    const double tol = 1e-9 * std::max(1.0, std::abs(integral_greedy));
    if (s < 0.0) {
      EXPECT_GE(integral_band, integral_greedy - tol);
    } else {
      EXPECT_LE(integral_band, integral_greedy + tol);
    }
    // And the LP objective is that integral, up to the constant z_L term and
    // the dt factor.
    double objective = 0.0;
    for (int k = 0; k < L; ++k) {
      for (int i = 0; i < stage.num_x; ++i) {
        objective += lp.cost(lp.cols.idxX(k, i))
                     * band[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)];
      }
    }
    double terminal = 0.0;
    for (int i = 0; i < stage.num_x; ++i) {
      terminal += 0.5 * q(i) * z_L[static_cast<std::size_t>(i)];
    }
    EXPECT_NEAR(stage.t_delta * (s * objective + terminal), integral_band,
                1e-9 * std::max(1.0, std::abs(integral_band)))
        << "s=" << s;
  }
}

TEST(barycentric_band_lp, rejects_malformed_input) {
  br::ReducedLpInput stage
      = makeInput(1.0, br::ObjectiveWeights::ProductCount);
  br::clampBlockRho(stage);
  const std::vector<double> z_L = makeXNext();

  br::BandLpInput tie_break = makeBandInput(stage, 2, z_L);
  br::ReducedLpInput with_tie_break = stage;
  with_tie_break.tie_break_eps = 1e-8;
  tie_break.stage = &with_tie_break;
  EXPECT_THROW(br::assembleBandLp(tie_break), std::invalid_argument);

  EXPECT_THROW(br::assembleBandLp(makeBandInput(stage, 0, z_L)),
               std::invalid_argument);
  EXPECT_THROW(br::assembleBandLp(makeBandInput(stage, 2, z_L, {stage.num_x})),
               std::invalid_argument);
  EXPECT_THROW(br::assembleBandLp(makeBandInput(stage, 2, z_L, {1, 1})),
               std::invalid_argument);

  br::BandLpInput short_terminal = makeBandInput(stage, 2, z_L);
  short_terminal.z_terminal.pop_back();
  EXPECT_THROW(br::assembleBandLp(short_terminal), std::invalid_argument);

  // A radius HiGHS would drop is refused by the band assembler too.
  br::ReducedLpInput tiny = stage;
  tiny.blocks[0].regions[0].vertices[0].rho(0) = 1e-12;
  EXPECT_THROW(br::assembleBandLp(makeBandInput(tiny, 2, z_L)),
               std::runtime_error);
}
