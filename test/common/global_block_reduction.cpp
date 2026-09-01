// Tier-1 verification of the block reduction for the GLOBAL affine
// approximator (docs/barycentric_block_reduction_context.md part VII, tier 1,
// applied to the second code path).
//
// Same discipline as test/common/barycentric_block_reduction.cpp and for the
// same reason: the product LP has one row per element of R_A x R_B x R_C, about
// 2e8 rows on the production geometry, so it cannot be built there and cannot
// serve as a runtime baseline. Equivalence is established once, on synthetic
// data of the right shape, against a product LP assembled here in the test.
//
// The global LP is a different animal from the barycentric one. It has a fixed
// 17 columns for the whole phase, [V(8) | v(1) | s(8)], one feasibility row per
// (area, vertex), and a single global modulus lift s. There is no per-region y,
// so the reduction needs no analogue of Lemma 5; what it does need is that the
// rows separate, which follows from A being block diagonal and Q_c, Q_r
// diagonal.
//
// As in the barycentric test, equivalence holds only under the product weights
// R/R_block. Production normalises to 1, which is a different objective with a
// different optimum, so that setting is checked for feasibility in the product
// LP instead -- the property that keeps the certificate valid.

#include <gtest/gtest.h>

#include <Highs.h>

#include <Eigen/Core>
#include <array>
#include <cmath>
#include <format>
#include <limits>
#include <random>
#include <vector>

#include "test_utils.hpp"
#include "util/global_block_reduction_lp.hpp"

namespace {

namespace gblp = hcpwa::util::global_block_lp;

constexpr int kBlocks = gblp::kBlockCount;
constexpr int kDim = gblp::kSpaceDim;
constexpr double kDt = 0.25;
constexpr double kSPin = 1e-6;
// Synthetic p, r and kappa do not satisfy the structural properties the real
// feasibility argument rests on, so the synthetic problem needs a box to stay
// bounded.
constexpr double kBox = 10.0;

// A coordinate partition of the same shape as the real one: 3 + 3 + 2. The
// coordinates are deliberately NOT contiguous per block, matching phase 0
// ({0,2,5}, {1,3,6}, {4,7}), so that any assumption of contiguity in the
// assembler shows up here.
const std::array<std::array<int, 3>, kBlocks> kCoords
    = {std::array<int, 3>{0, 2, 5}, std::array<int, 3>{1, 3, 6},
       std::array<int, 3>{4, 7, -1}};
constexpr std::array<int, kBlocks> kCoordCount = {3, 3, 2};
constexpr std::array<int, kBlocks> kRowsPerBlock = {4, 3, 3};

struct Instance {
  std::array<std::vector<gblp::GlobalBlockRow>, kBlocks> rows;
  std::vector<double> v_prev;  // [V_prev(8); v_prev(1)]
  int product_rows = 0;
};

Instance makeInstance(std::uint64_t seed) {
  std::mt19937_64 rng(seed);
  std::normal_distribution<double> normal(0.0, 1.0);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);

  Instance instance;
  instance.product_rows = 1;
  for (int b = 0; b < kBlocks; ++b) {
    for (int k = 0; k < kRowsPerBlock[b]; ++k) {
      gblp::GlobalBlockRow row;
      row.nu = Eigen::VectorXd(kCoordCount[b]);
      row.p = Eigen::VectorXd(kCoordCount[b]);
      row.r = Eigen::VectorXd(kCoordCount[b]);
      for (int d = 0; d < kCoordCount[b]; ++d) {
        row.nu(d) = normal(rng);
        row.p(d) = normal(rng);
        row.r(d) = uniform(rng);  // radius >= 0
      }
      row.kappa = normal(rng);
      instance.rows[b].push_back(std::move(row));
    }
    instance.product_rows *= kRowsPerBlock[b];
  }

  instance.v_prev.resize(kDim + 1);
  for (int i = 0; i <= kDim; ++i) {
    instance.v_prev[i] = normal(rng);
  }
  return instance;
}

std::array<gblp::GlobalBlockInput, kBlocks> makeInputs(
    const Instance& instance, bool product_weights) {
  std::array<gblp::GlobalBlockInput, kBlocks> inputs;
  for (int b = 0; b < kBlocks; ++b) {
    inputs[b].coords = kCoords[b];
    inputs[b].coord_count = kCoordCount[b];
    inputs[b].rows = instance.rows[b];
    inputs[b].objective_weight
        = product_weights ? static_cast<double>(instance.product_rows)
                                / static_cast<double>(kRowsPerBlock[b])
                          : 1.0;
  }
  return inputs;
}

struct DenseLp {
  std::vector<std::vector<double>> rows;
  std::vector<double> rhs;
  std::vector<double> objective;
  int num_cols = 0;
};

// The product form: one feasibility row per element of R_A x R_B x R_C, over
// the 17 columns [V | v | s], exactly as GlobalAffineApproximator assembles it
// today.
DenseLp buildProductLp(const Instance& instance, double sigma) {
  constexpr int kProductCols = 2 * kDim + 1;
  DenseLp lp;
  lp.num_cols = kProductCols;
  lp.objective.assign(kProductCols, 0.0);

  int vertex_count = 0;
  for (const auto& a : instance.rows[0]) {
    for (const auto& b : instance.rows[1]) {
      for (const auto& c : instance.rows[2]) {
        const std::array<const gblp::GlobalBlockRow*, kBlocks> parts
            = {&a, &b, &c};
        std::vector<double> row(kProductCols, 0.0);
        double kappa = 0.0;
        for (int blk = 0; blk < kBlocks; ++blk) {
          for (int d = 0; d < kCoordCount[blk]; ++d) {
            const int axis = kCoords[blk][d];
            row[axis] += sigma * parts[blk]->p(d);
            row[kDim + 1 + axis] += kDt * parts[blk]->r(d);
          }
          kappa += parts[blk]->kappa;
        }
        row[kDim] += -sigma;
        lp.rows.push_back(std::move(row));
        lp.rhs.push_back(-sigma * kappa);

        for (int blk = 0; blk < kBlocks; ++blk) {
          for (int d = 0; d < kCoordCount[blk]; ++d) {
            lp.objective[kCoords[blk][d]] += -sigma * parts[blk]->p(d);
          }
        }
        lp.objective[kDim] += sigma;
        ++vertex_count;
      }
    }
  }
  EXPECT_EQ(vertex_count, instance.product_rows);

  for (int sign_index = 0; sign_index < 2; ++sign_index) {
    const double sign = sign_index == 0 ? 1.0 : -1.0;
    for (int i = 0; i < kDim; ++i) {
      std::vector<double> row(kProductCols, 0.0);
      row[i] = sign;
      row[kDim + 1 + i] = -1.0;
      lp.rows.push_back(std::move(row));
      lp.rhs.push_back(0.0);
    }
  }
  for (int i = 0; i < kDim; ++i) {
    lp.objective[kDim + 1 + i] += kSPin;
  }
  return lp;
}

// The per-step right-hand side of the product form, for comparison against
// globalReducedLpRowUpper().
std::vector<double> productRowUpper(const Instance& instance, double sigma,
                                    const DenseLp& lp) {
  std::vector<double> upper = lp.rhs;
  int row_id = 0;
  for (const auto& a : instance.rows[0]) {
    for (const auto& b : instance.rows[1]) {
      for (const auto& c : instance.rows[2]) {
        const std::array<const gblp::GlobalBlockRow*, kBlocks> parts
            = {&a, &b, &c};
        double b_upd = instance.v_prev[kDim];
        for (int blk = 0; blk < kBlocks; ++blk) {
          for (int d = 0; d < kCoordCount[blk]; ++d) {
            b_upd += instance.v_prev[kCoords[blk][d]] * parts[blk]->nu(d);
          }
        }
        upper[row_id] = lp.rhs[row_id] - sigma * b_upd;
        ++row_id;
      }
    }
  }
  return upper;
}

struct SolveResult {
  bool optimal = false;
  double objective_value = 0.0;
  std::vector<double> solution;
};

SolveResult solveCsr(const std::vector<int>& starts,
                     const std::vector<int>& cols,
                     const std::vector<double>& values,
                     const std::vector<double>& row_lower,
                     const std::vector<double>& row_upper,
                     const std::vector<double>& objective,
                     std::vector<double> col_lower,
                     std::vector<double> col_upper) {
  Highs highs;
  highs.setOptionValue("solver", "simplex");
  highs.setOptionValue("presolve", "on");
  highs.setOptionValue("primal_feasibility_tolerance", 1e-9);
  highs.setOptionValue("dual_feasibility_tolerance", 1e-9);
  highs.setOptionValue("log_to_console", false);
  highs.changeObjectiveSense(ObjSense::kMinimize);

  const int n = static_cast<int>(objective.size());
  const int m = static_cast<int>(row_upper.size());
  if (highs.addCols(n, objective.data(), col_lower.data(), col_upper.data(), 0,
                    nullptr, nullptr, nullptr)
      != HighsStatus::kOk) {
    return {};
  }
  if (highs.addRows(m, row_lower.data(), row_upper.data(),
                    static_cast<int>(values.size()), starts.data(), cols.data(),
                    values.data())
      != HighsStatus::kOk) {
    return {};
  }
  if (highs.run() != HighsStatus::kOk) {
    return {};
  }
  SolveResult result;
  result.optimal = highs.getModelStatus() == HighsModelStatus::kOptimal;
  if (result.optimal) {
    result.objective_value = highs.getObjectiveValue();
    result.solution = highs.getSolution().col_value;
  }
  return result;
}

SolveResult solveDense(const DenseLp& lp, const std::vector<double>& upper) {
  std::vector<int> starts = {0};
  std::vector<int> cols;
  std::vector<double> values;
  for (const auto& row : lp.rows) {
    for (int c = 0; c < lp.num_cols; ++c) {
      if (row[c] != 0.0) {
        cols.push_back(c);
        values.push_back(row[c]);
      }
    }
    starts.push_back(static_cast<int>(values.size()));
  }
  std::vector<double> row_lower(lp.rows.size(),
                                -std::numeric_limits<double>::infinity());
  return solveCsr(starts, cols, values, row_lower, upper, lp.objective,
                  std::vector<double>(lp.num_cols, -kBox),
                  std::vector<double>(lp.num_cols, kBox));
}

SolveResult solveReduced(const gblp::GlobalReducedLp& lp,
                         const std::vector<double>& upper) {
  std::vector<double> objective(lp.objective.data(),
                                lp.objective.data() + lp.objective.size());
  std::vector<double> col_lower(lp.num_cols, -kBox);
  std::vector<double> col_upper(lp.num_cols, kBox);
  // mu is free; boxing it would change the problem.
  for (int b = 0; b < kBlocks; ++b) {
    col_lower[gblp::idxMu(b)] = -std::numeric_limits<double>::infinity();
    col_upper[gblp::idxMu(b)] = std::numeric_limits<double>::infinity();
  }
  return solveCsr(lp.starts, lp.cols, lp.values, lp.row_lower, upper, objective,
                  col_lower, col_upper);
}

}  // namespace

TEST(global_block_reduction, reduced_matches_product_with_product_weights) {
  for (double sigma : {-1.0, 1.0}) {
    const Instance instance = makeInstance(/*seed=*/11);

    const DenseLp product = buildProductLp(instance, sigma);
    const std::vector<double> product_upper
        = productRowUpper(instance, sigma, product);

    gblp::GlobalReducedLpOptions options;
    options.dt = kDt;
    options.sigma = sigma;
    options.s_pin_weight = kSPin;
    const gblp::GlobalReducedLp reduced = assembleGlobalReducedLp(
        makeInputs(instance, /*product_weights=*/true), options);
    const std::vector<double> reduced_upper
        = gblp::globalReducedLpRowUpper(reduced, instance.v_prev);

    // 36 product rows + 16 s rows against 10 block rows + 1 coupling + 16.
    EXPECT_EQ(product.rows.size(),
              static_cast<std::size_t>(instance.product_rows) + 2 * kDim)
        << "sigma=" << sigma;
    EXPECT_EQ(reduced_upper.size(),
              static_cast<std::size_t>(kRowsPerBlock[0] + kRowsPerBlock[1]
                                       + kRowsPerBlock[2])
                  + 1 + 2 * kDim)
        << "sigma=" << sigma;
    EXPECT_EQ(product.num_cols, 17) << "sigma=" << sigma;
    EXPECT_EQ(reduced.num_cols, 20) << "sigma=" << sigma;

    // The objective over the shared columns must be reproduced exactly, not
    // merely give the same optimum.
    double objective_gap = 0.0;
    for (int c = 0; c < product.num_cols; ++c) {
      objective_gap
          = std::max(objective_gap, std::abs(product.objective[c]
                                             - reduced.objective(c)));
    }

    const SolveResult product_result = solveDense(product, product_upper);
    const SolveResult reduced_result = solveReduced(reduced, reduced_upper);
    ASSERT_TRUE(product_result.optimal) << "sigma=" << sigma;
    ASSERT_TRUE(reduced_result.optimal) << "sigma=" << sigma;

    double solution_gap = 0.0;
    for (int c = 0; c < product.num_cols; ++c) {
      solution_gap = std::max(
          solution_gap,
          std::abs(product_result.solution[c] - reduced_result.solution[c]));
    }
    const double value_gap = std::abs(product_result.objective_value
                                      - reduced_result.objective_value);

    GTEST_COUT << std::format(
        "sigma={:+.0f}  rows {}/{}  cols {}/{}  |dc|_inf={:.2e}  "
        "|dvalue|={:.2e}  |d[V,v,s]|_inf={:.2e}\n",
        sigma, product.rows.size(), reduced_upper.size(), product.num_cols,
        reduced.num_cols, objective_gap, value_gap, solution_gap);

    EXPECT_LT(objective_gap, 1e-9) << "objectives disagree, sigma=" << sigma;
    EXPECT_NEAR(product_result.objective_value,
                reduced_result.objective_value, 1e-7)
        << "optimal values disagree, sigma=" << sigma;
    EXPECT_LT(solution_gap, 1e-6)
        << "optimal [V, v, s] disagree, sigma=" << sigma;
  }
}

TEST(global_block_reduction, normalised_weights_stay_feasible_in_product) {
  for (double sigma : {-1.0, 1.0}) {
    const Instance instance = makeInstance(/*seed=*/11);

    gblp::GlobalReducedLpOptions options;
    options.dt = kDt;
    options.sigma = sigma;
    options.s_pin_weight = kSPin;
    const gblp::GlobalReducedLp reduced = assembleGlobalReducedLp(
        makeInputs(instance, /*product_weights=*/false), options);
    const SolveResult reduced_result = solveReduced(
        reduced, gblp::globalReducedLpRowUpper(reduced, instance.v_prev));
    ASSERT_TRUE(reduced_result.optimal) << "sigma=" << sigma;

    // The normalised optimum solves a different problem, but it must still be
    // feasible in the product LP or it is not a certificate.
    const DenseLp product = buildProductLp(instance, sigma);
    const std::vector<double> product_upper
        = productRowUpper(instance, sigma, product);
    double worst = 0.0;
    for (std::size_t r = 0; r < product.rows.size(); ++r) {
      double lhs = 0.0;
      for (int c = 0; c < product.num_cols; ++c) {
        lhs += product.rows[r][c] * reduced_result.solution[c];
      }
      worst = std::max(worst, lhs - product_upper[r]);
    }
    EXPECT_LT(worst, 1e-6)
        << "normalised-weight optimum is infeasible in the product LP, sigma="
        << sigma;
  }
}
