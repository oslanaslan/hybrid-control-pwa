// Tier-1 verification of the block reduction (ctx part VII, tier 1).
//
// This is the ONLY correctness oracle available for the reduction. The
// production product LP has ~4e8 rows and is not solvable with the resources
// at hand, and there is no other implementation to diff against, so the
// equivalence has to be established on synthetic data of the right shape and
// then relied on. Accordingly this test lands before any production code
// changes, and the reduced LP it exercises is the same assembler that
// prepareLpMatrices() will later call -- not a parallel implementation.
//
// It is a port of notes/check_reduction_full.py from the paper repository.
// The original LP is assembled here exactly as written in step2_1_main_lp.md
// section 4.1a: rows (R) and (L) over the full Cartesian product
// R_A x R_B x R_C, an independent y_j in R^8 for every full region j, and 16
// rows (Y) per full region. Reference results from the Python version:
//
//                          lower (s=-1)      upper (s=+1)
//   rows, original/reduced   936 / 90         936 / 90
//   vars, original/reduced   159 / 43         159 / 43
//   ||c_orig - c_red||_inf   4.2e-12          4.2e-12
//   |value_orig-value_red|   1.7e-11          0.0
//   ||z*_orig - z*_red||_inf 2.8e-14          7.4e-14
//
// THE TRAP. The equivalence holds only with the product weights R/R_block on
// the objective (step 7 section 9). Production deliberately uses normalised
// weights of 1 (ctx VI.1), which is a DIFFERENT problem with a different z*.
// A version of this test written with normalised weights fails on z* and looks
// like a broken reduction. So the equivalence is asserted with R/R_block, and
// the normalised setting is checked for the property production actually
// relies on: that its optimum is still feasible in the product LP, hence still
// a valid certificate. Validity does not depend on the objective; only
// tightness does.

#include <gtest/gtest.h>

#include <Highs.h>

#include <Eigen/Core>
#include <array>
#include <cmath>
#include <format>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

#include "test_utils.hpp"
#include "util/block_reduction_lp.hpp"

namespace {

namespace blp = hcpwa::util::block_lp;

constexpr int kBlocks = blp::kBlockCount;
constexpr int kFullDim = 8;
constexpr double kDt = 0.25;
// Random phi and Psi do not satisfy the invariants 1^T phi = 5 and Psi 1 = 0
// that the feasibility argument of step 7 section 11 rests on, so the synthetic
// problem needs an artificial box on z to stay bounded. The Python reference
// uses the same box.
constexpr double kZBox = 5.0;

// Shapes, matching check_reduction_full.py exactly.
const std::array<std::vector<int>, kBlocks> kCoords
    = {std::vector<int>{0, 1, 2}, std::vector<int>{3, 4, 5},
       std::vector<int>{6, 7}};
constexpr std::array<int, kBlocks> kEta = {6, 5, 4};
constexpr std::array<int, kBlocks> kNumRegions = {3, 3, 2};
constexpr std::array<int, kBlocks> kVertsPerRegion = {3, 3, 2};

int etaTotal() { return kEta[0] + kEta[1] + kEta[2]; }

int etaOffset(int block) {
  int offset = 0;
  for (int b = 0; b < block; ++b) {
    offset += kEta[b];
  }
  return offset;
}

struct SyntheticRow {
  int block_region = 0;
  Eigen::VectorXd phi;  // length kEta[block]
  Eigen::VectorXd m;    // length |kCoords[block]|
  Eigen::VectorXd rho;  // length |kCoords[block]|, >= 0
  double g = 0.0;
};

struct SyntheticBlock {
  std::vector<Eigen::MatrixXd> psi;  // [j] -> coord_count x eta_block
  std::vector<SyntheticRow> rows;
};

struct SyntheticInstance {
  std::array<SyntheticBlock, kBlocks> blocks;
  Eigen::VectorXd x_next;  // known value at the right endpoint
  std::array<int, kBlocks> num_rows{};
  int product_rows = 0;
  int full_regions = 0;
};

SyntheticInstance makeInstance(std::uint64_t seed) {
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  std::normal_distribution<double> normal(0.0, 1.0);

  SyntheticInstance instance;
  for (int b = 0; b < kBlocks; ++b) {
    const int dim = static_cast<int>(kCoords[b].size());
    SyntheticBlock& block = instance.blocks[b];

    block.psi.reserve(kNumRegions[b]);
    for (int j = 0; j < kNumRegions[b]; ++j) {
      Eigen::MatrixXd psi(dim, kEta[b]);
      for (int r = 0; r < dim; ++r) {
        for (int c = 0; c < kEta[b]; ++c) {
          psi(r, c) = normal(rng);
        }
      }
      block.psi.push_back(std::move(psi));
    }

    for (int j = 0; j < kNumRegions[b]; ++j) {
      for (int v = 0; v < kVertsPerRegion[b]; ++v) {
        SyntheticRow row;
        row.block_region = j;
        row.phi = Eigen::VectorXd(kEta[b]);
        for (int c = 0; c < kEta[b]; ++c) {
          row.phi(c) = uniform(rng);
        }
        row.m = Eigen::VectorXd(dim);
        row.rho = Eigen::VectorXd(dim);
        for (int r = 0; r < dim; ++r) {
          row.m(r) = normal(rng);
          row.rho(r) = uniform(rng);  // rho >= 0 by construction
        }
        row.g = normal(rng);
        block.rows.push_back(std::move(row));
      }
    }
    instance.num_rows[b] = static_cast<int>(block.rows.size());
  }

  instance.product_rows
      = instance.num_rows[0] * instance.num_rows[1] * instance.num_rows[2];
  instance.full_regions
      = kNumRegions[0] * kNumRegions[1] * kNumRegions[2];

  instance.x_next = Eigen::VectorXd(etaTotal());
  for (int i = 0; i < etaTotal(); ++i) {
    instance.x_next(i) = normal(rng);
  }
  return instance;
}

// Decomposition of one block row's contribution, following blk_row_parts() in
// the Python reference.
struct RowParts {
  int block_region = 0;
  Eigen::VectorXd a_right;  // coefficients on z_block in (R), = -phi/dt
  double c_right = 0.0;     // phi^T x/dt + beta
  Eigen::VectorXd a_left;   // coefficients on z_block in (L), = Psi^T m - phi/dt
  double c_left = 0.0;      // phi^T x/dt + g
  Eigen::VectorXd rho;
};

RowParts rowParts(const SyntheticInstance& instance, int block,
                  const SyntheticRow& row, double s) {
  const SyntheticBlock& blk = instance.blocks[block];
  const int lo = etaOffset(block);
  const Eigen::VectorXd x_block = instance.x_next.segment(lo, kEta[block]);
  const Eigen::MatrixXd& psi = blk.psi[row.block_region];

  const Eigen::VectorXd q = psi * x_block;
  const double beta = q.dot(row.m) + s * row.rho.dot(q.cwiseAbs()) + row.g;
  const double phi_x = row.phi.dot(x_block);

  RowParts parts;
  parts.block_region = row.block_region;
  parts.a_right = -row.phi / kDt;
  parts.c_right = phi_x / kDt + beta;
  parts.a_left = psi.transpose() * row.m - row.phi / kDt;
  parts.c_left = phi_x / kDt + row.g;
  parts.rho = row.rho;
  return parts;
}

// A dense linear program, used only for the product form. At 936 x 159 the
// dense representation is a rounding error in memory and much easier to audit
// than CSR.
struct DenseLp {
  std::vector<std::vector<double>> rows;
  std::vector<double> rhs;
  std::vector<double> objective;
  std::vector<double> col_lower;
  std::vector<double> col_upper;
  int num_cols = 0;
};

// The original LP of step2_1_main_lp.md section 4.1a: rows over the full
// Cartesian product, one independent y_j in R^8 per full region.
DenseLp buildProductLp(const SyntheticInstance& instance, double s) {
  const int eta_total = etaTotal();
  const int num_cols = eta_total + kFullDim * instance.full_regions;

  DenseLp lp;
  lp.num_cols = num_cols;
  lp.objective.assign(num_cols, 0.0);
  lp.col_lower.assign(num_cols, 0.0);
  lp.col_upper.assign(num_cols, 0.0);
  for (int c = 0; c < eta_total; ++c) {
    lp.col_lower[c] = -kZBox;
    lp.col_upper[c] = kZBox;
  }
  for (int c = eta_total; c < num_cols; ++c) {
    lp.col_lower[c] = 0.0;
    lp.col_upper[c] = std::numeric_limits<double>::infinity();
  }

  std::array<std::vector<RowParts>, kBlocks> parts;
  for (int b = 0; b < kBlocks; ++b) {
    parts[b].reserve(instance.blocks[b].rows.size());
    for (const SyntheticRow& row : instance.blocks[b].rows) {
      parts[b].push_back(rowParts(instance, b, row, s));
    }
  }

  // y offset of a full region, in the region ordering (j_A, j_B, j_C).
  auto y_offset = [&](int ja, int jb, int jc) {
    const int index = (ja * kNumRegions[1] + jb) * kNumRegions[2] + jc;
    return eta_total + kFullDim * index;
  };

  for (const RowParts& pa : parts[0]) {
    for (const RowParts& pb : parts[1]) {
      for (const RowParts& pc : parts[2]) {
        const std::array<const RowParts*, kBlocks> p = {&pa, &pb, &pc};

        // (R): s F^R <= 0.
        std::vector<double> right(num_cols, 0.0);
        double right_const = 0.0;
        for (int b = 0; b < kBlocks; ++b) {
          const int lo = etaOffset(b);
          for (int c = 0; c < kEta[b]; ++c) {
            right[lo + c] += s * p[b]->a_right(c);
          }
          right_const += p[b]->c_right;
        }
        lp.rows.push_back(std::move(right));
        lp.rhs.push_back(-s * right_const);

        // Objective: sum of delta^R = -s sum F^R over the whole product.
        for (int b = 0; b < kBlocks; ++b) {
          const int lo = etaOffset(b);
          for (int c = 0; c < kEta[b]; ++c) {
            lp.objective[lo + c] += -s * p[b]->a_right(c);
          }
        }

        // (L): s(a^T z + c) + rho^T y_j <= 0, with y_j the full region's own
        // 8-vector.
        std::vector<double> left(num_cols, 0.0);
        double left_const = 0.0;
        const int yoff = y_offset(pa.block_region, pb.block_region,
                                  pc.block_region);
        for (int b = 0; b < kBlocks; ++b) {
          const int lo = etaOffset(b);
          for (int c = 0; c < kEta[b]; ++c) {
            left[lo + c] += s * p[b]->a_left(c);
          }
          left_const += p[b]->c_left;
          for (std::size_t t = 0; t < kCoords[b].size(); ++t) {
            left[yoff + kCoords[b][t]] += p[b]->rho(static_cast<int>(t));
          }
        }
        lp.rows.push_back(std::move(left));
        lp.rhs.push_back(-s * left_const);
      }
    }
  }

  // (Y): y_j >= +/- Psi_j z, 16 rows per full region.
  for (int ja = 0; ja < kNumRegions[0]; ++ja) {
    for (int jb = 0; jb < kNumRegions[1]; ++jb) {
      for (int jc = 0; jc < kNumRegions[2]; ++jc) {
        const int yoff = y_offset(ja, jb, jc);
        const std::array<int, kBlocks> region_ids = {ja, jb, jc};
        for (int b = 0; b < kBlocks; ++b) {
          const int lo = etaOffset(b);
          const Eigen::MatrixXd& psi
              = instance.blocks[b].psi[region_ids[b]];
          for (std::size_t t = 0; t < kCoords[b].size(); ++t) {
            for (double sign : {1.0, -1.0}) {
              std::vector<double> row(num_cols, 0.0);
              for (int c = 0; c < kEta[b]; ++c) {
                row[lo + c] = sign * psi(static_cast<int>(t), c);
              }
              row[yoff + kCoords[b][t]] = -1.0;
              lp.rows.push_back(std::move(row));
              lp.rhs.push_back(0.0);
            }
          }
        }
      }
    }
  }

  return lp;
}

std::array<blp::BlockInput, kBlocks> buildBlockInputs(
    const SyntheticInstance& instance, bool product_weights) {
  std::array<blp::BlockInput, kBlocks> blocks;
  for (int b = 0; b < kBlocks; ++b) {
    const int lo = etaOffset(b);
    const int dim = static_cast<int>(kCoords[b].size());
    blp::BlockInput& out = blocks[b];
    out.coord_count = dim;
    out.num_block_regions = kNumRegions[b];
    out.objective_weight
        = product_weights
              ? static_cast<double>(instance.product_rows)
                    / static_cast<double>(instance.num_rows[b])
              : 1.0;

    out.psi.reserve(kNumRegions[b]);
    for (int j = 0; j < kNumRegions[b]; ++j) {
      blp::BlockPsi psi;
      psi.rows.resize(dim);
      for (int r = 0; r < dim; ++r) {
        for (int c = 0; c < kEta[b]; ++c) {
          psi.rows[r].add(lo + c, instance.blocks[b].psi[j](r, c));
        }
      }
      out.psi.push_back(std::move(psi));
    }

    out.rows.reserve(instance.blocks[b].rows.size());
    for (const SyntheticRow& row : instance.blocks[b].rows) {
      blp::BlockRow block_row;
      block_row.block_region = row.block_region;
      for (int c = 0; c < kEta[b]; ++c) {
        block_row.phi.add(lo + c, row.phi(c));
      }
      block_row.m = row.m;
      block_row.rho = row.rho;
      block_row.g = row.g;
      out.rows.push_back(std::move(block_row));
    }
  }
  return blocks;
}

struct SolveResult {
  bool optimal = false;
  double objective_value = 0.0;
  std::vector<double> solution;
  int num_rows = 0;
  int num_cols = 0;
};

void configure(Highs& highs) {
  highs.setOptionValue("solver", "simplex");
  highs.setOptionValue("presolve", "on");
  highs.setOptionValue("primal_feasibility_tolerance", 1e-9);
  highs.setOptionValue("dual_feasibility_tolerance", 1e-9);
  highs.setOptionValue("log_to_console", false);
  highs.changeObjectiveSense(ObjSense::kMinimize);
}

SolveResult solveCsr(const std::vector<int>& starts,
                     const std::vector<int>& cols,
                     const std::vector<double>& values,
                     const std::vector<double>& row_lower,
                     const std::vector<double>& row_upper,
                     const std::vector<double>& objective,
                     const std::vector<double>& col_lower,
                     const std::vector<double>& col_upper) {
  Highs highs;
  configure(highs);

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
  result.num_rows = m;
  result.num_cols = n;
  result.optimal = highs.getModelStatus() == HighsModelStatus::kOptimal;
  if (result.optimal) {
    result.objective_value = highs.getObjectiveValue();
    result.solution = highs.getSolution().col_value;
  }
  return result;
}

SolveResult solveDense(const DenseLp& lp) {
  std::vector<int> starts = {0};
  std::vector<int> cols;
  std::vector<double> values;
  std::vector<double> row_lower(lp.rows.size(),
                                -std::numeric_limits<double>::infinity());
  for (const std::vector<double>& row : lp.rows) {
    for (int c = 0; c < lp.num_cols; ++c) {
      if (row[c] != 0.0) {
        cols.push_back(c);
        values.push_back(row[c]);
      }
    }
    starts.push_back(static_cast<int>(values.size()));
  }
  return solveCsr(starts, cols, values, row_lower, lp.rhs, lp.objective,
                  lp.col_lower, lp.col_upper);
}

SolveResult solveReduced(const blp::ReducedLp& lp,
                         const std::vector<double>& row_upper) {
  std::vector<double> objective(lp.objective.data(),
                                lp.objective.data() + lp.objective.size());
  std::vector<double> col_lower = lp.col_lower;
  std::vector<double> col_upper = lp.col_upper;
  // Same artificial box on z as the product form, for the same reason.
  for (int c = 0; c < lp.num_x; ++c) {
    col_lower[c] = -kZBox;
    col_upper[c] = kZBox;
  }
  return solveCsr(lp.starts, lp.cols, lp.values, lp.row_lower, row_upper,
                  objective, col_lower, col_upper);
}

double maxAbsDiff(const std::vector<double>& lhs,
                  const std::vector<double>& rhs, int count) {
  double worst = 0.0;
  for (int i = 0; i < count; ++i) {
    worst = std::max(worst, std::abs(lhs[i] - rhs[i]));
  }
  return worst;
}

}  // namespace

// The decisive test: both endpoints of the interval, the identification of y
// (Lemma 5), and both bound directions.
TEST(barycentric_block_reduction, reduced_matches_product_with_product_weights) {
  for (double s : {-1.0, 1.0}) {
    const SyntheticInstance instance = makeInstance(/*seed=*/7);

    const DenseLp product = buildProductLp(instance, s);
    const auto block_inputs
        = buildBlockInputs(instance, /*product_weights=*/true);
    blp::ReducedLpOptions options;
    options.dt = kDt;
    options.s = s;
    options.num_x = etaTotal();
    const blp::ReducedLp reduced = assembleReducedLp(block_inputs, options);
    const std::vector<double> row_upper = blp::reducedLpRowUpper(
        reduced, std::vector<double>(instance.x_next.data(),
                                     instance.x_next.data()
                                         + instance.x_next.size()));

    // Structural counts. These depend only on the shapes, so they pin the
    // assembler independently of any numerics.
    EXPECT_EQ(product.rows.size(), 936U) << "s=" << s;
    EXPECT_EQ(product.num_cols, 159) << "s=" << s;
    EXPECT_EQ(row_upper.size(), 90U) << "s=" << s;
    EXPECT_EQ(reduced.num_cols, 43) << "s=" << s;

    const SolveResult product_result = solveDense(product);
    const SolveResult reduced_result = solveReduced(reduced, row_upper);
    ASSERT_TRUE(product_result.optimal) << "product LP not optimal, s=" << s;
    ASSERT_TRUE(reduced_result.optimal) << "reduced LP not optimal, s=" << s;

    // The objective must be reproduced exactly, not just its optimal value
    // (step 7 section 9). Python reference: 4.2e-12.
    std::vector<double> reduced_objective(
        reduced.objective.data(),
        reduced.objective.data() + reduced.objective.size());
    const double objective_gap
        = maxAbsDiff(product.objective, reduced_objective, etaTotal());

    // Python reference: 1.7e-11 (lower), 0.0 (upper).
    const double value_gap = std::abs(product_result.objective_value
                                      - reduced_result.objective_value);

    // Python reference: 2.8e-14 (lower), 7.4e-14 (upper). Under the product
    // weights the two problems share an optimal face; a failure here that
    // survives the value check above would be genuine non-uniqueness rather
    // than an error in the reduction (step 7 section 12.2), but on this fixed
    // instance the solver picks the same vertex.
    const double z_gap = maxAbsDiff(product_result.solution,
                                    reduced_result.solution, etaTotal());

    GTEST_COUT << std::format(
        "s={:+.0f}  rows {}/{}  vars {}/{}  |dc|_inf={:.2e}  |dvalue|={:.2e}  "
        "|dz|_inf={:.2e}\n",
        s, product.rows.size(), row_upper.size(), product.num_cols,
        reduced.num_cols, objective_gap, value_gap, z_gap);

    EXPECT_LT(objective_gap, 1e-9) << "objective vectors disagree, s=" << s;
    EXPECT_NEAR(product_result.objective_value,
                reduced_result.objective_value, 1e-8)
        << "optimal values disagree, s=" << s;
    EXPECT_LT(z_gap, 1e-7) << "optimal z disagree, s=" << s;
  }
}

// What production actually relies on: with normalised weights the reduced
// optimum is a different point, but it is still feasible in the product LP,
// hence still a valid certificate (ctx VI.1).
TEST(barycentric_block_reduction, normalised_weights_stay_feasible_in_product) {
  for (double s : {-1.0, 1.0}) {
    const SyntheticInstance instance = makeInstance(/*seed=*/7);

    const auto block_inputs
        = buildBlockInputs(instance, /*product_weights=*/false);
    blp::ReducedLpOptions options;
    options.dt = kDt;
    options.s = s;
    options.num_x = etaTotal();
    const blp::ReducedLp reduced = assembleReducedLp(block_inputs, options);
    const std::vector<double> row_upper = blp::reducedLpRowUpper(
        reduced, std::vector<double>(instance.x_next.data(),
                                     instance.x_next.data()
                                         + instance.x_next.size()));
    const SolveResult reduced_result = solveReduced(reduced, row_upper);
    ASSERT_TRUE(reduced_result.optimal) << "reduced LP not optimal, s=" << s;

    // Rebuild the full point (z, y) using the canonical y of Lemma 5's proof,
    //   y_{j,r} = |(Psi_{blk(r), j_blk(r)} z_blk(r))_r|,
    // which is exactly the construction that makes the identification of y
    // safe (step 7 section 8).
    const DenseLp product = buildProductLp(instance, s);
    std::vector<double> point(product.num_cols, 0.0);
    for (int c = 0; c < etaTotal(); ++c) {
      point[c] = reduced_result.solution[c];
    }
    for (int ja = 0; ja < kNumRegions[0]; ++ja) {
      for (int jb = 0; jb < kNumRegions[1]; ++jb) {
        for (int jc = 0; jc < kNumRegions[2]; ++jc) {
          const int index = (ja * kNumRegions[1] + jb) * kNumRegions[2] + jc;
          const int yoff = etaTotal() + kFullDim * index;
          const std::array<int, kBlocks> region_ids = {ja, jb, jc};
          for (int b = 0; b < kBlocks; ++b) {
            const int lo = etaOffset(b);
            Eigen::VectorXd z_block(kEta[b]);
            for (int c = 0; c < kEta[b]; ++c) {
              z_block(c) = point[lo + c];
            }
            const Eigen::VectorXd q
                = instance.blocks[b].psi[region_ids[b]] * z_block;
            for (std::size_t t = 0; t < kCoords[b].size(); ++t) {
              point[yoff + kCoords[b][t]]
                  = std::abs(q(static_cast<int>(t)));
            }
          }
        }
      }
    }

    double worst_violation = 0.0;
    for (std::size_t r = 0; r < product.rows.size(); ++r) {
      double lhs = 0.0;
      for (int c = 0; c < product.num_cols; ++c) {
        lhs += product.rows[r][c] * point[c];
      }
      worst_violation = std::max(worst_violation, lhs - product.rhs[r]);
    }
    EXPECT_LT(worst_violation, 1e-6)
        << "normalised-weight optimum is infeasible in the product LP, s=" << s;
  }
}

// Pins the variable layout of ctx IV.6 so a later change to it is a test
// failure rather than a silent reindexing.
TEST(barycentric_block_reduction, layout_is_x_then_blockwise_y_then_mu) {
  const SyntheticInstance instance = makeInstance(/*seed=*/7);
  const auto block_inputs = buildBlockInputs(instance, true);
  blp::ReducedLpOptions options;
  options.dt = kDt;
  options.s = 1.0;
  options.num_x = etaTotal();
  const blp::ReducedLp lp = assembleReducedLp(block_inputs, options);

  EXPECT_EQ(lp.num_x, 15);
  EXPECT_EQ(lp.y_offset[0], 15);
  EXPECT_EQ(lp.y_offset[1], 15 + 3 * 3);
  EXPECT_EQ(lp.y_offset[2], 15 + 3 * 3 + 3 * 3);
  EXPECT_EQ(lp.mu_r_offset, 15 + 22);
  EXPECT_EQ(lp.mu_l_offset, 15 + 22 + 3);
  EXPECT_EQ(lp.num_cols, 15 + 22 + 6);

  // y >= 0 as a column bound, mu free, x free (the caller owns gauge fixing).
  for (int b = 0; b < kBlocks; ++b) {
    EXPECT_EQ(lp.col_lower[lp.idxY(b, 0, 0)], 0.0);
  }
  for (int b = 0; b < kBlocks; ++b) {
    EXPECT_EQ(lp.col_lower[lp.idxMuR(b)],
              -std::numeric_limits<double>::infinity());
    EXPECT_EQ(lp.col_lower[lp.idxMuL(b)],
              -std::numeric_limits<double>::infinity());
  }
  EXPECT_EQ(lp.col_lower[0], -std::numeric_limits<double>::infinity());

  EXPECT_THROW((void)lp.idxY(0, kNumRegions[0], 0), std::invalid_argument);
  EXPECT_THROW((void)lp.idxY(0, 0, 3), std::invalid_argument);
}

// The tie-break is mandatory (step 7 section 12.2b): the optimal face is
// genuinely wider than a point, and the march is greedy, so an arbitrary choice
// among optima at one step changes every step after it. What must hold is that
// the perturbed problem still returns a point feasible in the original one, and
// that the objective moves by at most O(epsilon).
TEST(barycentric_block_reduction, tie_break_keeps_feasibility_and_value) {
  constexpr double kEpsilon = 1e-8;
  for (double s : {-1.0, 1.0}) {
    const SyntheticInstance instance = makeInstance(/*seed=*/7);
    const auto block_inputs
        = buildBlockInputs(instance, /*product_weights=*/true);
    const std::vector<double> x_next(instance.x_next.data(),
                                     instance.x_next.data()
                                         + instance.x_next.size());

    blp::ReducedLpOptions plain;
    plain.dt = kDt;
    plain.s = s;
    plain.num_x = etaTotal();
    const blp::ReducedLp lp_plain = assembleReducedLp(block_inputs, plain);
    const SolveResult plain_result
        = solveReduced(lp_plain, blp::reducedLpRowUpper(lp_plain, x_next));
    ASSERT_TRUE(plain_result.optimal) << "s=" << s;

    blp::ReducedLpOptions tied = plain;
    tied.tie_break_epsilon = kEpsilon;
    const blp::ReducedLp lp_tied = assembleReducedLp(block_inputs, tied);
    const SolveResult tied_result
        = solveReduced(lp_tied, blp::reducedLpRowUpper(lp_tied, x_next));
    ASSERT_TRUE(tied_result.optimal) << "s=" << s;

    // One auxiliary column and two rows per x column.
    EXPECT_EQ(lp_tied.num_cols, lp_plain.num_cols + etaTotal()) << "s=" << s;
    EXPECT_EQ(lp_tied.row_upper.size(),
              lp_plain.row_upper.size() + 2U * etaTotal())
        << "s=" << s;
    EXPECT_EQ(lp_tied.tie_break_offset, lp_plain.num_cols) << "s=" << s;
    // Offsets of everything else are untouched, so a caller holding an idxY or
    // idxMu from the un-tied layout stays correct.
    EXPECT_EQ(lp_tied.y_offset, lp_plain.y_offset) << "s=" << s;
    EXPECT_EQ(lp_tied.mu_r_offset, lp_plain.mu_r_offset) << "s=" << s;

    // The tie-break point must remain feasible in the ORIGINAL product LP,
    // which is the property that keeps it a valid certificate.
    const DenseLp product = buildProductLp(instance, s);
    std::vector<double> point(product.num_cols, 0.0);
    for (int c = 0; c < etaTotal(); ++c) {
      point[c] = tied_result.solution[c];
    }
    for (int ja = 0; ja < kNumRegions[0]; ++ja) {
      for (int jb = 0; jb < kNumRegions[1]; ++jb) {
        for (int jc = 0; jc < kNumRegions[2]; ++jc) {
          const int index = (ja * kNumRegions[1] + jb) * kNumRegions[2] + jc;
          const int yoff = etaTotal() + kFullDim * index;
          const std::array<int, kBlocks> region_ids = {ja, jb, jc};
          for (int b = 0; b < kBlocks; ++b) {
            const int lo = etaOffset(b);
            Eigen::VectorXd z_block(kEta[b]);
            for (int c = 0; c < kEta[b]; ++c) {
              z_block(c) = point[lo + c];
            }
            const Eigen::VectorXd q
                = instance.blocks[b].psi[region_ids[b]] * z_block;
            for (std::size_t t = 0; t < kCoords[b].size(); ++t) {
              point[yoff + kCoords[b][t]] = std::abs(q(static_cast<int>(t)));
            }
          }
        }
      }
    }
    double worst_violation = 0.0;
    for (std::size_t r = 0; r < product.rows.size(); ++r) {
      double lhs = 0.0;
      for (int c = 0; c < product.num_cols; ++c) {
        lhs += product.rows[r][c] * point[c];
      }
      worst_violation = std::max(worst_violation, lhs - product.rhs[r]);
    }
    EXPECT_LT(worst_violation, 1e-6)
        << "tie-break optimum is infeasible in the product LP, s=" << s;

    // The true objective at the tie-break point, excluding the epsilon term.
    double tied_true_objective = 0.0;
    for (int c = 0; c < etaTotal(); ++c) {
      tied_true_objective += lp_plain.objective(c) * tied_result.solution[c];
    }
    const double slack = tied_true_objective - plain_result.objective_value;
    GTEST_COUT << std::format(
        "s={:+.0f}  tie-break objective slack = {:.3e} (epsilon={:.0e})\n", s,
        slack, kEpsilon);
    EXPECT_GT(slack, -1e-6) << "tie-break beat the unperturbed optimum, s=" << s;
    EXPECT_LT(slack, 1e-3) << "tie-break moved the objective too far, s=" << s;
  }
}
