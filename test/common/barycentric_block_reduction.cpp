#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

#include <Eigen/Core>
#include <Highs.h>

#include <barycentric_geometry_types.hpp>
#include <util/block_reduction_lp.hpp>

// Tier 1 of the block reduction: the reduced LP against the LP written out over
// the full product of regions and vertices.
//
// The product assembler lives here and only here. It is the specification the
// reduction is checked against, and it is only ever run on a fixture small
// enough to solve -- on production geometry it would be about 4e8 rows.
namespace {

using barycentric_affine_approximator::SparseVec;
namespace br = barycentric_affine_approximator::block_reduction;

constexpr double kDt = 0.25;
constexpr double kBox = 5.0;
constexpr int kSeed = 7;

// Shape of the fixture. Sizes are the ones the reduction arithmetic was checked
// against by hand: 936 rows / 159 columns in product form, 90 / 43 reduced.
const std::vector<int> kCoords = {3, 3, 2};
const std::vector<int> kEta = {6, 5, 4};
const std::vector<int> kNumRegions = {3, 3, 2};
const std::vector<int> kNumVertices = {3, 3, 2};
// Number of projection layers each block owns; only used to give phi a
// realistic total mass.
const std::vector<int> kLayerCount = {2, 2, 1};

int totalNumX() { return std::accumulate(kEta.begin(), kEta.end(), 0); }

int blockColBegin(int block) {
  return std::accumulate(kEta.begin(), kEta.begin() + block, 0);
}

// Builds a random but structurally faithful reduced-LP input:
//   * phi is nonnegative, supported on the block's own x columns and sums to
//     the block's layer count, exactly like a sum of barycentric weights;
//   * every Psi row sums to zero over the block's columns, which is what makes
//     a uniform shift of x leave the gradient alone;
//   * rho is comfortably above the HiGHS small_matrix_value.
// The first two properties are also what keep the fixture feasible: pushing z
// along s * 1 drives every residual down without limit.
br::ReducedLpInput makeInput(double sign_s, br::ObjectiveWeights weights,
                             double tie_break_eps = 0.0) {
  std::mt19937 rng(kSeed);
  std::uniform_real_distribution<double> sym(-1.0, 1.0);
  std::uniform_real_distribution<double> pos(0.1, 1.0);
  std::uniform_real_distribution<double> mid(-2.0, 2.0);
  std::uniform_real_distribution<double> radius(0.05, 1.0);

  br::ReducedLpInput input;
  input.num_x = totalNumX();
  input.t_delta = kDt;
  input.sign_s = sign_s;
  input.weights = weights;
  input.tie_break_eps = tie_break_eps;

  for (std::size_t b = 0; b < kCoords.size(); ++b) {
    const int begin = blockColBegin(static_cast<int>(b));
    const int eta = kEta[b];

    br::ReducedLpBlock block;
    block.coord_count = kCoords[b];
    for (int j = 0; j < kNumRegions[b]; ++j) {
      br::BlockRegionData region;
      for (int p = 0; p < block.coord_count; ++p) {
        std::vector<double> raw(static_cast<std::size_t>(eta));
        double mean = 0.0;
        for (int c = 0; c < eta; ++c) {
          raw[static_cast<std::size_t>(c)] = sym(rng);
          mean += raw[static_cast<std::size_t>(c)];
        }
        mean /= static_cast<double>(eta);
        SparseVec row;
        for (int c = 0; c < eta; ++c) {
          row.add(begin + c, raw[static_cast<std::size_t>(c)] - mean,
                  br::kAssembleEps);
        }
        region.psi_rows.push_back(std::move(row));
      }

      for (int k = 0; k < kNumVertices[b]; ++k) {
        br::BlockVertexData vertex;
        std::vector<double> weight(static_cast<std::size_t>(eta));
        double total = 0.0;
        for (int c = 0; c < eta; ++c) {
          weight[static_cast<std::size_t>(c)] = pos(rng);
          total += weight[static_cast<std::size_t>(c)];
        }
        const double scale = static_cast<double>(kLayerCount[b]) / total;
        for (int c = 0; c < eta; ++c) {
          vertex.phi.add(begin + c, weight[static_cast<std::size_t>(c)] * scale,
                         br::kAssembleEps);
        }
        vertex.m = Eigen::VectorXd(block.coord_count);
        vertex.rho = Eigen::VectorXd(block.coord_count);
        for (int p = 0; p < block.coord_count; ++p) {
          vertex.m(p) = mid(rng);
          vertex.rho(p) = radius(rng);
        }
        vertex.g = sym(rng);
        region.vertices.push_back(std::move(vertex));
      }
      block.regions.push_back(std::move(region));
    }
    input.blocks.push_back(std::move(block));
  }
  return input;
}

std::vector<double> makeXNext() {
  std::mt19937 rng(kSeed + 101);
  std::uniform_real_distribution<double> sym(-1.0, 1.0);
  std::vector<double> x(static_cast<std::size_t>(totalNumX()));
  for (double& value : x) {
    value = sym(rng);
  }
  return x;
}

// ---------------------------------------------------------------------------
// The product-form reference.
// ---------------------------------------------------------------------------

struct ProductLp {
  std::vector<int> starts;
  std::vector<int> index;
  std::vector<double> value;
  std::vector<double> row_lower;
  std::vector<double> row_upper;
  std::vector<double> col_lower;
  std::vector<double> col_upper;
  Eigen::RowVectorXd cost;
  int num_rows = 0;
  int num_cols = 0;
};

int productStateDim() {
  return std::accumulate(kCoords.begin(), kCoords.end(), 0);
}

// Global state coordinate of (block, p) in the reference: blocks concatenated.
int productCoord(int block, int p) {
  return std::accumulate(kCoords.begin(), kCoords.begin() + block, 0) + p;
}

int productRegionCount() {
  int total = 1;
  for (int count : kNumRegions) {
    total *= count;
  }
  return total;
}

// Product region j <-> (jA, jB, jC), enumerated the way the geometry code nests
// its loops: j = (jA * MB + jB) * MC + jC.
std::vector<int> decodeRegion(int j) {
  std::vector<int> out(kNumRegions.size());
  for (int b = static_cast<int>(kNumRegions.size()) - 1; b >= 0; --b) {
    out[static_cast<std::size_t>(b)] = j % kNumRegions[b];
    j /= kNumRegions[b];
  }
  return out;
}

std::vector<int> decodeVertex(int v) {
  std::vector<int> out(kNumVertices.size());
  for (int b = static_cast<int>(kNumVertices.size()) - 1; b >= 0; --b) {
    out[static_cast<std::size_t>(b)] = v % kNumVertices[b];
    v /= kNumVertices[b];
  }
  return out;
}

int productVertexCount() {
  int total = 1;
  for (int count : kNumVertices) {
    total *= count;
  }
  return total;
}

int productIdxY(int num_x, int region, int coord) {
  return num_x + region * productStateDim() + coord;
}

void appendRow(ProductLp& lp, const SparseVec& row) {
  for (std::size_t i = 0; i < row.cols.size(); ++i) {
    lp.index.push_back(row.cols[i]);
    lp.value.push_back(row.vals[i]);
  }
  lp.starts.push_back(static_cast<int>(lp.value.size()));
}

ProductLp buildProductLp(const br::ReducedLpInput& input,
                         const std::vector<double>& x_next) {
  const int num_blocks = static_cast<int>(input.blocks.size());
  const int state_dim = productStateDim();
  const int num_regions = productRegionCount();
  const double s = input.sign_s;
  const double dt = input.t_delta;

  ProductLp lp;
  lp.num_cols = input.num_x + state_dim * num_regions;
  lp.cost = Eigen::RowVectorXd::Zero(lp.num_cols);
  lp.starts.assign(1, 0);

  for (int j = 0; j < num_regions; ++j) {
    const std::vector<int> jb = decodeRegion(j);

    // Psi_j: the block rows stacked in block order.
    std::vector<SparseVec> psi(static_cast<std::size_t>(state_dim));
    for (int b = 0; b < num_blocks; ++b) {
      const br::BlockRegionData& region
          = input.blocks[static_cast<std::size_t>(b)]
                .regions[static_cast<std::size_t>(jb[static_cast<std::size_t>(b)])];
      for (int p = 0; p < kCoords[static_cast<std::size_t>(b)]; ++p) {
        psi[static_cast<std::size_t>(productCoord(b, p))] = region.psi_rows[
            static_cast<std::size_t>(p)];
      }
    }

    Eigen::VectorXd q_next = Eigen::VectorXd::Zero(state_dim);
    for (int r = 0; r < state_dim; ++r) {
      q_next(r) = psi[static_cast<std::size_t>(r)].dot(x_next);
    }

    for (int v = 0; v < productVertexCount(); ++v) {
      const std::vector<int> vb = decodeVertex(v);

      SparseVec phi;
      Eigen::VectorXd m = Eigen::VectorXd::Zero(state_dim);
      Eigen::VectorXd rho = Eigen::VectorXd::Zero(state_dim);
      double g = 0.0;
      for (int b = 0; b < num_blocks; ++b) {
        const br::BlockVertexData& vertex
            = input.blocks[static_cast<std::size_t>(b)]
                  .regions[static_cast<std::size_t>(jb[static_cast<std::size_t>(b)])]
                  .vertices[static_cast<std::size_t>(vb[static_cast<std::size_t>(b)])];
        for (std::size_t i = 0; i < vertex.phi.cols.size(); ++i) {
          phi.add(vertex.phi.cols[i], vertex.phi.vals[i], br::kAssembleEps);
        }
        for (int p = 0; p < kCoords[static_cast<std::size_t>(b)]; ++p) {
          m(productCoord(b, p)) = vertex.m(p);
          rho(productCoord(b, p)) = vertex.rho(p);
        }
        g += vertex.g;
      }

      const double phi_x = phi.dot(x_next);

      // (L): s a^T z + rho^T y_j <= -s b,  a = Psi_j^T m - phi/dt,
      //      b = phi^T x_next / dt + g.
      SparseVec a;
      for (int r = 0; r < state_dim; ++r) {
        const SparseVec& row = psi[static_cast<std::size_t>(r)];
        for (std::size_t i = 0; i < row.cols.size(); ++i) {
          a.add(row.cols[i], row.vals[i] * m(r), br::kAssembleEps);
        }
      }
      for (std::size_t i = 0; i < phi.cols.size(); ++i) {
        a.add(phi.cols[i], -phi.vals[i] / dt, br::kAssembleEps);
      }

      SparseVec left;
      for (std::size_t i = 0; i < a.cols.size(); ++i) {
        left.add(a.cols[i], s * a.vals[i], br::kAssembleEps);
      }
      for (int r = 0; r < state_dim; ++r) {
        left.add(productIdxY(input.num_x, j, r), rho(r), 0.0);
      }
      appendRow(lp, left);
      lp.row_lower.push_back(-std::numeric_limits<double>::infinity());
      lp.row_upper.push_back(-s * (phi_x / dt + g));

      // (R): -(s/dt) phi^T z <= -s (phi^T x/dt + beta).
      const double beta
          = q_next.dot(m) + s * rho.dot(q_next.cwiseAbs()) + g;
      SparseVec right;
      for (std::size_t i = 0; i < phi.cols.size(); ++i) {
        right.add(phi.cols[i], -s * phi.vals[i] / dt, br::kAssembleEps);
      }
      appendRow(lp, right);
      lp.row_lower.push_back(-std::numeric_limits<double>::infinity());
      lp.row_upper.push_back(-s * (phi_x / dt + beta));

      for (std::size_t i = 0; i < phi.cols.size(); ++i) {
        lp.cost(phi.cols[i]) += (s / dt) * phi.vals[i];
      }
    }

    for (int r = 0; r < state_dim; ++r) {
      for (int sign_index = 0; sign_index < 2; ++sign_index) {
        const double scale = sign_index == 0 ? 1.0 : -1.0;
        SparseVec row;
        const SparseVec& psi_row = psi[static_cast<std::size_t>(r)];
        for (std::size_t i = 0; i < psi_row.cols.size(); ++i) {
          row.add(psi_row.cols[i], scale * psi_row.vals[i], br::kAssembleEps);
        }
        row.add(productIdxY(input.num_x, j, r), -1.0, 0.0);
        appendRow(lp, row);
        lp.row_lower.push_back(-std::numeric_limits<double>::infinity());
        lp.row_upper.push_back(0.0);
      }
    }
  }

  lp.num_rows = static_cast<int>(lp.row_upper.size());
  lp.col_lower.assign(static_cast<std::size_t>(lp.num_cols), -kBox);
  lp.col_upper.assign(static_cast<std::size_t>(lp.num_cols), kBox);
  for (int j = 0; j < num_regions; ++j) {
    for (int r = 0; r < state_dim; ++r) {
      const int col = productIdxY(input.num_x, j, r);
      lp.col_lower[static_cast<std::size_t>(col)] = 0.0;
      lp.col_upper[static_cast<std::size_t>(col)]
          = std::numeric_limits<double>::infinity();
    }
  }
  return lp;
}

// ---------------------------------------------------------------------------

struct LpResult {
  bool optimal = false;
  double objective = 0.0;
  std::vector<double> col_value;
};

LpResult solveLp(const std::vector<int>& starts, const std::vector<int>& index,
                 const std::vector<double>& value,
                 const std::vector<double>& row_lower,
                 const std::vector<double>& row_upper,
                 const std::vector<double>& col_lower,
                 const std::vector<double>& col_upper,
                 const Eigen::RowVectorXd& cost) {
  Highs highs;
  highs.setOptionValue("solver", "simplex");
  highs.setOptionValue("presolve", "on");
  highs.setOptionValue("primal_feasibility_tolerance", 1e-9);
  highs.setOptionValue("dual_feasibility_tolerance", 1e-9);
  highs.setOptionValue("small_matrix_value", br::kSmallMatrixValue);
  highs.setOptionValue("log_to_console", false);
  highs.changeObjectiveSense(ObjSense::kMinimize);

  const int n = static_cast<int>(cost.size());
  const int m = static_cast<int>(row_upper.size());
  if (highs.addCols(n, cost.data(), col_lower.data(), col_upper.data(), 0,
                    nullptr, nullptr, nullptr)
      != HighsStatus::kOk) {
    return LpResult{};
  }
  if (highs.addRows(m, row_lower.data(), row_upper.data(),
                    static_cast<int>(value.size()), starts.data(),
                    index.data(), value.data())
      != HighsStatus::kOk) {
    return LpResult{};
  }
  if (highs.run() != HighsStatus::kOk) {
    return LpResult{};
  }

  LpResult result;
  result.optimal = highs.getModelStatus() == HighsModelStatus::kOptimal;
  result.objective = highs.getInfo().objective_function_value;
  result.col_value = highs.getSolution().col_value;
  return result;
}

LpResult solveProduct(const ProductLp& lp) {
  return solveLp(lp.starts, lp.index, lp.value, lp.row_lower, lp.row_upper,
                 lp.col_lower, lp.col_upper, lp.cost);
}

LpResult solveReduced(const br::ReducedLpMatrices& lp,
                      const std::vector<double>& row_upper) {
  std::vector<double> col_lower = lp.col_lower;
  std::vector<double> col_upper = lp.col_upper;
  for (int k = 0; k < lp.cols.num_x; ++k) {
    col_lower[static_cast<std::size_t>(k)] = -kBox;
    col_upper[static_cast<std::size_t>(k)] = kBox;
  }
  return solveLp(lp.starts, lp.col_index, lp.value, lp.row_lower, row_upper,
                 col_lower, col_upper, lp.cost);
}

// Largest violation of a row of the LP by the given column values.
double maxRowViolation(const std::vector<int>& starts,
                       const std::vector<int>& index,
                       const std::vector<double>& value,
                       const std::vector<double>& row_upper,
                       const std::vector<double>& col_value) {
  double worst = 0.0;
  for (std::size_t r = 0; r + 1 < starts.size(); ++r) {
    double lhs = 0.0;
    for (int k = starts[r]; k < starts[r + 1]; ++k) {
      lhs += value[static_cast<std::size_t>(k)]
             * col_value[static_cast<std::size_t>(index[
                 static_cast<std::size_t>(k)])];
    }
    worst = std::max(worst, lhs - row_upper[r]);
  }
  return worst;
}

// Fills the product y block with |Psi_j z| so that a z taken from the reduced
// LP can be tested for feasibility in the product LP.
std::vector<double> liftToProduct(const br::ReducedLpInput& input,
                                  const ProductLp& lp,
                                  const std::vector<double>& z) {
  std::vector<double> col(static_cast<std::size_t>(lp.num_cols), 0.0);
  for (int k = 0; k < input.num_x; ++k) {
    col[static_cast<std::size_t>(k)] = z[static_cast<std::size_t>(k)];
  }
  const int state_dim = productStateDim();
  for (int j = 0; j < productRegionCount(); ++j) {
    const std::vector<int> jb = decodeRegion(j);
    for (int b = 0; b < static_cast<int>(input.blocks.size()); ++b) {
      const br::BlockRegionData& region
          = input.blocks[static_cast<std::size_t>(b)]
                .regions[static_cast<std::size_t>(jb[static_cast<std::size_t>(b)])];
      for (int p = 0; p < kCoords[static_cast<std::size_t>(b)]; ++p) {
        const double q = region.psi_rows[static_cast<std::size_t>(p)].dot(z);
        col[static_cast<std::size_t>(
            productIdxY(input.num_x, j, productCoord(b, p)))] = std::abs(q);
      }
    }
  }
  (void)state_dim;
  return col;
}

// Brute-force worst residual over the whole product, the definition that
// worstReducedResidual() claims to compute in linear time.
br::WorstResidual bruteForceWorst(const br::ReducedLpInput& input,
                                  const std::vector<double>& x_next,
                                  const std::vector<double>& z) {
  const int state_dim = productStateDim();
  const double s = input.sign_s;
  const double dt = input.t_delta;
  br::WorstResidual out;
  out.left = -std::numeric_limits<double>::infinity();
  out.right = -std::numeric_limits<double>::infinity();

  for (int j = 0; j < productRegionCount(); ++j) {
    const std::vector<int> jb = decodeRegion(j);
    Eigen::VectorXd q_left = Eigen::VectorXd::Zero(state_dim);
    Eigen::VectorXd q_right = Eigen::VectorXd::Zero(state_dim);
    for (int b = 0; b < static_cast<int>(input.blocks.size()); ++b) {
      const br::BlockRegionData& region
          = input.blocks[static_cast<std::size_t>(b)]
                .regions[static_cast<std::size_t>(jb[static_cast<std::size_t>(b)])];
      for (int p = 0; p < kCoords[static_cast<std::size_t>(b)]; ++p) {
        q_left(productCoord(b, p))
            = region.psi_rows[static_cast<std::size_t>(p)].dot(z);
        q_right(productCoord(b, p))
            = region.psi_rows[static_cast<std::size_t>(p)].dot(x_next);
      }
    }

    for (int v = 0; v < productVertexCount(); ++v) {
      const std::vector<int> vb = decodeVertex(v);
      SparseVec phi;
      Eigen::VectorXd m = Eigen::VectorXd::Zero(state_dim);
      Eigen::VectorXd rho = Eigen::VectorXd::Zero(state_dim);
      double g = 0.0;
      for (int b = 0; b < static_cast<int>(input.blocks.size()); ++b) {
        const br::BlockVertexData& vertex
            = input.blocks[static_cast<std::size_t>(b)]
                  .regions[static_cast<std::size_t>(jb[static_cast<std::size_t>(b)])]
                  .vertices[static_cast<std::size_t>(vb[static_cast<std::size_t>(b)])];
        for (std::size_t i = 0; i < vertex.phi.cols.size(); ++i) {
          phi.add(vertex.phi.cols[i], vertex.phi.vals[i], br::kAssembleEps);
        }
        for (int p = 0; p < kCoords[static_cast<std::size_t>(b)]; ++p) {
          m(productCoord(b, p)) = vertex.m(p);
          rho(productCoord(b, p)) = vertex.rho(p);
        }
        g += vertex.g;
      }
      const double slope = (phi.dot(x_next) - phi.dot(z)) / dt;
      const double f_left
          = slope + q_left.dot(m) + s * rho.dot(q_left.cwiseAbs()) + g;
      const double f_right
          = slope + q_right.dot(m) + s * rho.dot(q_right.cwiseAbs()) + g;
      out.left = std::max(out.left, s * f_left);
      out.right = std::max(out.right, s * f_right);
    }
  }
  return out;
}

std::vector<double> reducedZ(const br::ReducedLpMatrices& lp,
                             const LpResult& result) {
  return std::vector<double>(
      result.col_value.begin(),
      result.col_value.begin() + lp.cols.num_x);
}

}  // namespace

TEST(barycentric_block_reduction, layout_sizes_match_the_hand_check) {
  br::ReducedLpInput input
      = makeInput(1.0, br::ObjectiveWeights::ProductCount);
  br::clampBlockRho(input);
  const br::ReducedLpMatrices reduced = assembleReducedLp(input);

  EXPECT_EQ(reduced.rows.num_rows, 90);
  EXPECT_EQ(reduced.cols.num_cols, 43);

  const ProductLp product = buildProductLp(input, makeXNext());
  EXPECT_EQ(product.num_rows, 936);
  EXPECT_EQ(product.num_cols, 159);
}

TEST(barycentric_block_reduction, tie_break_block_only_appears_when_enabled) {
  br::ReducedLpInput input
      = makeInput(1.0, br::ObjectiveWeights::ProductCount);
  br::clampBlockRho(input);
  const br::ReducedLpMatrices plain = assembleReducedLp(input);
  EXPECT_EQ(plain.cols.u_offset, -1);
  EXPECT_EQ(plain.rows.tie_break_offset, -1);

  input.tie_break_eps = 1e-8;
  const br::ReducedLpMatrices lifted = assembleReducedLp(input);
  // The [x | y | mu] prefix must not move when the tie-break is switched on.
  EXPECT_EQ(lifted.cols.mu_l_offset, plain.cols.mu_l_offset);
  EXPECT_EQ(lifted.cols.num_cols, plain.cols.num_cols + plain.cols.num_x);
  EXPECT_EQ(lifted.rows.num_rows, plain.rows.num_rows + 2 * plain.cols.num_x);
}

TEST(barycentric_block_reduction, reduction_reproduces_the_product_lp) {
  const std::vector<double> x_next = makeXNext();

  for (const double s : {1.0, -1.0}) {
    br::ReducedLpInput input
        = makeInput(s, br::ObjectiveWeights::ProductCount);
    br::clampBlockRho(input);

    const br::ReducedLpMatrices reduced = assembleReducedLp(input);
    const std::vector<double> row_upper = br::updateReducedLpRowUpper(
        input, reduced.rows, reduced.row_upper, x_next);
    const ProductLp product = buildProductLp(input, x_next);

    // The objective is the same linear form. Only the x block can be nonzero.
    for (int k = 0; k < input.num_x; ++k) {
      EXPECT_NEAR(reduced.cost(k), product.cost(k), 1e-9)
          << "s=" << s << " column " << k;
    }
    for (int k = input.num_x; k < reduced.cols.num_cols; ++k) {
      EXPECT_EQ(reduced.cost(k), 0.0) << "s=" << s << " column " << k;
    }

    const LpResult product_result = solveProduct(product);
    const LpResult reduced_result = solveReduced(reduced, row_upper);
    ASSERT_TRUE(product_result.optimal) << "s=" << s;
    ASSERT_TRUE(reduced_result.optimal) << "s=" << s;

    EXPECT_NEAR(product_result.objective, reduced_result.objective, 1e-8)
        << "s=" << s;

    const std::vector<double> z_reduced = reducedZ(reduced, reduced_result);
    for (int k = 0; k < input.num_x; ++k) {
      EXPECT_NEAR(z_reduced[static_cast<std::size_t>(k)],
                  product_result.col_value[static_cast<std::size_t>(k)], 1e-7)
          << "s=" << s << " column " << k;
    }

    // Cross-feasibility: the reduced optimum, lifted with y = |Psi_j z|, must
    // satisfy every row of the product LP. This is the statement the reduction
    // actually makes, independently of which optimal vertex HiGHS picked.
    const std::vector<double> lifted
        = liftToProduct(input, product, z_reduced);
    EXPECT_LE(maxRowViolation(product.starts, product.index, product.value,
                              product.row_upper, lifted),
              1e-7)
        << "s=" << s;
  }
}

TEST(barycentric_block_reduction, unit_weights_still_produce_a_valid_bound) {
  const std::vector<double> x_next = makeXNext();

  for (const double s : {1.0, -1.0}) {
    br::ReducedLpInput input = makeInput(s, br::ObjectiveWeights::Unit);
    br::clampBlockRho(input);
    const br::ReducedLpMatrices reduced = assembleReducedLp(input);
    const std::vector<double> row_upper = br::updateReducedLpRowUpper(
        input, reduced.rows, reduced.row_upper, x_next);
    const LpResult reduced_result = solveReduced(reduced, row_upper);
    ASSERT_TRUE(reduced_result.optimal) << "s=" << s;

    // Validity does not depend on the objective, so the solution of the
    // unit-weight LP must still be feasible for the product LP. Its objective
    // value is a different number and is deliberately not compared.
    const ProductLp product = buildProductLp(input, x_next);
    const std::vector<double> lifted
        = liftToProduct(input, product, reducedZ(reduced, reduced_result));
    EXPECT_LE(maxRowViolation(product.starts, product.index, product.value,
                              product.row_upper, lifted),
              1e-7)
        << "s=" << s;
  }
}

TEST(barycentric_block_reduction, worst_residual_matches_brute_force) {
  const std::vector<double> x_next = makeXNext();

  for (const double s : {1.0, -1.0}) {
    br::ReducedLpInput input
        = makeInput(s, br::ObjectiveWeights::ProductCount);
    br::clampBlockRho(input);
    const br::ReducedLpMatrices reduced = assembleReducedLp(input);
    const std::vector<double> row_upper = br::updateReducedLpRowUpper(
        input, reduced.rows, reduced.row_upper, x_next);
    const LpResult result = solveReduced(reduced, row_upper);
    ASSERT_TRUE(result.optimal) << "s=" << s;

    const std::vector<double> z = reducedZ(reduced, result);
    const br::WorstResidual fast = br::worstReducedResidual(input, x_next, z);
    const br::WorstResidual slow = bruteForceWorst(input, x_next, z);
    EXPECT_NEAR(fast.left, slow.left, 1e-9) << "s=" << s;
    EXPECT_NEAR(fast.right, slow.right, 1e-9) << "s=" << s;
    // The scale is the size of the terms F cancels between, so it must be at
    // least as large as F itself.
    EXPECT_GE(fast.scale, std::abs(fast.worst())) << "s=" << s;
    // And the solution really is a bound.
    EXPECT_LE(fast.worst(), 1e-7) << "s=" << s;
  }
}

TEST(barycentric_block_reduction, tie_break_keeps_feasibility_and_repeats) {
  const std::vector<double> x_next = makeXNext();
  const double eps = 1e-8;

  for (const double s : {1.0, -1.0}) {
    br::ReducedLpInput plain_input
        = makeInput(s, br::ObjectiveWeights::ProductCount);
    br::clampBlockRho(plain_input);
    const br::ReducedLpMatrices plain = assembleReducedLp(plain_input);
    const LpResult plain_result = solveReduced(
        plain, br::updateReducedLpRowUpper(plain_input, plain.rows,
                                           plain.row_upper, x_next));
    ASSERT_TRUE(plain_result.optimal) << "s=" << s;

    br::ReducedLpInput input
        = makeInput(s, br::ObjectiveWeights::ProductCount, eps);
    br::clampBlockRho(input);
    const br::ReducedLpMatrices lp = assembleReducedLp(input);
    const std::vector<double> row_upper = br::updateReducedLpRowUpper(
        input, lp.rows, lp.row_upper, x_next);
    const LpResult first = solveReduced(lp, row_upper);
    const LpResult second = solveReduced(lp, row_upper);
    ASSERT_TRUE(first.optimal) << "s=" << s;
    ASSERT_TRUE(second.optimal) << "s=" << s;

    const std::vector<double> z_first = reducedZ(lp, first);
    const std::vector<double> z_second = reducedZ(lp, second);
    for (int k = 0; k < lp.cols.num_x; ++k) {
      EXPECT_EQ(z_first[static_cast<std::size_t>(k)],
                z_second[static_cast<std::size_t>(k)])
          << "s=" << s << " column " << k;
    }

    // The regularized solution is still feasible, and the unregularized part of
    // its objective is worse than the true optimum by at most O(eps * ||z||_1).
    const br::WorstResidual worst
        = br::worstReducedResidual(input, x_next, z_first);
    EXPECT_LE(worst.worst(), 1e-6) << "s=" << s;

    double unregularized = 0.0;
    for (int k = 0; k < lp.cols.num_x; ++k) {
      unregularized += lp.cost(k) * z_first[static_cast<std::size_t>(k)];
    }
    EXPECT_GE(unregularized, plain_result.objective - 1e-6) << "s=" << s;
    EXPECT_LE(unregularized,
              plain_result.objective
                  + eps * static_cast<double>(lp.cols.num_x) * kBox + 1e-6)
        << "s=" << s;
  }
}

TEST(barycentric_block_reduction, tiny_rho_is_rejected_rather_than_dropped) {
  br::ReducedLpInput input
      = makeInput(1.0, br::ObjectiveWeights::ProductCount);
  br::clampBlockRho(input);
  // A radius HiGHS would silently drop. Dropping it weakens row (L), so the
  // assembler must refuse instead.
  input.blocks[0].regions[0].vertices[0].rho(0) = 1e-12;
  EXPECT_THROW(assembleReducedLp(input), std::runtime_error);

  input.blocks[0].regions[0].vertices[0].rho(0) = -1.0;
  EXPECT_THROW(br::clampBlockRho(input), std::runtime_error);
}
