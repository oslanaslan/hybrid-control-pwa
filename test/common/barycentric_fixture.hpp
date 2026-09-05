#ifndef HCPWA_TEST_BARYCENTRIC_FIXTURE_HPP
#define HCPWA_TEST_BARYCENTRIC_FIXTURE_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <random>
#include <vector>

#include <Eigen/Core>
#include <Highs.h>

#include <barycentric_geometry_types.hpp>
#include <util/block_reduction_lp.hpp>

// The synthetic block-reduced LP input the assembler tests run on, and a bare
// HiGHS runner for it. Shared by barycentric_block_reduction.cpp (one-step LP
// against its product form) and barycentric_band_lp.cpp (band LP against the
// one-step LP), so that both argue about the same object.
namespace barycentric_test_fixture {

using barycentric_affine_approximator::SparseVec;
namespace br = barycentric_affine_approximator::block_reduction;

constexpr double kDt = 0.25;
constexpr double kBox = 5.0;
constexpr int kSeed = 7;

// Shape of the fixture. Sizes are the ones the reduction arithmetic was checked
// against by hand: 936 rows / 159 columns in product form, 90 / 43 reduced.
inline const std::vector<int> kCoords = {3, 3, 2};
inline const std::vector<int> kEta = {6, 5, 4};
inline const std::vector<int> kNumRegions = {3, 3, 2};
inline const std::vector<int> kNumVertices = {3, 3, 2};
// Number of projection layers each block owns; only used to give phi a
// realistic total mass.
inline const std::vector<int> kLayerCount = {2, 2, 1};

inline int totalNumX() { return std::accumulate(kEta.begin(), kEta.end(), 0); }

inline int blockColBegin(int block) {
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
inline br::ReducedLpInput makeInput(double sign_s, br::ObjectiveWeights weights,
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

inline std::vector<double> makeXNext() {
  std::mt19937 rng(kSeed + 101);
  std::uniform_real_distribution<double> sym(-1.0, 1.0);
  std::vector<double> x(static_cast<std::size_t>(totalNumX()));
  for (double& value : x) {
    value = sym(rng);
  }
  return x;
}

struct LpResult {
  bool optimal = false;
  HighsModelStatus status = HighsModelStatus::kNotset;
  double objective = 0.0;
  std::vector<double> col_value;
};

inline LpResult solveLp(const std::vector<int>& starts, const std::vector<int>& index,
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
  result.status = highs.getModelStatus();
  result.optimal = result.status == HighsModelStatus::kOptimal;
  result.objective = highs.getInfo().objective_function_value;
  result.col_value = highs.getSolution().col_value;
  return result;
}

}  // namespace barycentric_test_fixture

#endif  // HCPWA_TEST_BARYCENTRIC_FIXTURE_HPP
