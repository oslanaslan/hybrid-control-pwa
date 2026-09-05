#ifndef HCPWA_UTIL_BAND_LP_HPP
#define HCPWA_UTIL_BAND_LP_HPP

#include <Eigen/Core>
#include <vector>

#include "util/block_reduction_lp.hpp"

// The band LP of one node (r, i, theta): the block-reduced residual rows on
// every segment [t_k, t_{k+1}] of the time grid t_0 < ... < t_L = theta at
// once, with z_0 .. z_{L-1} unknown and z_L given (paper, section "Ленточная
// задача линейного программирования", problem lp_horizon).
//
// Node values are linear in time on each segment, so on segment k the slope is
// dz_k = (z_{k+1} - z_k) / dt and the residual is checked at both ends of the
// segment with that slope (lemma on interval endpoints): at t_k with the
// gradient of z_k, at t_{k+1} with the gradient of z_{k+1}. Stage k therefore
// carries the same rows as the one-step LP of block_reduction_lp.hpp,
//
//   (L_k)  s a^T z_k + (s/dt) phi^T z_{k+1} + rho^T y_k - muL_k <= -s g
//   (R_k) -(s/dt) phi^T z_k + s (Psi^T m + phi/dt)^T z_{k+1} + rho^T y_{k+1}
//                                                        - muR_k <= -s g
//   (Y_k)  +-Psi z_k - y_k <= 0,  sum_b muR_{k,b} <= 0,  sum_b muL_{k,b} <= 0
//
// and at the last stage z_L is data and its terms move to the right-hand side,
// where they read exactly as assembleReducedLp + updateReducedLpRowUpper. The
// one-step LP is the L = 1 case of this one, and is kept as its reference.
//
// The objective is the integral of the estimate (lemma on quadrature):
//   min  s * sum_{k<L} omega_k q^T z_k,  omega_0 = 1/2, omega_k = 1,
// q being the node weights; the z_L term and the dt factor are constants.
//
// Naming style follows barycentric_affine_approximator.
// NOLINTBEGIN(readability-identifier-naming)

namespace barycentric_affine_approximator {
namespace block_reduction {

struct BandLpInput {
  // Block data of the phase, shared by every stage. Must have been through
  // clampBlockRho, and must not ask for a tie-break: the integral objective is
  // strictly increasing in every node value, so there is nothing to break.
  const ReducedLpInput* stage = nullptr;
  // L >= 1.
  int num_stages = 0;
  // z_L, num_x entries.
  std::vector<double> z_terminal;
  // q, num_x entries: the integral of each node's hat function over Omega.
  Eigen::VectorXd node_weights;
  // x columns fixed to zero at every stage k < L: the gauge normalisation and,
  // at level r = 0, the whole plane of the third block. z_terminal must obey
  // the same pins; the caller checks that, this file only applies them.
  std::vector<int> pinned_columns;
};

// Columns: [ stage 0 | stage 1 | ... | stage L-1 ], each stage laid out as the
// one-step ReducedLpColLayout without its tie-break block.
struct BandLpColLayout {
  ReducedLpColLayout stage;
  int stage_cols = 0;
  int num_stages = 0;
  int num_cols = 0;

  int idxX(int k, int column) const;
  int idxYBlock(int k, int block, int block_region, int p) const;
  int idxMuR(int k, int block) const;
  int idxMuL(int k, int block) const;
};

// Rows: [ stage 0 | stage 1 | ... | stage L-1 ], each stage in the order of
// ReducedLpRowLayout: (L,R) pairs per block, then the |.| rows, then the two
// coupling rows. Row ids are closed formulas, as before.
struct BandLpRowLayout {
  ReducedLpRowLayout stage;
  int stage_rows = 0;
  int num_stages = 0;
  int num_rows = 0;

  int rowLeft(int k, int block, int block_region, int vertex) const;
  int rowRight(int k, int block, int block_region, int vertex) const;
  int rowAbs(int k, int block, int block_region, int p, bool positive) const;
  int rowSumR(int k) const;
  int rowSumL(int k) const;
};

struct BandLpMatrices {
  BandLpColLayout cols;
  BandLpRowLayout rows;
  // CSR of the constraint matrix.
  std::vector<int> starts;
  std::vector<int> col_index;
  std::vector<double> value;
  std::vector<double> row_lower;
  std::vector<double> row_upper;
  // Unnormalised: s * omega_k * q on the x columns of stage k, zero elsewhere.
  // Divide by cost_scale before handing it to a solver -- the weights carry
  // N^6 -- and multiply the reported objective back by it.
  Eigen::RowVectorXd cost;
  double cost_scale = 1.0;
  // y >= 0, pinned x columns fixed at zero, everything else free.
  std::vector<double> col_lower;
  std::vector<double> col_upper;
};

// omega_k of the trapezoid rule on L segments: 1/2 at both ends, 1 inside.
double trapezoidWeight(int k, int num_stages);

// Builds the band LP. Throws on a malformed input, on a tie-break request, on
// a rho that clampBlockRho did not normalise, or on an all-zero objective.
BandLpMatrices assembleBandLp(const BandLpInput& input);

// dt * sum_{k=0}^{L} omega_k q^T z_k for z = (z_0, ..., z_L): the integral of
// the estimate over [t_0, theta] x Omega. Zero when there is no segment.
double bandIntegral(const Eigen::VectorXd& node_weights, double t_delta,
                    const std::vector<std::vector<double>>& z);

// The exact worst s * F of every stage, stage k being worstReducedResidual with
// z_k on the left and z_{k+1} on the right. z has L + 1 entries.
std::vector<WorstResidual> bandStageResiduals(
    const ReducedLpInput& stage, const std::vector<std::vector<double>>& z);

}  // namespace block_reduction
}  // namespace barycentric_affine_approximator

// NOLINTEND(readability-identifier-naming)

#endif  // HCPWA_UTIL_BAND_LP_HPP
