#ifndef HCPWA_UTIL_GLOBAL_BLOCK_REDUCTION_LP_HPP
#define HCPWA_UTIL_GLOBAL_BLOCK_REDUCTION_LP_HPP

#include <Eigen/Core>
#include <array>
#include <vector>

#include "util/block_reduction_lp.hpp"

// Block reduction for the GLOBAL affine approximator's main-step LP.
//
// This is a sibling of util/block_reduction_lp.hpp, not a reuse of it. The two
// LPs share the geometric fact that makes the reduction work -- the row index
// set is the full Cartesian product R_A x R_B x R_C
// (notes/step7_block_reduction.md, Lemma 1) -- but nothing else:
//
//                      barycentric                  global affine
//   value function     sum of 5 plane functions     one affine function
//   columns            eta + 3M_A+3M_B+2M_C + 6     8 + 1 + 8 + 3  (constant!)
//   row families       (R), (L), (Y)                one, plus the s rows
//   modulus lift       y per block region           one global s
//   Lemma 5 needed     yes                          no
//
// Because the modulus lift s is already a single global vector, the step that
// needed the most care in the barycentric reduction -- identifying y down to
// one vector per block region, which is only safe because rho enters row (L)
// with a plus for both bound directions -- does not arise here at all.
//
// WHAT MAKES THE ROWS SEPARABLE. One feasibility row reads
//   sigma * (p(n)^T V - v) + dt * r(n)^T s <= -sigma * kappa(n)
// with
//   p(n) = (A dt - I + dt Q_c) n + dt (f + q_c),
//   r(n) = Q_r n + q_r,
//   kappa(n) = dt * (g_n^T n + g_0).
// A is block diagonal with respect to the coordinate partition and Q_c, Q_r are
// diagonal, so the matrix acting on n is block diagonal and both p and r are
// DIRECT SUMS over the three blocks: component index rho depends only on the
// coordinates of block(rho). kappa is additively separable because g is. Hence
// p^T V = sum_block p_block^T V_block, r^T s = sum_block r_block^T s_block, and
// the row splits.
//
// The one term that does NOT split is the scalar v, which appears once per row
// regardless of block. It is assigned in full to block A (kVCarryingBlock).
// That choice is arbitrary but must be consistent between the row assembly, the
// objective and the per-step right-hand-side update, so it is named once here
// and referred to everywhere else.

namespace hcpwa::util::global_block_lp {

using hcpwa::util::block_lp::kBlockCount;
using hcpwa::util::block_lp::SparseRow;

constexpr int kSpaceDim = 8;

// Columns: [ V (8) | v (1) | s (8) | mu (3) ].
constexpr int kVOffset = 0;
constexpr int kConstOffset = kSpaceDim;
constexpr int kSOffset = kSpaceDim + 1;
constexpr int kMuOffset = 2 * kSpaceDim + 1;
constexpr int kNumCols = 2 * kSpaceDim + 1 + kBlockCount;

// The block that carries the scalar v term. See the note above.
constexpr int kVCarryingBlock = 0;

inline int idxV(int coordinate) { return kVOffset + coordinate; }
inline int idxConst() { return kConstOffset; }
inline int idxS(int coordinate) { return kSOffset + coordinate; }
inline int idxMu(int block) { return kMuOffset + block; }

// One row of one block: a vertex of that block's polytope, with the per-vertex
// data already restricted to the block's own coordinates.
struct GlobalBlockRow {
  Eigen::VectorXd nu;  // block vertex, length coord_count
  Eigen::VectorXd p;   // coeffP restricted, length coord_count
  Eigen::VectorXd r;   // radiusR restricted, length coord_count
  double kappa = 0.0;  // dt * g, this block's share
};

struct GlobalBlockInput {
  std::array<int, 3> coords{};  // state coordinates, ascending
  int coord_count = 0;          // 3 for A and B, 2 for C
  std::vector<GlobalBlockRow> rows;

  // See block_lp::BlockInput::objective_weight. R/R_block reproduces the
  // product objective exactly and is what the equivalence test uses; 1 is the
  // normalised choice taken in production.
  double objective_weight = 1.0;
};

struct GlobalReducedLpOptions {
  double dt = 1.0;
  double sigma = 1.0;  // +1 upper bound, -1 lower bound
  double s_pin_weight = 1e-6;
  double coefficient_eps = 1e-9;
};

// Per-row data needed to refresh a right-hand side once the previous step's
// value function (V_prev, v_prev) is known.
struct GlobalRhsTerm {
  int row_id = 0;
  int block = 0;
  Eigen::VectorXd nu;  // block vertex, length coord_count
  bool carries_v = false;
};

struct GlobalReducedLp {
  std::vector<int> starts;
  std::vector<int> cols;
  std::vector<double> values;
  std::vector<double> row_lower;
  std::vector<double> row_upper;  // static part, -sigma * kappa
  std::vector<double> col_lower;
  std::vector<double> col_upper;
  Eigen::RowVectorXd objective;
  std::vector<GlobalRhsTerm> rhs_terms;

  std::array<int, kBlockCount> coord_count{};
  std::array<std::array<int, 3>, kBlockCount> coords{};
  std::array<int, kBlockCount> num_block_rows{};
  int num_cols = kNumCols;
  double dt = 1.0;
  double sigma = 1.0;
};

GlobalReducedLp assembleGlobalReducedLp(
    const std::array<GlobalBlockInput, kBlockCount>& blocks,
    const GlobalReducedLpOptions& options);

// Row upper bounds given the previous step's affine value function, packed as
// v_prev = [V_prev (8); v_prev (1)] exactly as the approximator stores it.
// Cost is O(R_A + R_B + R_C).
std::vector<double> globalReducedLpRowUpper(
    const GlobalReducedLp& lp, const std::vector<double>& v_prev);

}  // namespace hcpwa::util::global_block_lp

#endif  // HCPWA_UTIL_GLOBAL_BLOCK_REDUCTION_LP_HPP
