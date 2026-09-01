#ifndef HCPWA_UTIL_BLOCK_REDUCTION_LP_HPP
#define HCPWA_UTIL_BLOCK_REDUCTION_LP_HPP

#include <Eigen/Core>
#include <array>
#include <vector>

// Assembler for the block-reduced main-step LP of
//   notes/step7_block_reduction.md, section 10
// (referred to below as "step 7"), and
//   docs/barycentric_block_reduction_context.md, part II.4
// (referred to below as "ctx").
//
// This module is deliberately independent of BarycentricAffineApproximator. It
// is driven by the Tier-1 equivalence test on synthetic data
// (test/common/barycentric_block_reduction.cpp) before any production caller
// exists, so that when prepareLpMatrices() is switched over to it, the code it
// switches to is already under test. That ordering is the whole point of
// ctx part IX step 1: the production product LP is not solvable at production
// scale, so a synthetic equivalence test is the only correctness oracle there
// is.
//
// SCALE CONVENTION. Both row families are written at the plain F scale of
// step 7 section 10, i.e. row (R) here is the (R) of step 2.1 divided by dt.
// The pre-existing approximator writes (R) at the dt*F scale and (L) at the
// plain F scale; that mixture is legal (mu^R and mu^L never appear in the same
// row, see step 7 section 10) but it makes the code hard to read against the
// derivation. What is NOT legal is scaling only some rows of one family, since
// then mu_block is compared against inconsistent units. One scale per family.
//
// COLUMN INDICES ARE ALWAYS GLOBAL. phi and Psi are supplied in the caller's
// global x column space. The three blocks have disjoint column supports
// (step 7 section 2.4) but this module never relies on that: it simply emits
// whatever columns it is given. Keeping rows blockwise while columns stay
// global is the discipline that ctx part VIII flags as the second-highest risk
// in this work.

namespace hcpwa::util::block_lp {

// A, B, C in phase 0; A', B', C' in phase 1.
constexpr int kBlockCount = 3;

// Sparse row over global columns.
//
// This mirrors barycentric_affine_approximator::SparseVec. It is duplicated
// rather than shared because step 1 of the work order must not modify
// production code; once prepareLpMatrices() is ported (step 6), SparseVec
// should become an alias of this type and the duplicate should go away.
struct SparseRow {
  std::vector<int> cols;
  std::vector<double> vals;

  // Adds value at col, merging duplicates. Entries whose magnitude falls to
  // eps or below are dropped, including entries that cancel on merge.
  void add(int col, double value, double eps = 1e-12);
  double dot(const std::vector<double>& x) const;
};

// Psi_{block, j_block}: one row per block coordinate, over global columns.
// Rows are indexed by the block's coordinates in increasing order of the state
// coordinate number (step 7 section 2.5).
struct BlockPsi {
  std::vector<SparseRow> rows;
};

// One row a = (j_block, nu_block) of one block's row set R_block.
struct BlockRow {
  int block_region = 0;  // j_block, indexes BlockInput::psi
  SparseRow phi;         // phi_block(a), global columns
  Eigen::VectorXd m;     // m_block(a), length coord_count
  Eigen::VectorXd rho;   // rho_block(a) >= 0, length coord_count
  double g = 0.0;        // g_block(a)
};

// Everything the reduction needs about one block.
struct BlockInput {
  int coord_count = 0;       // |I_block|: 3 for A and B, 2 for C
  int num_block_regions = 0;  // M_block
  std::vector<BlockPsi> psi;  // size num_block_regions
  std::vector<BlockRow> rows;  // size R_block

  // Weight on this block's contribution to the objective.
  //
  // Two settings matter and they are NOT interchangeable:
  //   R / R_block : reproduces the product objective exactly, hence the same
  //                 optimal value AND the same optimal z (step 7 section 9).
  //                 This is the only setting under which the reduced LP is
  //                 equivalent to the product LP, so it is what the Tier-1
  //                 equivalence assertion must use.
  //   1           : normalised. A different problem with a different z*, which
  //                 ctx VI.1 adopts deliberately for production because the
  //                 product objective silently weights each block by the row
  //                 count of the other two (a 26x spread) and was never
  //                 computable at production scale anyway. Validity of the
  //                 bound does not depend on the objective -- any feasible
  //                 point is a certificate -- only tightness does.
  double objective_weight = 1.0;
};

struct ReducedLpOptions {
  double dt = 1.0;
  double s = 1.0;   // +1 upper bound, -1 lower bound (step 2.1)
  int num_x = 0;    // eta: number of global x columns
  double coefficient_eps = 1e-12;

  // Deterministic tie-break. Adds epsilon * ||z||_1 to the objective, via one
  // auxiliary column t_i >= |z_i| per x column and two rows apiece.
  //
  // This is not optional polish. The optimal face of this LP is genuinely more
  // than a point: in a run over 120 synthetic instances one case showed two
  // optima with ||z1 - z2||_inf = 7.66 at objective values agreeing to 5e-13,
  // and the real objective is structured (c_z proportional to Phi^T 1), so
  // degeneracy is expected to be MORE common there, not less. Without a rule
  // for choosing among optima, two runs march differently from the first tied
  // step onward and stop being comparable.
  //
  // The cost is that the returned point is optimal for the perturbed
  // objective, so the bound is looser by O(epsilon). It stays a valid
  // certificate: validity needs only feasibility.
  //
  // Zero disables it, which is what the equivalence test uses so that its row
  // and column counts stay comparable to the product form.
  double tie_break_epsilon = 0.0;
};

// Which endpoint of [t_{k-1}, t_k] a residual row controls. Left is the
// endpoint with the unknown z, so its modulus is lifted onto y; Right is the
// endpoint with the known x_next, so its modulus is arithmetic.
enum class RowKind { Left, Right };

// Per-row data needed to recompute a right-hand side after x_next changes.
struct RhsTerm {
  int row_id = 0;
  RowKind kind = RowKind::Left;
  int block = 0;
  int block_region = 0;  // Right only: selects q[block][block_region]
  SparseRow phi;
  Eigen::VectorXd m;    // Right only
  Eigen::VectorXd rho;  // Right only
  double g = 0.0;       // Right only
};

// Variable layout:
//   [ x (num_x) | y_A | y_B | y_C | muR (3) | muL (3) ]
// with y_block holding coord_count entries per block region.
struct ReducedLp {
  // Row-major sparse matrix, HiGHS addRows format.
  std::vector<int> starts;
  std::vector<int> cols;
  std::vector<double> values;
  std::vector<double> row_lower;
  // Static part of the upper bound. Right rows are fully dynamic and sit at 0
  // here; Left rows carry their -s*g term. Use reducedLpRowUpper() to obtain
  // the bounds for a given x_next.
  std::vector<double> row_upper;

  // x is left free: gauge fixing and any artificial box are the caller's
  // business. y >= 0, mu free.
  std::vector<double> col_lower;
  std::vector<double> col_upper;

  Eigen::RowVectorXd objective;
  std::vector<RhsTerm> rhs_terms;

  int num_cols = 0;
  int num_x = 0;
  std::array<int, kBlockCount> y_offset{};
  std::array<int, kBlockCount> coord_count{};
  std::array<int, kBlockCount> num_block_regions{};
  std::array<int, kBlockCount> num_block_rows{};
  int mu_r_offset = 0;
  int mu_l_offset = 0;
  // Start of the tie-break auxiliary columns, or -1 when disabled. They sit
  // after mu so that adding them leaves every other offset unchanged.
  int tie_break_offset = -1;

  // Psi copied in so that per-step q = Psi x can be recomputed without the
  // caller having to keep BlockInput alive.
  std::array<std::vector<BlockPsi>, kBlockCount> psi;
  double dt = 1.0;
  double s = 1.0;

  int idxY(int block, int block_region, int local_dim) const;
  int idxMuR(int block) const;
  int idxMuL(int block) const;
};

// Builds the reduced LP. Throws std::invalid_argument on shape mismatches.
ReducedLp assembleReducedLp(const std::array<BlockInput, kBlockCount>& blocks,
                            const ReducedLpOptions& options);

// Row upper bounds for a given right-endpoint value x_next. Cost is
// O(R_A + R_B + R_C) plus O(M_A + M_B + M_C) for the q vectors -- this is the
// only place the modulus is arithmetic (step 7 section 13).
std::vector<double> reducedLpRowUpper(const ReducedLp& lp,
                                      const std::vector<double>& x_next);

}  // namespace hcpwa::util::block_lp

#endif  // HCPWA_UTIL_BLOCK_REDUCTION_LP_HPP
