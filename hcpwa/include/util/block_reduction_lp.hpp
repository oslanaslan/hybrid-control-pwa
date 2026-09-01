#ifndef HCPWA_UTIL_BLOCK_REDUCTION_LP_HPP
#define HCPWA_UTIL_BLOCK_REDUCTION_LP_HPP

#include <Eigen/Core>
#include <vector>

#include "barycentric_geometry_types.hpp"

// Block reduction of the barycentric residual LP.
//
// The full 8D region of one phase is a direct product of three coordinate
// blocks, Omega^(j) = P_A x P_B x P_C, and the residual is additively separable
// over those blocks:
//
//   F(j, nu) = F_A(j_A, nu_A) + F_B(j_B, nu_B) + F_C(j_C, nu_C).
//
// Enumerating the product costs prod_b R_b rows. Lifting the three maxima into
// epigraph variables mu_b and coupling them by sum_b mu_b <= 0 costs
// sum_b R_b rows instead, and is an *equality* of the projections onto z, not a
// relaxation: see the derivation in the .cpp next to assembleReducedLp().
//
// This file owns the row builder, the per-step RHS update and the exact worst
// residual together on purpose. All three have to agree on the row scale, and
// they used to live 350 lines apart.
//
// Naming style follows barycentric_affine_approximator.
// NOLINTBEGIN(readability-identifier-naming)

namespace barycentric_affine_approximator {
namespace block_reduction {

// Pruning threshold used while a row is being built. It is far below the HiGHS
// small_matrix_value so that nothing is dropped here that HiGHS would have
// kept; the point of pruning at all is to keep the sparse rows short.
constexpr double kAssembleEps = 1e-12;

// Mirrors the HiGHS "small_matrix_value" option. HiGHS silently drops matrix
// entries at or below it, which for a rho coefficient would *weaken* row (L)
// and destroy the bound property, so emitted rho entries are checked against
// this value.
constexpr double kSmallMatrixValue = 1e-9;

// Residual row upper bounds this close to zero are snapped to zero, matching
// the behaviour the product assembler has always had.
constexpr double kRhsSnapEps = kEps;

// One vertex of one block-region: everything the two residual rows need.
struct BlockVertexData {
  // Value selector over the *global* x columns: V_b(nu_b) = phi^T x.
  SparseVec phi;
  // m_b = A_b nu + f_b + c_b(nu), indexed by the block's own coordinates.
  Eigen::VectorXd m;
  // Disturbance box radius at nu, indexed by the block's own coordinates.
  // Must be nonnegative: it enters row (L) with a plus sign for both
  // approximation directions, and that sign is what makes the y-lift valid.
  Eigen::VectorXd rho;
  // The block's share of the affine term g_i(nu).
  double g = 0.0;
};

// One block-region j_b. Psi_{b,j} depends on j_b only, which is exactly why the
// auxiliary y can be indexed per block-region instead of per product region.
struct BlockRegionData {
  // coord_count rows, each a sparse vector over the global x columns.
  std::vector<SparseVec> psi_rows;
  std::vector<BlockVertexData> vertices;
};

struct ReducedLpBlock {
  int coord_count = 0;
  std::vector<BlockRegionData> regions;
};

// Weights of the per-block objective terms.
enum class ObjectiveWeights {
  // w_b = R / R_b with R = prod_b R_b. Reproduces the product objective
  // coefficient for coefficient, and keeps the LP bounded: under a gauge shift
  // of block b the objective moves by -(1/dt) * w_b * R_b * theta_b, and with
  // these weights w_b * R_b = R for *every* block, so the shift is bounded by
  // the same feasibility row that bounds the shift itself.
  ProductCount,
  // w_b = 1. Better conditioned, but w_b * R_b then differs across blocks and
  // raising one block's gauge while lowering another's drives the objective
  // down along a feasible ray. Kept as an option; treat a kUnbounded model
  // status under it as this, not as a solver hiccup.
  Unit,
};

struct ReducedLpInput {
  int num_x = 0;
  double t_delta = 0.0;
  // s of step 2.1: +1 for the upper bound, -1 for the lower one.
  double sign_s = 1.0;
  ObjectiveWeights weights = ObjectiveWeights::ProductCount;
  // Weak epsilon * ||z||_1 regularization for a reproducible tie-break. The u
  // block and its rows are emitted only when this is positive.
  double tie_break_eps = 0.0;
  std::vector<ReducedLpBlock> blocks;
};

// Column layout:
//   [ x (num_x) | y_A | y_B | y_C | muR (B) | muL (B) | u (num_x) ]
// The u block sits last so the [x | y | mu] prefix does not depend on whether
// the tie-break is enabled.
struct ReducedLpColLayout {
  int num_x = 0;
  std::vector<int> coord_count;
  std::vector<int> num_regions;
  // Absolute column of y_{b,0,0}.
  std::vector<int> y_offset;
  int mu_r_offset = 0;
  int mu_l_offset = 0;
  // -1 when the tie-break block is absent.
  int u_offset = -1;
  int num_cols = 0;

  int numBlocks() const { return static_cast<int>(coord_count.size()); }
  int idxX(int column) const;
  int idxYBlock(int block, int block_region, int p) const;
  int idxMuR(int block) const;
  int idxMuL(int block) const;
  int idxU(int k) const;
};

// Row layout. Every row id is a closed formula of (block, region, vertex) so
// that no caller has to reconstruct an id from the emission order.
//   [ (L,R) pairs per block | abs rows per block | 2 coupling | tie-break ]
struct ReducedLpRowLayout {
  std::vector<int> coord_count;
  std::vector<int> num_regions;
  // sum of R_{b'} over b' < b, where R_b counts (region, vertex) pairs.
  std::vector<int> pair_offset;
  // per block, per region: pairs of that block before that region.
  std::vector<std::vector<int>> region_pair_offset;
  // sum of coord_count_{b'} * num_regions_{b'} over b' < b.
  std::vector<int> y_offset_local;
  int total_pairs = 0;
  int abs_offset = 0;
  int coupling_offset = 0;
  int tie_break_offset = -1;
  int num_x = 0;
  int num_rows = 0;

  int numBlocks() const { return static_cast<int>(coord_count.size()); }
  int rowLeft(int block, int block_region, int vertex) const;
  int rowRight(int block, int block_region, int vertex) const;
  int rowAbs(int block, int block_region, int p, bool positive) const;
  int rowSumR() const { return coupling_offset; }
  int rowSumL() const { return coupling_offset + 1; }
  int rowTieBreak(int k, bool positive) const;
};

struct ReducedLpMatrices {
  ReducedLpColLayout cols;
  ReducedLpRowLayout rows;
  // CSR of the constraint matrix.
  std::vector<int> starts;
  std::vector<int> col_index;
  std::vector<double> value;
  std::vector<double> row_lower;
  // The x_next-independent part of the upper bounds. Feed it to
  // updateReducedLpRowUpper() as base_upper.
  std::vector<double> row_upper;
  Eigen::RowVectorXd cost;
  // y >= 0 and u >= 0 only. Gauge fixing of the x block is the caller's job.
  std::vector<double> col_lower;
  std::vector<double> col_upper;
};

struct RhoClampStats {
  // Small negatives, which are numerical noise around an exact zero radius.
  int negatives_clamped = 0;
  // Magnitudes at or below kSmallMatrixValue, snapped to exactly zero. rho is
  // affine in nu with coefficients of order 1e-1 to 1e1, so a magnitude that
  // small means the vertex sits on the branch line where the two bounds
  // coincide and the radius is exactly zero. Snapping is what keeps HiGHS from
  // dropping the entry silently instead.
  int tiny_snapped = 0;
};

// Normalises the disturbance radii. Call this once on the input: the row
// builder, the RHS update and the residual check all read rho exactly as
// given, so this has to happen in one place. Throws on a radius that is
// negative beyond numerical noise.
RhoClampStats clampBlockRho(ReducedLpInput& input, double tolerance = kEps);

ReducedLpColLayout makeReducedLpColLayout(const ReducedLpInput& input);
ReducedLpRowLayout makeReducedLpRowLayout(const ReducedLpInput& input);

// Builds the static part of the reduced LP. Throws if rho is negative beyond
// numerical noise, or if a nonzero rho entry would be dropped by HiGHS.
ReducedLpMatrices assembleReducedLp(const ReducedLpInput& input);

// Recomputes the residual row upper bounds for a new x_next. base_upper is the
// row_upper returned by assembleReducedLp().
std::vector<double> updateReducedLpRowUpper(
    const ReducedLpInput& input, const ReducedLpRowLayout& rows,
    const std::vector<double>& base_upper, const std::vector<double>& x_next);

struct WorstResidual {
  double left = 0.0;
  double right = 0.0;
  double worst() const { return left > right ? left : right; }
};

// Exact worst s * F over the whole product of regions and vertices, computed
// as sum_b max_{a in R_b} s * F_b(a). Separability makes this the true global
// maximum, not a sample of it, at O(sum_b R_b) cost.
WorstResidual worstReducedResidual(const ReducedLpInput& input,
                                   const std::vector<double>& x_next,
                                   const std::vector<double>& z);

}  // namespace block_reduction
}  // namespace barycentric_affine_approximator

// NOLINTEND(readability-identifier-naming)

#endif  // HCPWA_UTIL_BLOCK_REDUCTION_LP_HPP
