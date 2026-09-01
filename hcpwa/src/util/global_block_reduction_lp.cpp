#include "util/global_block_reduction_lp.hpp"

#include <cmath>
#include <format>
#include <limits>
#include <stdexcept>

namespace hcpwa::util::global_block_lp {

namespace {

constexpr const char* kBlockNames[kBlockCount] = {"A", "B", "C"};

void appendRow(GlobalReducedLp& lp, const SparseRow& row, double eps) {
  for (std::size_t i = 0; i < row.cols.size(); ++i) {
    if (std::abs(row.vals[i]) <= eps) {
      continue;
    }
    lp.cols.push_back(row.cols[i]);
    lp.values.push_back(row.vals[i]);
  }
  lp.starts.push_back(static_cast<int>(lp.values.size()));
}

void validate(const std::array<GlobalBlockInput, kBlockCount>& blocks,
              const GlobalReducedLpOptions& options) {
  if (!(options.dt > 0.0)) {
    throw std::invalid_argument("assembleGlobalReducedLp: dt must be positive");
  }
  if (options.sigma != 1.0 && options.sigma != -1.0) {
    throw std::invalid_argument(
        "assembleGlobalReducedLp: sigma must be +1 or -1");
  }

  std::array<int, kSpaceDim> coordinate_uses{};
  for (int b = 0; b < kBlockCount; ++b) {
    const GlobalBlockInput& block = blocks[b];
    if (block.coord_count < 2 || block.coord_count > 3) {
      throw std::invalid_argument(std::format(
          "assembleGlobalReducedLp: block {} has coord_count {}",
          kBlockNames[b], block.coord_count));
    }
    for (int d = 0; d < block.coord_count; ++d) {
      const int axis = block.coords[d];
      if (axis < 0 || axis >= kSpaceDim) {
        throw std::invalid_argument(std::format(
            "assembleGlobalReducedLp: block {} coordinate {} is {}",
            kBlockNames[b], d, axis));
      }
      if (d > 0 && axis <= block.coords[d - 1]) {
        throw std::invalid_argument(std::format(
            "assembleGlobalReducedLp: block {} coordinates are not strictly "
            "ascending",
            kBlockNames[b]));
      }
      ++coordinate_uses[axis];
    }
    if (block.rows.empty()) {
      throw std::invalid_argument(std::format(
          "assembleGlobalReducedLp: block {} has no rows", kBlockNames[b]));
    }
    for (std::size_t k = 0; k < block.rows.size(); ++k) {
      const GlobalBlockRow& row = block.rows[k];
      if (row.nu.size() != block.coord_count
          || row.p.size() != block.coord_count
          || row.r.size() != block.coord_count) {
        throw std::invalid_argument(std::format(
            "assembleGlobalReducedLp: block {} row {} has nu/p/r of size "
            "{}/{}/{}, expected {}",
            kBlockNames[b], k, row.nu.size(), row.p.size(), row.r.size(),
            block.coord_count));
      }
    }
  }

  // The three blocks must partition the eight state coordinates. If they did
  // not, p and r would not be direct sums and the rows would not separate.
  for (int axis = 0; axis < kSpaceDim; ++axis) {
    if (coordinate_uses[axis] != 1) {
      throw std::invalid_argument(std::format(
          "assembleGlobalReducedLp: coordinate {} is claimed by {} blocks, "
          "expected exactly 1",
          axis, coordinate_uses[axis]));
    }
  }
}

}  // namespace

GlobalReducedLp assembleGlobalReducedLp(
    const std::array<GlobalBlockInput, kBlockCount>& blocks,
    const GlobalReducedLpOptions& options) {
  validate(blocks, options);

  const double sigma = options.sigma;
  const double dt = options.dt;
  const double eps = options.coefficient_eps;
  const double inf = std::numeric_limits<double>::infinity();

  GlobalReducedLp lp;
  lp.dt = dt;
  lp.sigma = sigma;
  lp.num_cols = kNumCols;
  lp.objective = Eigen::RowVectorXd::Zero(kNumCols);
  lp.col_lower.assign(kNumCols, -inf);
  lp.col_upper.assign(kNumCols, inf);
  for (int b = 0; b < kBlockCount; ++b) {
    lp.coord_count[b] = blocks[b].coord_count;
    lp.coords[b] = blocks[b].coords;
    lp.num_block_rows[b] = static_cast<int>(blocks[b].rows.size());
  }

  lp.starts.push_back(0);

  for (int b = 0; b < kBlockCount; ++b) {
    const GlobalBlockInput& block = blocks[b];
    const bool carries_v = (b == kVCarryingBlock);

    for (const GlobalBlockRow& row : block.rows) {
      // sigma * (p_block^T V_block - v [block A only])
      //   + dt * r_block^T s_block - mu_block  <=  -sigma * kappa_block
      SparseRow lp_row;
      for (int d = 0; d < block.coord_count; ++d) {
        const int axis = block.coords[d];
        lp_row.add(idxV(axis), sigma * row.p(d), eps);
        // The s coefficient carries the SAME sign for both bound directions;
        // only the p and v terms and the right-hand side flip. This mirrors the
        // rho sign in the barycentric row (L) and is what keeps the modulus
        // lift valid in both directions.
        lp_row.add(idxS(axis), dt * row.r(d), eps);
      }
      if (carries_v) {
        lp_row.add(idxConst(), -sigma, eps);
      }
      lp_row.add(idxMu(b), -1.0, eps);
      appendRow(lp, lp_row, eps);
      lp.row_lower.push_back(-inf);
      lp.row_upper.push_back(-sigma * row.kappa);

      GlobalRhsTerm term;
      term.row_id = static_cast<int>(lp.row_upper.size() - 1);
      term.block = b;
      term.nu = row.nu;
      term.carries_v = carries_v;
      lp.rhs_terms.push_back(std::move(term));

      // Objective. The product form sums -sigma*p over every row of the
      // product; here each block's sum is taken once and weighted.
      for (int d = 0; d < block.coord_count; ++d) {
        lp.objective(idxV(block.coords[d]))
            += -sigma * block.objective_weight * row.p(d);
      }
      // The v coefficient is accumulated over EVERY block's rows, divided by
      // the block count, rather than over the carrying block's rows alone.
      // Both forms give sigma*R under the product weights R/R_block, so the
      // equivalence is unaffected; but the per-carrying-block form also
      // multiplies by R_block, which under production's normalised weights
      // makes the objective depend on which block happens to carry v -- and
      // block row counts span a factor of ~26. This form keeps
      // kVCarryingBlock genuinely arbitrary, as the header claims.
      lp.objective(idxConst())
          += sigma * block.objective_weight / static_cast<double>(kBlockCount);
    }
  }

  // The single coupling row: sum of the three epigraph scalars <= 0. This is
  // what replaces the R_A * R_B * R_C individual rows.
  {
    SparseRow row;
    for (int b = 0; b < kBlockCount; ++b) {
      row.add(idxMu(b), 1.0, eps);
    }
    appendRow(lp, row, eps);
    lp.row_lower.push_back(-inf);
    lp.row_upper.push_back(0.0);
  }

  // s >= +/- V, unchanged from the product form: s is global, so these 16 rows
  // are not affected by the reduction at all.
  for (int sign_index = 0; sign_index < 2; ++sign_index) {
    const double sign = sign_index == 0 ? 1.0 : -1.0;
    for (int i = 0; i < kSpaceDim; ++i) {
      SparseRow row;
      row.add(idxV(i), sign, eps);
      row.add(idxS(i), -1.0, eps);
      appendRow(lp, row, eps);
      lp.row_lower.push_back(-inf);
      lp.row_upper.push_back(0.0);
    }
  }

  // s-pinning regularisation, a constant of the formulation rather than
  // something accumulated over rows, so it is identical in both forms.
  for (int i = 0; i < kSpaceDim; ++i) {
    lp.objective(idxS(i)) += options.s_pin_weight;
  }

  return lp;
}

std::vector<double> globalReducedLpRowUpper(
    const GlobalReducedLp& lp, const std::vector<double>& v_prev) {
  if (v_prev.size() != static_cast<std::size_t>(kSpaceDim + 1)) {
    throw std::invalid_argument(std::format(
        "globalReducedLpRowUpper: v_prev has size {}, expected {}",
        v_prev.size(), kSpaceDim + 1));
  }

  std::vector<double> row_upper = lp.row_upper;
  for (const GlobalRhsTerm& term : lp.rhs_terms) {
    // b_upd is the previous step's value function at this vertex. It splits the
    // same way the rows do: the linear part restricted to the block, and the
    // constant only on the block that carries v.
    double b_upd = term.carries_v ? v_prev[kSpaceDim] : 0.0;
    for (int d = 0; d < lp.coord_count[term.block]; ++d) {
      b_upd += v_prev[lp.coords[term.block][d]] * term.nu(d);
    }
    row_upper[term.row_id] = lp.row_upper[term.row_id] - lp.sigma * b_upd;
  }
  return row_upper;
}

}  // namespace hcpwa::util::global_block_lp
