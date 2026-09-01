#include "util/block_reduction_lp.hpp"

#include <cmath>
#include <format>
#include <limits>
#include <stdexcept>

namespace hcpwa::util::block_lp {

namespace {

constexpr const char* kBlockNames[kBlockCount] = {"A", "B", "C"};

// Appends one sparse row to the CSR arrays, skipping negligible coefficients.
void appendRow(ReducedLp& lp, const SparseRow& row, double eps) {
  for (std::size_t i = 0; i < row.cols.size(); ++i) {
    if (std::abs(row.vals[i]) <= eps) {
      continue;
    }
    lp.cols.push_back(row.cols[i]);
    lp.values.push_back(row.vals[i]);
  }
  lp.starts.push_back(static_cast<int>(lp.values.size()));
}

void validate(const std::array<BlockInput, kBlockCount>& blocks,
              const ReducedLpOptions& options) {
  if (!(options.dt > 0.0)) {
    throw std::invalid_argument("assembleReducedLp: dt must be positive");
  }
  if (options.s != 1.0 && options.s != -1.0) {
    throw std::invalid_argument("assembleReducedLp: s must be +1 or -1");
  }
  if (options.num_x <= 0) {
    throw std::invalid_argument("assembleReducedLp: num_x must be positive");
  }

  for (int b = 0; b < kBlockCount; ++b) {
    const BlockInput& block = blocks[b];
    if (block.coord_count <= 0) {
      throw std::invalid_argument(std::format(
          "assembleReducedLp: block {} has coord_count {}", kBlockNames[b],
          block.coord_count));
    }
    if (block.num_block_regions <= 0) {
      throw std::invalid_argument(std::format(
          "assembleReducedLp: block {} has no block regions", kBlockNames[b]));
    }
    if (static_cast<int>(block.psi.size()) != block.num_block_regions) {
      throw std::invalid_argument(std::format(
          "assembleReducedLp: block {} has {} Psi entries but {} block regions",
          kBlockNames[b], block.psi.size(), block.num_block_regions));
    }
    for (int j = 0; j < block.num_block_regions; ++j) {
      if (static_cast<int>(block.psi[j].rows.size()) != block.coord_count) {
        throw std::invalid_argument(std::format(
            "assembleReducedLp: block {} region {} has {} Psi rows, expected {}",
            kBlockNames[b], j, block.psi[j].rows.size(), block.coord_count));
      }
    }
    if (block.rows.empty()) {
      throw std::invalid_argument(std::format(
          "assembleReducedLp: block {} has no rows", kBlockNames[b]));
    }
    for (std::size_t r = 0; r < block.rows.size(); ++r) {
      const BlockRow& row = block.rows[r];
      if (row.block_region < 0
          || row.block_region >= block.num_block_regions) {
        throw std::invalid_argument(std::format(
            "assembleReducedLp: block {} row {} references block region {}",
            kBlockNames[b], r, row.block_region));
      }
      if (row.m.size() != block.coord_count
          || row.rho.size() != block.coord_count) {
        throw std::invalid_argument(std::format(
            "assembleReducedLp: block {} row {} has m/rho of size {}/{}, "
            "expected {}",
            kBlockNames[b], r, row.m.size(), row.rho.size(),
            block.coord_count));
      }
    }
  }
}

}  // namespace

void SparseRow::add(int col, double value, double eps) {
  if (std::abs(value) <= eps) {
    return;
  }
  for (std::size_t i = 0; i < cols.size(); ++i) {
    if (cols[i] == col) {
      vals[i] += value;
      if (std::abs(vals[i]) <= eps) {
        cols.erase(cols.begin() + static_cast<std::ptrdiff_t>(i));
        vals.erase(vals.begin() + static_cast<std::ptrdiff_t>(i));
      }
      return;
    }
  }
  cols.push_back(col);
  vals.push_back(value);
}

double SparseRow::dot(const std::vector<double>& x) const {
  double result = 0.0;
  for (std::size_t i = 0; i < cols.size(); ++i) {
    const int col = cols[i];
    if (col < 0 || col >= static_cast<int>(x.size())) {
      throw std::runtime_error("SparseRow::dot: column outside vector size");
    }
    result += vals[i] * x[col];
  }
  return result;
}

int ReducedLp::idxY(int block, int block_region, int local_dim) const {
  if (block < 0 || block >= kBlockCount) {
    throw std::invalid_argument("ReducedLp::idxY: bad block");
  }
  if (block_region < 0 || block_region >= num_block_regions[block]) {
    throw std::invalid_argument("ReducedLp::idxY: bad block region");
  }
  if (local_dim < 0 || local_dim >= coord_count[block]) {
    throw std::invalid_argument("ReducedLp::idxY: bad local dimension");
  }
  return y_offset[block] + block_region * coord_count[block] + local_dim;
}

int ReducedLp::idxMuR(int block) const {
  if (block < 0 || block >= kBlockCount) {
    throw std::invalid_argument("ReducedLp::idxMuR: bad block");
  }
  return mu_r_offset + block;
}

int ReducedLp::idxMuL(int block) const {
  if (block < 0 || block >= kBlockCount) {
    throw std::invalid_argument("ReducedLp::idxMuL: bad block");
  }
  return mu_l_offset + block;
}

ReducedLp assembleReducedLp(const std::array<BlockInput, kBlockCount>& blocks,
                            const ReducedLpOptions& options) {
  validate(blocks, options);

  const double s = options.s;
  const double dt = options.dt;
  const double eps = options.coefficient_eps;
  const double inf = std::numeric_limits<double>::infinity();

  ReducedLp lp;
  lp.num_x = options.num_x;
  lp.dt = dt;
  lp.s = s;

  int offset = options.num_x;
  for (int b = 0; b < kBlockCount; ++b) {
    lp.coord_count[b] = blocks[b].coord_count;
    lp.num_block_regions[b] = blocks[b].num_block_regions;
    lp.num_block_rows[b] = static_cast<int>(blocks[b].rows.size());
    lp.y_offset[b] = offset;
    offset += blocks[b].coord_count * blocks[b].num_block_regions;
    lp.psi[b] = blocks[b].psi;
  }
  lp.mu_r_offset = offset;
  lp.mu_l_offset = offset + kBlockCount;
  offset += 2 * kBlockCount;
  const bool tie_break = options.tie_break_epsilon > 0.0;
  if (tie_break) {
    lp.tie_break_offset = offset;
    offset += options.num_x;
  }
  lp.num_cols = offset;

  lp.objective = Eigen::RowVectorXd::Zero(lp.num_cols);
  lp.col_lower.assign(lp.num_cols, -inf);
  lp.col_upper.assign(lp.num_cols, inf);
  // y_block >= 0. The two (Y) rows give y >= +/- Psi z; this bound supplies
  // y >= 0 without extra rows.
  for (int b = 0; b < kBlockCount; ++b) {
    for (int j = 0; j < lp.num_block_regions[b]; ++j) {
      for (int p = 0; p < lp.coord_count[b]; ++p) {
        lp.col_lower[lp.idxY(b, j, p)] = 0.0;
      }
    }
  }

  lp.starts.push_back(0);

  // Two rows -- one (R), one (L) -- per block row, both carrying an RhsTerm.
  std::size_t total_block_rows = 0;
  for (int b = 0; b < kBlockCount; ++b) {
    total_block_rows += blocks[b].rows.size();
  }
  lp.rhs_terms.reserve(2 * total_block_rows);

  for (int b = 0; b < kBlockCount; ++b) {
    const BlockInput& block = blocks[b];

    for (const BlockRow& row : block.rows) {
      const BlockPsi& psi = block.psi[row.block_region];

      // rho is the radius of the uncertainty box. Small negatives can only be
      // numerical noise; larger ones mean the block did not resolve the
      // disturbance branches correctly.
      Eigen::VectorXd rho = row.rho;
      for (int p = 0; p < block.coord_count; ++p) {
        if (rho(p) < -1e-9) {
          throw std::runtime_error(std::format(
              "assembleReducedLp: block {} has negative uncertainty radius {}",
              kBlockNames[b], rho(p)));
        }
        if (rho(p) < 0.0) {
          rho(p) = 0.0;
        }
      }

      // ---- Row (R_block), step 7 section 10 ----
      //   -(s/dt) phi^T z_block - muR_block
      //       <= -s ( phi^T x_block/dt + beta_block ),
      //   beta_block = q^T m + s rho^T |q| + g,  q = Psi_{block,j} x_block.
      // The right-hand side is entirely dynamic, so the static part is 0.
      SparseRow right_row;
      for (std::size_t k = 0; k < row.phi.cols.size(); ++k) {
        right_row.add(row.phi.cols[k], -s * row.phi.vals[k] / dt, eps);
      }
      right_row.add(lp.idxMuR(b), -1.0, eps);
      appendRow(lp, right_row, eps);
      lp.row_lower.push_back(-inf);
      lp.row_upper.push_back(0.0);
      {
        RhsTerm term;
        term.row_id = static_cast<int>(lp.row_upper.size() - 1);
        term.kind = RowKind::Right;
        term.block = b;
        term.block_region = row.block_region;
        term.phi = row.phi;
        term.m = row.m;
        term.rho = rho;
        term.g = row.g;
        lp.rhs_terms.push_back(std::move(term));
      }

      // ---- Row (L_block), step 7 section 10 ----
      //   s a^T z_block + rho^T y_{block,j} - muL_block <= -s b,
      //   a = Psi^T m - phi/dt,   b = phi^T x_block/dt + g.
      //
      // rho enters with a PLUS for both signs s; only a and the right-hand
      // side flip. Lemma 5 of step 7 section 8 -- the identification of y down
      // to one vector per block region -- rests on exactly that fact, because
      // it is what makes lowering y only loosen this row. If an s ever appears
      // on rho here, block-level y becomes invalid and Lemma 5 must be
      // re-derived.
      SparseRow a_row;
      for (int p = 0; p < block.coord_count; ++p) {
        const SparseRow& psi_p = psi.rows[p];
        for (std::size_t k = 0; k < psi_p.cols.size(); ++k) {
          a_row.add(psi_p.cols[k], psi_p.vals[k] * row.m(p), eps);
        }
      }
      for (std::size_t k = 0; k < row.phi.cols.size(); ++k) {
        a_row.add(row.phi.cols[k], -row.phi.vals[k] / dt, eps);
      }

      SparseRow left_row;
      for (std::size_t k = 0; k < a_row.cols.size(); ++k) {
        left_row.add(a_row.cols[k], s * a_row.vals[k], eps);
      }
      for (int p = 0; p < block.coord_count; ++p) {
        left_row.add(lp.idxY(b, row.block_region, p), rho(p), eps);
      }
      left_row.add(lp.idxMuL(b), -1.0, eps);
      appendRow(lp, left_row, eps);
      lp.row_lower.push_back(-inf);
      // Static part of -s b; the -s phi^T x_next / dt part is added per step.
      lp.row_upper.push_back(-s * row.g);
      {
        RhsTerm term;
        term.row_id = static_cast<int>(lp.row_upper.size() - 1);
        term.kind = RowKind::Left;
        term.block = b;
        term.block_region = row.block_region;
        term.phi = row.phi;
        lp.rhs_terms.push_back(std::move(term));
      }

      // ---- Objective, step 7 section 9 ----
      //   delta^R = -s F^R = (s/dt) phi^T z + const,  minimised.
      // The weight is the caller's decision; see BlockInput::objective_weight.
      for (std::size_t k = 0; k < row.phi.cols.size(); ++k) {
        lp.objective(row.phi.cols[k])
            += block.objective_weight * (s / dt) * row.phi.vals[k];
      }
    }

    // ---- Rows (Y_block): y_{block,j} >= +/- Psi_{block,j} z_block ----
    for (int j = 0; j < block.num_block_regions; ++j) {
      const BlockPsi& psi = block.psi[j];
      for (int sign_index = 0; sign_index < 2; ++sign_index) {
        const double sign = sign_index == 0 ? 1.0 : -1.0;
        for (int p = 0; p < block.coord_count; ++p) {
          SparseRow row;
          const SparseRow& psi_p = psi.rows[p];
          for (std::size_t k = 0; k < psi_p.cols.size(); ++k) {
            row.add(psi_p.cols[k], sign * psi_p.vals[k], eps);
          }
          row.add(lp.idxY(b, j, p), -1.0, eps);
          appendRow(lp, row, eps);
          lp.row_lower.push_back(-inf);
          lp.row_upper.push_back(0.0);
        }
      }
    }
  }

  // ---- Tie-break rows: t_i >= |z_i| ----
  //
  // With epsilon * sum(t) in the objective these pick a specific point out of
  // the optimal face instead of whichever vertex the solver happened to reach.
  // The optimal face here is genuinely wider than a point (see
  // ReducedLpOptions::tie_break_epsilon), so without this two runs diverge at
  // the first tie and every later step of the march differs.
  if (tie_break) {
    for (int i = 0; i < options.num_x; ++i) {
      const int t_col = lp.tie_break_offset + i;
      lp.col_lower[t_col] = 0.0;
      lp.objective(t_col) += options.tie_break_epsilon;
      for (int sign_index = 0; sign_index < 2; ++sign_index) {
        const double sign = sign_index == 0 ? 1.0 : -1.0;
        SparseRow row;
        row.add(i, sign, eps);
        row.add(t_col, -1.0, eps);
        appendRow(lp, row, eps);
        lp.row_lower.push_back(-inf);
        lp.row_upper.push_back(0.0);
      }
    }
  }

  // ---- Rows (Sigma): the only coupling between the three blocks ----
  // sum_block muR_block <= 0 and sum_block muL_block <= 0. Each row involves
  // one family only, so the two families never mix scales (step 7 section 10).
  for (int family = 0; family < 2; ++family) {
    SparseRow row;
    for (int b = 0; b < kBlockCount; ++b) {
      row.add(family == 0 ? lp.idxMuR(b) : lp.idxMuL(b), 1.0, eps);
    }
    appendRow(lp, row, eps);
    lp.row_lower.push_back(-inf);
    lp.row_upper.push_back(0.0);
  }

  return lp;
}

std::vector<double> reducedLpRowUpper(const ReducedLp& lp,
                                      const std::vector<double>& x_next) {
  if (x_next.size() != static_cast<std::size_t>(lp.num_x)) {
    throw std::invalid_argument(std::format(
        "reducedLpRowUpper: x_next has size {}, expected {}", x_next.size(),
        lp.num_x));
  }

  const double s = lp.s;
  const double dt = lp.dt;

  // q_{block,j} = Psi_{block,j} x_block, computed once per block region and
  // shared by every Right row of that block region.
  std::array<std::vector<Eigen::VectorXd>, kBlockCount> q;
  for (int b = 0; b < kBlockCount; ++b) {
    q[b].resize(lp.psi[b].size());
    for (std::size_t j = 0; j < lp.psi[b].size(); ++j) {
      q[b][j] = Eigen::VectorXd::Zero(lp.coord_count[b]);
      for (int p = 0; p < lp.coord_count[b]; ++p) {
        q[b][j](p) = lp.psi[b][j].rows[p].dot(x_next);
      }
    }
  }

  std::vector<double> row_upper = lp.row_upper;
  for (const RhsTerm& term : lp.rhs_terms) {
    const double phi_x = term.phi.dot(x_next);
    if (term.kind == RowKind::Left) {
      row_upper[term.row_id] = lp.row_upper[term.row_id] - s * phi_x / dt;
    } else {
      const Eigen::VectorXd& qj = q[term.block][term.block_region];
      const double beta
          = qj.dot(term.m) + s * term.rho.dot(qj.cwiseAbs()) + term.g;
      row_upper[term.row_id] = -s * (phi_x / dt + beta);
    }
  }
  return row_upper;
}

}  // namespace hcpwa::util::block_lp
