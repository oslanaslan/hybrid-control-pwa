#include "util/block_reduction_lp.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>

// NOLINTBEGIN(readability-identifier-naming)

namespace barycentric_affine_approximator {
namespace block_reduction {

namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();

void checkBlocks(const ReducedLpInput& input) {
  if (input.blocks.empty()) {
    throw std::invalid_argument("block_reduction: no blocks");
  }
  if (input.num_x <= 0) {
    throw std::invalid_argument("block_reduction: num_x must be positive");
  }
  if (!(input.t_delta > 0.0)) {
    throw std::invalid_argument("block_reduction: t_delta must be positive");
  }
  if (input.sign_s != 1.0 && input.sign_s != -1.0) {
    throw std::invalid_argument("block_reduction: sign_s must be +-1");
  }
  for (std::size_t b = 0; b < input.blocks.size(); ++b) {
    const ReducedLpBlock& block = input.blocks[b];
    if (block.coord_count <= 0) {
      throw std::invalid_argument("block_reduction: empty coordinate block "
                                  + std::to_string(b));
    }
    if (block.regions.empty()) {
      throw std::invalid_argument("block_reduction: block "
                                  + std::to_string(b) + " has no regions");
    }
    for (std::size_t j = 0; j < block.regions.size(); ++j) {
      const BlockRegionData& region = block.regions[j];
      if (static_cast<int>(region.psi_rows.size()) != block.coord_count) {
        throw std::invalid_argument(
            "block_reduction: Psi row count does not match coord_count in "
            "block " + std::to_string(b) + " region " + std::to_string(j));
      }
      if (region.vertices.empty()) {
        throw std::invalid_argument(
            "block_reduction: block " + std::to_string(b) + " region "
            + std::to_string(j) + " has no vertices");
      }
      for (const BlockVertexData& vertex : region.vertices) {
        if (vertex.m.size() != block.coord_count
            || vertex.rho.size() != block.coord_count) {
          throw std::invalid_argument(
              "block_reduction: m/rho size does not match coord_count in block "
              + std::to_string(b));
        }
      }
    }
  }
}

// Psi_{b,j} x for one block-region.
Eigen::VectorXd psiTimes(const BlockRegionData& region,
                         const std::vector<double>& x) {
  const int rows = static_cast<int>(region.psi_rows.size());
  Eigen::VectorXd out(rows);
  for (int p = 0; p < rows; ++p) {
    out(p) = region.psi_rows[p].dot(x);
  }
  return out;
}

// Number of (region, vertex) pairs of one block.
int blockPairCount(const ReducedLpBlock& block) {
  int total = 0;
  for (const BlockRegionData& region : block.regions) {
    total += static_cast<int>(region.vertices.size());
  }
  return total;
}

// Appends one sparse row to the CSR arrays. Pruning has already happened while
// the row was assembled, with a per-quantity epsilon, so nothing is dropped
// here: dropping a rho coefficient at emission time would silently weaken the
// constraint.
void appendRow(std::vector<int>& starts, std::vector<int>& cols,
               std::vector<double>& values, const SparseVec& row) {
  for (std::size_t i = 0; i < row.cols.size(); ++i) {
    cols.push_back(row.cols[i]);
    values.push_back(row.vals[i]);
  }
  starts.push_back(static_cast<int>(values.size()));
}

}  // namespace

int ReducedLpColLayout::idxX(int column) const {
  if (column < 0 || column >= num_x) {
    throw std::invalid_argument("ReducedLpColLayout::idxX: bad column");
  }
  return column;
}

int ReducedLpColLayout::idxYBlock(int block, int block_region, int p) const {
  if (block < 0 || block >= numBlocks()) {
    throw std::invalid_argument("ReducedLpColLayout::idxYBlock: bad block");
  }
  if (block_region < 0 || block_region >= num_regions[block]) {
    throw std::invalid_argument("ReducedLpColLayout::idxYBlock: bad region");
  }
  if (p < 0 || p >= coord_count[block]) {
    throw std::invalid_argument("ReducedLpColLayout::idxYBlock: bad coordinate");
  }
  return y_offset[block] + block_region * coord_count[block] + p;
}

int ReducedLpColLayout::idxMuR(int block) const {
  if (block < 0 || block >= numBlocks()) {
    throw std::invalid_argument("ReducedLpColLayout::idxMuR: bad block");
  }
  return mu_r_offset + block;
}

int ReducedLpColLayout::idxMuL(int block) const {
  if (block < 0 || block >= numBlocks()) {
    throw std::invalid_argument("ReducedLpColLayout::idxMuL: bad block");
  }
  return mu_l_offset + block;
}

int ReducedLpColLayout::idxU(int k) const {
  if (u_offset < 0) {
    throw std::invalid_argument("ReducedLpColLayout::idxU: no tie-break block");
  }
  if (k < 0 || k >= num_x) {
    throw std::invalid_argument("ReducedLpColLayout::idxU: bad column");
  }
  return u_offset + k;
}

int ReducedLpRowLayout::rowLeft(int block, int block_region,
                                int vertex) const {
  if (block < 0 || block >= numBlocks()) {
    throw std::invalid_argument("ReducedLpRowLayout::rowLeft: bad block");
  }
  const std::vector<int>& offsets = region_pair_offset[block];
  if (block_region < 0
      || block_region + 1 >= static_cast<int>(offsets.size())) {
    throw std::invalid_argument("ReducedLpRowLayout::rowLeft: bad region");
  }
  const int count = offsets[block_region + 1] - offsets[block_region];
  if (vertex < 0 || vertex >= count) {
    throw std::invalid_argument("ReducedLpRowLayout::rowLeft: bad vertex");
  }
  return 2 * (pair_offset[block] + offsets[block_region] + vertex);
}

int ReducedLpRowLayout::rowRight(int block, int block_region,
                                 int vertex) const {
  return rowLeft(block, block_region, vertex) + 1;
}

int ReducedLpRowLayout::rowAbs(int block, int block_region, int p,
                               bool positive) const {
  if (block < 0 || block >= numBlocks()) {
    throw std::invalid_argument("ReducedLpRowLayout::rowAbs: bad block");
  }
  if (block_region < 0 || block_region >= num_regions[block]) {
    throw std::invalid_argument("ReducedLpRowLayout::rowAbs: bad region");
  }
  if (p < 0 || p >= coord_count[block]) {
    throw std::invalid_argument("ReducedLpRowLayout::rowAbs: bad coordinate");
  }
  const int local
      = y_offset_local[block] + block_region * coord_count[block] + p;
  return abs_offset + 2 * local + (positive ? 0 : 1);
}

int ReducedLpRowLayout::rowTieBreak(int k, bool positive) const {
  if (tie_break_offset < 0) {
    throw std::invalid_argument(
        "ReducedLpRowLayout::rowTieBreak: no tie-break rows");
  }
  if (k < 0 || k >= num_x) {
    throw std::invalid_argument("ReducedLpRowLayout::rowTieBreak: bad column");
  }
  return tie_break_offset + 2 * k + (positive ? 0 : 1);
}

void clampBlockRho(ReducedLpInput& input, double tolerance) {
  for (std::size_t b = 0; b < input.blocks.size(); ++b) {
    for (std::size_t j = 0; j < input.blocks[b].regions.size(); ++j) {
      for (BlockVertexData& vertex : input.blocks[b].regions[j].vertices) {
        for (int p = 0; p < vertex.rho.size(); ++p) {
          if (vertex.rho(p) < -tolerance) {
            throw std::runtime_error(
                "block_reduction: negative uncertainty radius in block "
                + std::to_string(b) + " region " + std::to_string(j));
          }
          if (vertex.rho(p) < 0.0) {
            vertex.rho(p) = 0.0;
          }
        }
      }
    }
  }
}

ReducedLpColLayout makeReducedLpColLayout(const ReducedLpInput& input) {
  checkBlocks(input);
  const int num_blocks = static_cast<int>(input.blocks.size());

  ReducedLpColLayout layout;
  layout.num_x = input.num_x;
  layout.coord_count.resize(num_blocks);
  layout.num_regions.resize(num_blocks);
  layout.y_offset.resize(num_blocks);

  int cursor = input.num_x;
  for (int b = 0; b < num_blocks; ++b) {
    layout.coord_count[b] = input.blocks[b].coord_count;
    layout.num_regions[b] = static_cast<int>(input.blocks[b].regions.size());
    layout.y_offset[b] = cursor;
    cursor += layout.coord_count[b] * layout.num_regions[b];
  }
  layout.mu_r_offset = cursor;
  cursor += num_blocks;
  layout.mu_l_offset = cursor;
  cursor += num_blocks;
  if (input.tie_break_eps > 0.0) {
    layout.u_offset = cursor;
    cursor += input.num_x;
  }
  layout.num_cols = cursor;
  return layout;
}

ReducedLpRowLayout makeReducedLpRowLayout(const ReducedLpInput& input) {
  checkBlocks(input);
  const int num_blocks = static_cast<int>(input.blocks.size());

  ReducedLpRowLayout rows;
  rows.num_x = input.num_x;
  rows.coord_count.resize(num_blocks);
  rows.num_regions.resize(num_blocks);
  rows.pair_offset.resize(num_blocks);
  rows.region_pair_offset.resize(num_blocks);
  rows.y_offset_local.resize(num_blocks);

  int pair_cursor = 0;
  int y_cursor = 0;
  for (int b = 0; b < num_blocks; ++b) {
    const ReducedLpBlock& block = input.blocks[b];
    rows.coord_count[b] = block.coord_count;
    rows.num_regions[b] = static_cast<int>(block.regions.size());
    rows.pair_offset[b] = pair_cursor;
    rows.y_offset_local[b] = y_cursor;

    rows.region_pair_offset[b].resize(block.regions.size() + 1);
    int local = 0;
    for (std::size_t j = 0; j < block.regions.size(); ++j) {
      rows.region_pair_offset[b][j] = local;
      local += static_cast<int>(block.regions[j].vertices.size());
    }
    rows.region_pair_offset[b].back() = local;

    pair_cursor += local;
    y_cursor += block.coord_count * static_cast<int>(block.regions.size());
  }

  rows.total_pairs = pair_cursor;
  rows.abs_offset = 2 * pair_cursor;
  rows.coupling_offset = rows.abs_offset + 2 * y_cursor;
  int cursor = rows.coupling_offset + 2;
  if (input.tie_break_eps > 0.0) {
    rows.tie_break_offset = cursor;
    cursor += 2 * input.num_x;
  }
  rows.num_rows = cursor;
  return rows;
}

ReducedLpMatrices assembleReducedLp(const ReducedLpInput& input) {
  checkBlocks(input);

  ReducedLpMatrices out;
  out.cols = makeReducedLpColLayout(input);
  out.rows = makeReducedLpRowLayout(input);

  const int num_blocks = static_cast<int>(input.blocks.size());
  const double s = input.sign_s;
  const double dt = input.t_delta;

  out.starts.assign(1, 0);
  out.row_lower.reserve(out.rows.num_rows);
  out.row_upper.reserve(out.rows.num_rows);
  out.cost = Eigen::RowVectorXd::Zero(out.cols.num_cols);

  // Objective weight of block b. See ObjectiveWeights: with ProductCount the
  // reduced objective reproduces the product objective coefficient for
  // coefficient, because a block-A row occurs in exactly R_B * R_C = R / R_A of
  // the product (region, vertex) pairs.
  std::vector<double> weight(num_blocks, 1.0);
  if (input.weights == ObjectiveWeights::ProductCount) {
    double product = 1.0;
    std::vector<double> pairs(num_blocks);
    for (int b = 0; b < num_blocks; ++b) {
      pairs[b] = static_cast<double>(blockPairCount(input.blocks[b]));
      product *= pairs[b];
    }
    for (int b = 0; b < num_blocks; ++b) {
      weight[b] = product / pairs[b];
    }
  }

  // ---- Residual rows, block by block. ----------------------------------
  //
  // Both groups are written on the plain residual F, with no dt factor left in
  // the row. updateReducedLpRowUpper() below must use the same scale; that is
  // why the two live in one file.
  //
  //   (L_b)  s a_b(a)^T z + rho_b(a)^T y_{b,j} - muL_b <= -s b_b(a)
  //   (R_b) -(s/dt) phi_b(a)^T z          - muR_b <= -s (phi_b^T x/dt + beta_b)
  //
  // with a_b = Psi_{b,j}^T m_b - phi_b/dt and b_b = phi_b^T x/dt + g_b. rho
  // enters (L) with a plus for BOTH directions, because s^2 = 1; that sign is
  // what makes the y-lift monotone and lets y be shared across all product
  // regions with the same j_b.
  for (int b = 0; b < num_blocks; ++b) {
    const ReducedLpBlock& block = input.blocks[b];
    for (int j = 0; j < static_cast<int>(block.regions.size()); ++j) {
      const BlockRegionData& region = block.regions[j];
      for (int k = 0; k < static_cast<int>(region.vertices.size()); ++k) {
        const BlockVertexData& vertex = region.vertices[k];

        // a_b = Psi_{b,j}^T m_b - phi_b / dt.
        SparseVec a;
        for (int p = 0; p < block.coord_count; ++p) {
          const SparseVec& psi_row = region.psi_rows[p];
          for (std::size_t i = 0; i < psi_row.cols.size(); ++i) {
            a.add(psi_row.cols[i], psi_row.vals[i] * vertex.m(p),
                  kAssembleEps);
          }
        }
        for (std::size_t i = 0; i < vertex.phi.cols.size(); ++i) {
          a.add(vertex.phi.cols[i], -vertex.phi.vals[i] / dt, kAssembleEps);
        }

        SparseVec left_row;
        for (std::size_t i = 0; i < a.cols.size(); ++i) {
          left_row.add(a.cols[i], s * a.vals[i], kAssembleEps);
        }
        for (int p = 0; p < block.coord_count; ++p) {
          const double rho = vertex.rho(p);
          if (rho < 0.0) {
            throw std::runtime_error(
                "assembleReducedLp: negative rho; call clampBlockRho first");
          }
          if (rho > 0.0 && rho <= kSmallMatrixValue) {
            // HiGHS would drop this entry, which removes a term that only ever
            // tightens (L). The result would silently stop being a bound.
            throw std::runtime_error(
                "assembleReducedLp: rho = " + std::to_string(rho)
                + " is at or below the HiGHS small_matrix_value");
          }
          left_row.add(out.cols.idxYBlock(b, j, p), rho, 0.0);
        }
        left_row.add(out.cols.idxMuL(b), -1.0, 0.0);

        const int expected_left = out.rows.rowLeft(b, j, k);
        appendRow(out.starts, out.col_index, out.value, left_row);
        out.row_lower.push_back(-kInf);
        // Fixed part of -s b_b; the -s phi_b^T x_next / dt part is per step.
        out.row_upper.push_back(-s * vertex.g);
        if (static_cast<int>(out.row_upper.size()) - 1 != expected_left) {
          throw std::runtime_error("assembleReducedLp: left row id mismatch");
        }

        SparseVec right_row;
        for (std::size_t i = 0; i < vertex.phi.cols.size(); ++i) {
          right_row.add(vertex.phi.cols[i], -s * vertex.phi.vals[i] / dt,
                        kAssembleEps);
        }
        right_row.add(out.cols.idxMuR(b), -1.0, 0.0);

        const int expected_right = out.rows.rowRight(b, j, k);
        appendRow(out.starts, out.col_index, out.value, right_row);
        out.row_lower.push_back(-kInf);
        // Entirely dynamic; recomputed by updateReducedLpRowUpper().
        out.row_upper.push_back(0.0);
        if (static_cast<int>(out.row_upper.size()) - 1 != expected_right) {
          throw std::runtime_error("assembleReducedLp: right row id mismatch");
        }

        // Objective: the l1 norm of the residual at the endpoint with the known
        // value, summed over the whole product. Multiply by the weight before
        // accumulating so the small barycentric coefficients are not lost
        // against a weight of order 1e8.
        const double factor = (s / dt) * weight[b];
        for (std::size_t i = 0; i < vertex.phi.cols.size(); ++i) {
          out.cost(vertex.phi.cols[i]) += factor * vertex.phi.vals[i];
        }
      }
    }
  }

  // ---- Absolute value rows: y_{b,j} >= |Psi_{b,j} z|. -------------------
  for (int b = 0; b < num_blocks; ++b) {
    const ReducedLpBlock& block = input.blocks[b];
    for (int j = 0; j < static_cast<int>(block.regions.size()); ++j) {
      const BlockRegionData& region = block.regions[j];
      for (int p = 0; p < block.coord_count; ++p) {
        for (int sign_index = 0; sign_index < 2; ++sign_index) {
          const bool positive = sign_index == 0;
          const double scale = positive ? 1.0 : -1.0;
          SparseVec row;
          const SparseVec& psi_row = region.psi_rows[p];
          for (std::size_t i = 0; i < psi_row.cols.size(); ++i) {
            row.add(psi_row.cols[i], scale * psi_row.vals[i], kAssembleEps);
          }
          row.add(out.cols.idxYBlock(b, j, p), -1.0, 0.0);

          const int expected = out.rows.rowAbs(b, j, p, positive);
          appendRow(out.starts, out.col_index, out.value, row);
          out.row_lower.push_back(-kInf);
          out.row_upper.push_back(0.0);
          if (static_cast<int>(out.row_upper.size()) - 1 != expected) {
            throw std::runtime_error("assembleReducedLp: abs row id mismatch");
          }
        }
      }
    }
  }

  // ---- Coupling rows: sum_b mu_b <= 0 for each group. -------------------
  // This is the whole reduction. max over the product of a separable function
  // equals the sum of the per-block maxima, so bounding each block maximum by
  // mu_b and asking the three to sum to at most zero is exactly the original
  // family of constraints, not a relaxation of it.
  for (int group = 0; group < 2; ++group) {
    SparseVec row;
    for (int b = 0; b < num_blocks; ++b) {
      row.add(group == 0 ? out.cols.idxMuR(b) : out.cols.idxMuL(b), 1.0, 0.0);
    }
    const int expected = group == 0 ? out.rows.rowSumR() : out.rows.rowSumL();
    appendRow(out.starts, out.col_index, out.value, row);
    out.row_lower.push_back(-kInf);
    out.row_upper.push_back(0.0);
    if (static_cast<int>(out.row_upper.size()) - 1 != expected) {
      throw std::runtime_error("assembleReducedLp: coupling row id mismatch");
    }
  }

  // ---- Tie-break rows: u >= |z|, objective gains eps * sum u. -----------
  if (input.tie_break_eps > 0.0) {
    for (int k = 0; k < input.num_x; ++k) {
      for (int sign_index = 0; sign_index < 2; ++sign_index) {
        const bool positive = sign_index == 0;
        SparseVec row;
        row.add(k, positive ? 1.0 : -1.0, 0.0);
        row.add(out.cols.idxU(k), -1.0, 0.0);

        const int expected = out.rows.rowTieBreak(k, positive);
        appendRow(out.starts, out.col_index, out.value, row);
        out.row_lower.push_back(-kInf);
        out.row_upper.push_back(0.0);
        if (static_cast<int>(out.row_upper.size()) - 1 != expected) {
          throw std::runtime_error(
              "assembleReducedLp: tie-break row id mismatch");
        }
      }
      out.cost(out.cols.idxU(k)) += input.tie_break_eps;
    }
  }

  if (static_cast<int>(out.row_upper.size()) != out.rows.num_rows) {
    throw std::runtime_error("assembleReducedLp: row count mismatch");
  }

  // Column bounds owned here: y and u are moduli and must be nonnegative. The
  // x block is left free; gauge fixing and any box are the caller's business.
  out.col_lower.assign(out.cols.num_cols, -kInf);
  out.col_upper.assign(out.cols.num_cols, kInf);
  for (int b = 0; b < num_blocks; ++b) {
    for (int j = 0; j < out.cols.num_regions[b]; ++j) {
      for (int p = 0; p < out.cols.coord_count[b]; ++p) {
        out.col_lower[out.cols.idxYBlock(b, j, p)] = 0.0;
      }
    }
  }
  if (out.cols.u_offset >= 0) {
    for (int k = 0; k < input.num_x; ++k) {
      out.col_lower[out.cols.idxU(k)] = 0.0;
    }
  }

  return out;
}

std::vector<double> updateReducedLpRowUpper(
    const ReducedLpInput& input, const ReducedLpRowLayout& rows,
    const std::vector<double>& base_upper, const std::vector<double>& x_next) {
  checkBlocks(input);
  if (static_cast<int>(base_upper.size()) != rows.num_rows) {
    throw std::invalid_argument(
        "updateReducedLpRowUpper: base_upper has the wrong length");
  }
  if (static_cast<int>(x_next.size()) != input.num_x) {
    throw std::invalid_argument(
        "updateReducedLpRowUpper: x_next has the wrong length");
  }

  const double s = input.sign_s;
  const double dt = input.t_delta;
  std::vector<double> upper = base_upper;

  for (int b = 0; b < static_cast<int>(input.blocks.size()); ++b) {
    const ReducedLpBlock& block = input.blocks[b];
    for (int j = 0; j < static_cast<int>(block.regions.size()); ++j) {
      const BlockRegionData& region = block.regions[j];
      // q_{b,j} = Psi_{b,j} x_next, shared by every vertex of this region.
      const Eigen::VectorXd q = psiTimes(region, x_next);
      for (int k = 0; k < static_cast<int>(region.vertices.size()); ++k) {
        const BlockVertexData& vertex = region.vertices[k];
        const double phi_x = vertex.phi.dot(x_next);

        // Same scale as the row builder above: plain F, no dt factor.
        const int row_left = rows.rowLeft(b, j, k);
        double left = base_upper[row_left] - s * phi_x / dt;
        if (std::abs(left) <= kRhsSnapEps) {
          left = 0.0;
        }
        upper[row_left] = left;

        const double beta
            = q.dot(vertex.m) + s * vertex.rho.dot(q.cwiseAbs()) + vertex.g;
        double right = -s * (phi_x / dt + beta);
        if (std::abs(right) <= kRhsSnapEps) {
          right = 0.0;
        }
        upper[rows.rowRight(b, j, k)] = right;
      }
    }
  }
  return upper;
}

WorstResidual worstReducedResidual(const ReducedLpInput& input,
                                   const std::vector<double>& x_next,
                                   const std::vector<double>& z) {
  checkBlocks(input);
  if (static_cast<int>(x_next.size()) != input.num_x
      || static_cast<int>(z.size()) != input.num_x) {
    throw std::invalid_argument(
        "worstReducedResidual: x_next and z must have num_x entries");
  }

  const double s = input.sign_s;
  const double dt = input.t_delta;
  WorstResidual out;

  for (const ReducedLpBlock& block : input.blocks) {
    double best_left = -kInf;
    double best_right = -kInf;
    for (const BlockRegionData& region : block.regions) {
      const Eigen::VectorXd q_left = psiTimes(region, z);
      const Eigen::VectorXd q_right = psiTimes(region, x_next);
      for (const BlockVertexData& vertex : region.vertices) {
        const double slope
            = (vertex.phi.dot(x_next) - vertex.phi.dot(z)) / dt;
        const double f_left = slope + q_left.dot(vertex.m)
                              + s * vertex.rho.dot(q_left.cwiseAbs())
                              + vertex.g;
        const double f_right = slope + q_right.dot(vertex.m)
                               + s * vertex.rho.dot(q_right.cwiseAbs())
                               + vertex.g;
        best_left = std::max(best_left, s * f_left);
        best_right = std::max(best_right, s * f_right);
      }
    }
    // Lemma 3: the maximum of a separable function over a product is the sum of
    // the per-factor maxima, so this is the exact global worst case.
    out.left += best_left;
    out.right += best_right;
  }
  return out;
}

}  // namespace block_reduction
}  // namespace barycentric_affine_approximator

// NOLINTEND(readability-identifier-naming)
