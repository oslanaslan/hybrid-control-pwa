#include "util/band_lp.hpp"

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

void appendRow(std::vector<int>& starts, std::vector<int>& cols,
               std::vector<double>& values, const SparseVec& row) {
  for (std::size_t i = 0; i < row.cols.size(); ++i) {
    cols.push_back(row.cols[i]);
    values.push_back(row.vals[i]);
  }
  starts.push_back(static_cast<int>(values.size()));
}

void checkStage(const BandLpRowLayout& rows, int k) {
  if (k < 0 || k >= rows.num_stages) {
    throw std::invalid_argument("BandLpRowLayout: bad stage");
  }
}

void checkStage(const BandLpColLayout& cols, int k) {
  if (k < 0 || k >= cols.num_stages) {
    throw std::invalid_argument("BandLpColLayout: bad stage");
  }
}

void checkInput(const BandLpInput& input) {
  if (input.stage == nullptr) {
    throw std::invalid_argument("assembleBandLp: no stage data");
  }
  if (input.stage->tie_break_eps != 0.0) {
    throw std::invalid_argument(
        "assembleBandLp: the band LP has no tie-break; the integral objective "
        "is strictly increasing in every node value");
  }
  if (input.num_stages < 1) {
    throw std::invalid_argument("assembleBandLp: num_stages must be >= 1");
  }
  const int num_x = input.stage->num_x;
  if (static_cast<int>(input.z_terminal.size()) != num_x) {
    throw std::invalid_argument("assembleBandLp: z_terminal must have num_x "
                                "entries");
  }
  if (static_cast<int>(input.node_weights.size()) != num_x) {
    throw std::invalid_argument(
        "assembleBandLp: node_weights must have num_x entries");
  }
  std::vector<bool> seen(static_cast<std::size_t>(num_x), false);
  for (const int column : input.pinned_columns) {
    if (column < 0 || column >= num_x) {
      throw std::invalid_argument("assembleBandLp: pinned column "
                                  + std::to_string(column)
                                  + " outside the x block");
    }
    if (seen[static_cast<std::size_t>(column)]) {
      throw std::invalid_argument("assembleBandLp: pinned column "
                                  + std::to_string(column) + " listed twice");
    }
    seen[static_cast<std::size_t>(column)] = true;
  }
}

}  // namespace

int BandLpColLayout::idxX(int k, int column) const {
  checkStage(*this, k);
  return k * stage_cols + stage.idxX(column);
}

int BandLpColLayout::idxYBlock(int k, int block, int block_region,
                               int p) const {
  checkStage(*this, k);
  return k * stage_cols + stage.idxYBlock(block, block_region, p);
}

int BandLpColLayout::idxMuR(int k, int block) const {
  checkStage(*this, k);
  return k * stage_cols + stage.idxMuR(block);
}

int BandLpColLayout::idxMuL(int k, int block) const {
  checkStage(*this, k);
  return k * stage_cols + stage.idxMuL(block);
}

int BandLpRowLayout::rowLeft(int k, int block, int block_region,
                             int vertex) const {
  checkStage(*this, k);
  return k * stage_rows + stage.rowLeft(block, block_region, vertex);
}

int BandLpRowLayout::rowRight(int k, int block, int block_region,
                              int vertex) const {
  checkStage(*this, k);
  return k * stage_rows + stage.rowRight(block, block_region, vertex);
}

int BandLpRowLayout::rowAbs(int k, int block, int block_region, int p,
                            bool positive) const {
  checkStage(*this, k);
  return k * stage_rows + stage.rowAbs(block, block_region, p, positive);
}

int BandLpRowLayout::rowSumR(int k) const {
  checkStage(*this, k);
  return k * stage_rows + stage.rowSumR();
}

int BandLpRowLayout::rowSumL(int k) const {
  checkStage(*this, k);
  return k * stage_rows + stage.rowSumL();
}

double trapezoidWeight(int k, int num_stages) {
  if (k < 0 || k > num_stages) {
    throw std::invalid_argument("trapezoidWeight: k outside 0..L");
  }
  return (k == 0 || k == num_stages) ? 0.5 : 1.0;
}

BandLpMatrices assembleBandLp(const BandLpInput& input) {
  checkInput(input);
  const ReducedLpInput& stage = *input.stage;
  const int L = input.num_stages;
  const int num_blocks = static_cast<int>(stage.blocks.size());
  const double s = stage.sign_s;
  const double dt = stage.t_delta;

  BandLpMatrices out;
  // The stage layouts are the one-step layouts (they validate the input), and
  // the stage sizes are theirs, so the rows of stage L-1 are, column for
  // column, the rows of the one-step LP.
  out.cols.stage = makeReducedLpColLayout(stage);
  out.cols.stage_cols = out.cols.stage.num_cols;
  out.cols.num_stages = L;
  out.cols.num_cols = L * out.cols.stage_cols;
  out.rows.stage = makeReducedLpRowLayout(stage);
  out.rows.stage_rows = out.rows.stage.num_rows;
  out.rows.num_stages = L;
  out.rows.num_rows = L * out.rows.stage_rows;

  out.starts.assign(1, 0);
  out.row_lower.reserve(static_cast<std::size_t>(out.rows.num_rows));
  out.row_upper.reserve(static_cast<std::size_t>(out.rows.num_rows));

  auto emit = [&out](const SparseVec& row, double upper, int expected_id,
                     const char* what) {
    appendRow(out.starts, out.col_index, out.value, row);
    out.row_lower.push_back(-kInf);
    out.row_upper.push_back(upper);
    if (static_cast<int>(out.row_upper.size()) - 1 != expected_id) {
      throw std::runtime_error(std::string("assembleBandLp: ") + what
                               + " row id mismatch");
    }
  };

  for (int k = 0; k < L; ++k) {
    const bool terminal = (k + 1 == L);
    const int z_offset = k * out.cols.stage_cols;
    const int z_next_offset = (k + 1) * out.cols.stage_cols;

    // ---- Residual rows of stage k, block by block. ---------------------
    for (int b = 0; b < num_blocks; ++b) {
      const ReducedLpBlock& block = stage.blocks[static_cast<std::size_t>(b)];
      for (int j = 0; j < static_cast<int>(block.regions.size()); ++j) {
        const BlockRegionData& region
            = block.regions[static_cast<std::size_t>(j)];
        // Psi_{b,j} z_L, shared by every vertex of the region at the last
        // stage, where z_L is data.
        Eigen::VectorXd psi_terminal;
        if (terminal) {
          psi_terminal = psiTimes(region, input.z_terminal);
        }
        for (int v = 0; v < static_cast<int>(region.vertices.size()); ++v) {
          const BlockVertexData& vertex
              = region.vertices[static_cast<std::size_t>(v)];
          TerminalRowUpper rhs;
          if (terminal) {
            rhs = terminalRowUpper(vertex, psi_terminal,
                                   vertex.phi.dot(input.z_terminal), s, dt);
          }

          // (L_k): s a^T z_k + (s/dt) phi^T z_{k+1} + rho^T y_k - muL_k.
          SparseVec left
              = leftRowZCoefficients(block, region, vertex, s, dt, z_offset);
          if (!terminal) {
            for (std::size_t i = 0; i < vertex.phi.cols.size(); ++i) {
              left.add(z_next_offset + vertex.phi.cols[i],
                       s * vertex.phi.vals[i] / dt, kAssembleEps);
            }
          }
          addRhoTerms(left, vertex, block.coord_count,
                      out.cols.idxYBlock(k, b, j, 0));
          left.add(out.cols.idxMuL(k, b), -1.0, 0.0);
          emit(left, terminal ? rhs.left : -s * vertex.g,
               out.rows.rowLeft(k, b, j, v), "left");

          // (R_k): -(s/dt) phi^T z_k + s (Psi^T m + phi/dt)^T z_{k+1}
          //        + rho^T y_{k+1} - muR_k.
          SparseVec right;
          for (std::size_t i = 0; i < vertex.phi.cols.size(); ++i) {
            right.add(z_offset + vertex.phi.cols[i],
                      -s * vertex.phi.vals[i] / dt, kAssembleEps);
          }
          if (!terminal) {
            SparseVec psi_m;
            for (int p = 0; p < block.coord_count; ++p) {
              const SparseVec& psi_row
                  = region.psi_rows[static_cast<std::size_t>(p)];
              for (std::size_t i = 0; i < psi_row.cols.size(); ++i) {
                psi_m.add(psi_row.cols[i], psi_row.vals[i] * vertex.m(p),
                          kAssembleEps);
              }
            }
            for (std::size_t i = 0; i < psi_m.cols.size(); ++i) {
              right.add(z_next_offset + psi_m.cols[i], s * psi_m.vals[i],
                        kAssembleEps);
            }
            for (std::size_t i = 0; i < vertex.phi.cols.size(); ++i) {
              right.add(z_next_offset + vertex.phi.cols[i],
                        s * vertex.phi.vals[i] / dt, kAssembleEps);
            }
            // The modulus of the gradient at z_{k+1} is the next stage's y:
            // one column serves the right end of this segment and the left
            // end of the next.
            addRhoTerms(right, vertex, block.coord_count,
                        out.cols.idxYBlock(k + 1, b, j, 0));
          }
          right.add(out.cols.idxMuR(k, b), -1.0, 0.0);
          emit(right, terminal ? rhs.right : -s * vertex.g,
               out.rows.rowRight(k, b, j, v), "right");
        }
      }
    }

    // ---- |.| rows of stage k: y_k >= |Psi z_k|. -----------------------
    for (int b = 0; b < num_blocks; ++b) {
      const ReducedLpBlock& block = stage.blocks[static_cast<std::size_t>(b)];
      for (int j = 0; j < static_cast<int>(block.regions.size()); ++j) {
        const BlockRegionData& region
            = block.regions[static_cast<std::size_t>(j)];
        for (int p = 0; p < block.coord_count; ++p) {
          for (int sign_index = 0; sign_index < 2; ++sign_index) {
            const bool positive = sign_index == 0;
            const double scale = positive ? 1.0 : -1.0;
            SparseVec row;
            const SparseVec& psi_row
                = region.psi_rows[static_cast<std::size_t>(p)];
            for (std::size_t i = 0; i < psi_row.cols.size(); ++i) {
              row.add(z_offset + psi_row.cols[i], scale * psi_row.vals[i],
                      kAssembleEps);
            }
            row.add(out.cols.idxYBlock(k, b, j, p), -1.0, 0.0);
            emit(row, 0.0, out.rows.rowAbs(k, b, j, p, positive), "abs");
          }
        }
      }
    }

    // ---- Coupling rows of stage k: sum_b mu_b <= 0 for each end. ------
    for (int group = 0; group < 2; ++group) {
      SparseVec row;
      for (int b = 0; b < num_blocks; ++b) {
        row.add(group == 0 ? out.cols.idxMuR(k, b) : out.cols.idxMuL(k, b),
                1.0, 0.0);
      }
      emit(row, 0.0, group == 0 ? out.rows.rowSumR(k) : out.rows.rowSumL(k),
           "coupling");
    }
  }

  if (static_cast<int>(out.row_upper.size()) != out.rows.num_rows) {
    throw std::runtime_error("assembleBandLp: row count mismatch");
  }

  // ---- Objective: s * omega_k q on the x columns of stage k. -----------
  out.cost = Eigen::RowVectorXd::Zero(out.cols.num_cols);
  for (int k = 0; k < L; ++k) {
    const double omega = trapezoidWeight(k, L);
    for (int i = 0; i < stage.num_x; ++i) {
      out.cost(out.cols.idxX(k, i)) = s * omega * input.node_weights(i);
    }
  }
  out.cost_scale = out.cost.cwiseAbs().maxCoeff();
  if (!(out.cost_scale > 0.0) || !std::isfinite(out.cost_scale)) {
    throw std::runtime_error(
        "assembleBandLp: node weights are all zero or not finite");
  }

  // ---- Column bounds. -------------------------------------------------
  out.col_lower.assign(static_cast<std::size_t>(out.cols.num_cols), -kInf);
  out.col_upper.assign(static_cast<std::size_t>(out.cols.num_cols), kInf);
  for (int k = 0; k < L; ++k) {
    for (int b = 0; b < num_blocks; ++b) {
      for (int j = 0; j < out.cols.stage.num_regions[static_cast<std::size_t>(b)];
           ++j) {
        for (int p = 0; p < out.cols.stage.coord_count[static_cast<std::size_t>(b)];
             ++p) {
          out.col_lower[static_cast<std::size_t>(
              out.cols.idxYBlock(k, b, j, p))] = 0.0;
        }
      }
    }
    for (const int column : input.pinned_columns) {
      const int idx = out.cols.idxX(k, column);
      out.col_lower[static_cast<std::size_t>(idx)] = 0.0;
      out.col_upper[static_cast<std::size_t>(idx)] = 0.0;
    }
  }

  return out;
}

double bandIntegral(const Eigen::VectorXd& node_weights, double t_delta,
                    const std::vector<std::vector<double>>& z) {
  if (z.empty()) {
    throw std::invalid_argument("bandIntegral: no node values");
  }
  const int L = static_cast<int>(z.size()) - 1;
  if (L == 0) {
    return 0.0;
  }
  double sum = 0.0;
  for (int k = 0; k <= L; ++k) {
    const std::vector<double>& z_k = z[static_cast<std::size_t>(k)];
    if (static_cast<int>(z_k.size()) != node_weights.size()) {
      throw std::invalid_argument(
          "bandIntegral: node value vector " + std::to_string(k)
          + " does not match the weights");
    }
    double dot = 0.0;
    for (int i = 0; i < node_weights.size(); ++i) {
      dot += node_weights(i) * z_k[static_cast<std::size_t>(i)];
    }
    sum += trapezoidWeight(k, L) * dot;
  }
  return t_delta * sum;
}

std::vector<WorstResidual> bandStageResiduals(
    const ReducedLpInput& stage, const std::vector<std::vector<double>>& z) {
  if (z.size() < 2) {
    throw std::invalid_argument("bandStageResiduals: need at least one stage");
  }
  std::vector<WorstResidual> out;
  out.reserve(z.size() - 1);
  for (std::size_t k = 0; k + 1 < z.size(); ++k) {
    out.push_back(worstReducedResidual(stage, z[k + 1], z[k]));
  }
  return out;
}

}  // namespace block_reduction
}  // namespace barycentric_affine_approximator

// NOLINTEND(readability-identifier-naming)
