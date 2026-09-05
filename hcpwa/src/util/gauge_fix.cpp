#include "util/gauge_fix.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <format>
#include <limits>
#include <stdexcept>
#include <string>

// NOLINTBEGIN(readability-identifier-naming)

namespace barycentric_affine_approximator {

namespace {

// Relative threshold below which a singular value counts as zero. Same as
// the gauge diagnostic that first measured the kernel.
constexpr double kRankTol = 1e-10;

std::vector<int> sortedUnique(std::vector<int> ids) {
  std::sort(ids.begin(), ids.end());
  ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
  return ids;
}

// Is the line {coordinate axis = c} a union of edges of the layer's
// triangulation? Exactly when no triangle straddles it.
bool lineIsEdgeUnion(const ProjectionLayer& layer, int axis, double c) {
  for (const TriangleBasis& basis : layer.bases) {
    double lo = std::numeric_limits<double>::infinity();
    double hi = -std::numeric_limits<double>::infinity();
    for (const int id : basis.vertex_ids) {
      const double value
          = layer.unique_vertices[static_cast<std::size_t>(id)](axis);
      lo = std::min(lo, value);
      hi = std::max(hi, value);
    }
    if (lo < c - kGeomEps && hi > c + kGeomEps) {
      return false;
    }
  }
  return true;
}

// Index of `value` in `axes`, or -1.
int axisPosition(const std::array<int, 2>& axes, int value) {
  if (axes[0] == value) {
    return 0;
  }
  if (axes[1] == value) {
    return 1;
  }
  return -1;
}

int numericalRank(const Eigen::MatrixXd& m, double tol_relative) {
  if (m.rows() == 0 || m.cols() == 0) {
    return 0;
  }
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(m);
  const Eigen::VectorXd sv = svd.singularValues();
  int rank = 0;
  for (int i = 0; i < sv.size(); ++i) {
    if (sv(i) > sv(0) * tol_relative) {
      ++rank;
    }
  }
  return rank;
}

Eigen::MatrixXd restrictRows(const Eigen::MatrixXd& m,
                             const std::vector<int>& rows) {
  Eigen::MatrixXd out(static_cast<int>(rows.size()), m.cols());
  for (std::size_t i = 0; i < rows.size(); ++i) {
    out.row(static_cast<int>(i)) = m.row(rows[i]);
  }
  return out;
}

}  // namespace

std::vector<int> GaugeFix::pinsForLevel(int switch_cnt) const {
  if (switch_cnt < 0) {
    throw std::invalid_argument("GaugeFix::pinsForLevel: negative level");
  }
  std::vector<int> pins = line_pins;
  const std::vector<int>& extra
      = switch_cnt == 0 ? group_c_columns : constant_pins;
  pins.insert(pins.end(), extra.begin(), extra.end());
  return sortedUnique(std::move(pins));
}

GaugeFix buildGaugeFix(const PhaseGeometry& geometry,
                       const BarycentricVarLayout& layout, int phase) {
  const auto axes = projectionAxesForPhase(phase);
  GaugeFix gauge;

  for (int b = 0; b < kBlockCount; ++b) {
    const BlockGeometry& block = geometry.blocks[static_cast<std::size_t>(b)];
    if (block.layer_count == 1) {
      if (gauge.group_c_layer >= 0) {
        throw std::runtime_error(
            "buildGaugeFix: two single-plane blocks; expected exactly one");
      }
      gauge.group_c_block = b;
      gauge.group_c_layer = block.layer_ids[0];
      for (int k = 0; k < layout.eta_s[static_cast<std::size_t>(
               gauge.group_c_layer)]; ++k) {
        gauge.group_c_columns.push_back(layout.idxX(gauge.group_c_layer, k));
      }
      continue;
    }
    if (block.layer_count != 2) {
      throw std::runtime_error(std::format(
          "buildGaugeFix: block {} of phase {} owns {} planes; the kernel "
          "argument covers one or two",
          b, phase, block.layer_count));
    }

    GaugeFix::Pair pair;
    pair.block = b;
    pair.layer_plus = block.layer_ids[0];
    pair.layer_minus = block.layer_ids[1];
    const auto& axes_plus = axes[static_cast<std::size_t>(pair.layer_plus)];
    const auto& axes_minus = axes[static_cast<std::size_t>(pair.layer_minus)];
    int shared_pos_plus = -1;
    int shared_pos_minus = -1;
    int shared_count = 0;
    for (int i = 0; i < 2; ++i) {
      const int pos = axisPosition(axes_minus, axes_plus[i]);
      if (pos >= 0) {
        ++shared_count;
        shared_pos_plus = i;
        shared_pos_minus = pos;
      }
    }
    if (shared_count != 1) {
      throw std::runtime_error(std::format(
          "buildGaugeFix: planes {} and {} of phase {} share {} coordinates; "
          "the kernel argument needs exactly one",
          pair.layer_plus, pair.layer_minus, phase, shared_count));
    }
    pair.shared_coord = axes_plus[shared_pos_plus];
    const int other_pos_plus = 1 - shared_pos_plus;

    const ProjectionLayer& plus
        = geometry.layers[static_cast<std::size_t>(pair.layer_plus)];
    const ProjectionLayer& minus
        = geometry.layers[static_cast<std::size_t>(pair.layer_minus)];

    // Candidate values of the shared coordinate: vertex coordinates present
    // in both planes. A line that is a union of edges has its edge endpoints
    // among the vertices, so nothing is missed by restricting to these.
    std::vector<double> candidates;
    for (const Eigen::Vector2d& v : plus.unique_vertices) {
      const double c = v(shared_pos_plus);
      bool in_minus = false;
      for (const Eigen::Vector2d& w : minus.unique_vertices) {
        if (std::abs(w(shared_pos_minus) - c) <= kGeomEps) {
          in_minus = true;
          break;
        }
      }
      if (!in_minus) {
        continue;
      }
      bool seen = false;
      for (const double known : candidates) {
        if (std::abs(known - c) <= kGeomEps) {
          seen = true;
          break;
        }
      }
      if (!seen) {
        candidates.push_back(c);
      }
    }
    for (const double c : candidates) {
      if (lineIsEdgeUnion(plus, shared_pos_plus, c)
          && lineIsEdgeUnion(minus, shared_pos_minus, c)) {
        pair.xi.push_back(c);
      }
    }
    std::sort(pair.xi.begin(), pair.xi.end());
    if (pair.xi.size() < 2) {
      throw std::runtime_error(std::format(
          "buildGaugeFix: planes {} and {} of phase {} share only {} edge "
          "lines; the two sides of the square alone give two",
          pair.layer_plus, pair.layer_minus, phase, pair.xi.size()));
    }

    // Line pins: the nodes of the first plane on the boundary line
    // {other coordinate = 0} at the values of Xi. The boundary is a union of
    // edges of any triangulation of the square, and every line of Xi runs
    // across the whole square, so each such node exists.
    for (const double c : pair.xi) {
      int found = -1;
      for (int k = 0; k < static_cast<int>(plus.unique_vertices.size()); ++k) {
        const Eigen::Vector2d& v = plus.unique_vertices[static_cast<std::size_t>(k)];
        if (std::abs(v(other_pos_plus)) <= kGeomEps
            && std::abs(v(shared_pos_plus) - c) <= kGeomEps) {
          if (found >= 0) {
            throw std::runtime_error(std::format(
                "buildGaugeFix: plane {} of phase {} has two nodes at the "
                "boundary point with shared coordinate {}",
                pair.layer_plus, phase, c));
          }
          found = k;
        }
      }
      if (found < 0) {
        throw std::runtime_error(std::format(
            "buildGaugeFix: plane {} of phase {} has no node on the boundary "
            "at shared coordinate {}, although that line is a union of edges",
            pair.layer_plus, phase, c));
      }
      gauge.line_pins.push_back(layout.idxX(pair.layer_plus, found));
    }
    if (layout.eta_s[static_cast<std::size_t>(pair.layer_minus)] == 0) {
      throw std::runtime_error("buildGaugeFix: empty projection plane");
    }
    gauge.constant_pins.push_back(layout.idxX(pair.layer_minus, 0));
    gauge.pairs.push_back(std::move(pair));
  }

  if (gauge.pairs.size() != 2 || gauge.group_c_layer < 0) {
    throw std::runtime_error(std::format(
        "buildGaugeFix: phase {} has {} plane pairs and {} single planes; "
        "expected two and one",
        phase, gauge.pairs.size(), gauge.group_c_layer < 0 ? 0 : 1));
  }
  gauge.line_pins = sortedUnique(std::move(gauge.line_pins));
  gauge.constant_pins = sortedUnique(std::move(gauge.constant_pins));
  gauge.group_c_columns = sortedUnique(std::move(gauge.group_c_columns));
  return gauge;
}

GaugeFixReport verifyGaugeFix(const GaugeFix& gauge,
                              const block_reduction::ReducedLpInput& input,
                              const BarycentricVarLayout& layout) {
  const int num_x = layout.num_x;
  if (input.num_x != num_x) {
    throw std::invalid_argument("verifyGaugeFix: input and layout disagree "
                                "on num_x");
  }

  // Every row that depends on z: the value selectors of all block vertices
  // and the gradient rows of all block cells.
  long long total_rows = 0;
  for (const auto& block : input.blocks) {
    for (const auto& region : block.regions) {
      total_rows += static_cast<long long>(region.psi_rows.size());
      total_rows += static_cast<long long>(region.vertices.size());
    }
  }
  Eigen::MatrixXd A
      = Eigen::MatrixXd::Zero(static_cast<int>(total_rows), num_x);
  int r = 0;
  for (const auto& block : input.blocks) {
    for (const auto& region : block.regions) {
      for (const auto& psi : region.psi_rows) {
        for (std::size_t i = 0; i < psi.cols.size(); ++i) {
          A(r, psi.cols[i]) = psi.vals[i];
        }
        ++r;
      }
      for (const auto& vertex : region.vertices) {
        for (std::size_t i = 0; i < vertex.phi.cols.size(); ++i) {
          A(r, vertex.phi.cols[i]) = vertex.phi.vals[i];
        }
        ++r;
      }
    }
  }

  Eigen::BDCSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinV);
  const Eigen::VectorXd sv = svd.singularValues();
  GaugeFixReport report;
  report.sigma_max = sv(0);
  const double tol = sv(0) * kRankTol;
  int rank = 0;
  for (int i = 0; i < sv.size(); ++i) {
    if (sv(i) > tol) {
      ++rank;
      report.sigma_min_kept = sv(i);
    } else {
      report.sigma_max_dropped = std::max(report.sigma_max_dropped, sv(i));
    }
  }
  report.kernel_dim = num_x - rank;
  const Eigen::MatrixXd K0 = svd.matrixV().rightCols(report.kernel_dim);

  // What the remark predicts: one kernel function per pair, of dimension |Xi|.
  int predicted = 0;
  std::string xi_text;
  for (const GaugeFix::Pair& pair : gauge.pairs) {
    predicted += static_cast<int>(pair.xi.size());
    xi_text += std::format(" planes {}/{}: |Xi| = {} (", pair.layer_plus,
                           pair.layer_minus, pair.xi.size());
    for (const double c : pair.xi) {
      xi_text += std::format("{:.6g} ", c);
    }
    xi_text += ")";
  }
  if (report.kernel_dim != predicted) {
    throw std::runtime_error(std::format(
        "verifyGaugeFix: the kernel of the value and gradient rows has "
        "dimension {} but the shared edge lines predict {}:{}; singular "
        "values: max {:.3e}, smallest kept {:.3e}, largest dropped {:.3e}",
        report.kernel_dim, predicted, xi_text, report.sigma_max,
        report.sigma_min_kept, report.sigma_max_dropped));
  }

  // The kernel lives on the two plane pairs; the C plane's hat functions are
  // linearly independent, so it must not appear there.
  double on_c = 0.0;
  for (const int col : gauge.group_c_columns) {
    on_c = std::max(on_c, K0.row(col).cwiseAbs().maxCoeff());
  }
  if (on_c > 1e-8) {
    throw std::runtime_error(std::format(
        "verifyGaugeFix: a kernel direction has a component {:.3e} on the "
        "block-C plane, which should carry none",
        on_c));
  }

  // The constant exchanges of the coupled form: +1 on the second plane of a
  // pair, -1 on the C plane. Not kernel directions of the rows -- V_A moves
  // by the constant -- but null directions of the coupled LP.
  Eigen::MatrixXd K(num_x, report.kernel_dim + 2);
  K.leftCols(report.kernel_dim) = K0;
  for (std::size_t i = 0; i < gauge.pairs.size(); ++i) {
    Eigen::VectorXd d = Eigen::VectorXd::Zero(num_x);
    const int s_minus = gauge.pairs[i].layer_minus;
    for (int k = 0; k < layout.eta_s[static_cast<std::size_t>(s_minus)]; ++k) {
      d(layout.idxX(s_minus, k)) = 1.0;
    }
    for (const int col : gauge.group_c_columns) {
      d(col) = -1.0;
    }
    K.col(report.kernel_dim + static_cast<int>(i)) = d;
  }

  report.line_rank = numericalRank(restrictRows(K0, gauge.line_pins), kRankTol);
  std::vector<int> both = gauge.line_pins;
  both.insert(both.end(), gauge.constant_pins.begin(),
              gauge.constant_pins.end());
  both = sortedUnique(std::move(both));
  report.full_rank = numericalRank(restrictRows(K, both), kRankTol);

  if (report.line_rank != report.kernel_dim) {
    throw std::runtime_error(std::format(
        "verifyGaugeFix: {} line pins meet the {}-dimensional kernel with "
        "rank {} only; some kernel direction survives the normalisation:{}",
        gauge.line_pins.size(), report.kernel_dim, report.line_rank, xi_text));
  }
  if (report.full_rank != report.kernel_dim + 2) {
    throw std::runtime_error(std::format(
        "verifyGaugeFix: line and constant pins ({} columns) meet the kernel "
        "extended by the two constant exchanges with rank {} instead of {}",
        both.size(), report.full_rank, report.kernel_dim + 2));
  }
  return report;
}

}  // namespace barycentric_affine_approximator

// NOLINTEND(readability-identifier-naming)
