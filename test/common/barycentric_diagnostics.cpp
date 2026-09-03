// Diagnostics for the geometric growth of the barycentric value function.
//
// These are experiments, not regression tests, and they are DISABLED_ so the
// normal suite does not pay for them. Run one with
//   ./test_common --gtest_also_run_disabled_tests \
//                 --gtest_filter='common.DISABLED_<name>'
//
// Each one answers a specific question that the production logs could not:
//   gauge_nullspace     -- does the constraint system have a null direction the
//                          gauge fixing does not remove?
//   step_amplification  -- does one backward step of the reduced LP multiply the
//                          solution, and does halving dt halve the effect?
//   border_vs_candidates-- does the border LP's output stay under the candidates
//                          it is supposed to be bounded by?

#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <Highs.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <format>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include <barycentric_affine_approximator.hpp>
#include <hcpwa.hpp>
#include "cddwrap/cdd.hpp"
#include "utility.hpp"

namespace {

namespace baa = barycentric_affine_approximator;
namespace br = barycentric_affine_approximator::block_reduction;

// Production parameters, identical to test/common/barycentric_affine_solver.cpp
// so the numbers here are directly comparable to the run logs.
constexpr double kN = 160;

baa::BarycentricAffineApproximator makeApproximator(int t_split_count) {
  const baa::SystemParams system_params{
      kN,   0.5,  0.017, 0.0055, 0.6,  0.4,  0.8,  0.2,  0.6, 0.4,
      0.7,  0.3,  0.44,  0.23,   0.09, 0.23, 0.46, 0.27, 0.11, 0.27};
  return baa::BarycentricAffineApproximator(
      /*t_max=*/1200.0, t_split_count, /*tau_min=*/10, /*tau_max=*/50,
      system_params, /*highs_verbose=*/false, baa::ApproximationMode::Lower);
}

// The candidate null direction of hypothesis 2: within one coordinate block two
// projection planes share a coordinate, and the globally affine function
// alpha * n_shared can be added to one plane's nodes and subtracted from the
// other's. Both are exactly representable, so V is unchanged.
//
// Returns one direction per (block, shared coordinate), in the x basis.
struct NullCandidate {
  std::vector<double> d;
  int block = 0;
  int layer_plus = 0;
  int layer_minus = 0;
  int shared_coord = 0;
};

std::vector<NullCandidate> buildNullCandidates(
    const baa::PhaseGeometry& geometry, const baa::BarycentricVarLayout& layout,
    int phase) {
  const auto axes = baa::projectionAxesForPhase(phase);
  std::vector<NullCandidate> out;

  for (int b = 0; b < baa::kBlockCount; ++b) {
    const baa::BlockGeometry& block = geometry.blocks[static_cast<std::size_t>(b)];
    for (int l1 = 0; l1 < block.layer_count; ++l1) {
      for (int l2 = l1 + 1; l2 < block.layer_count; ++l2) {
        const int s1 = block.layer_ids[static_cast<std::size_t>(l1)];
        const int s2 = block.layer_ids[static_cast<std::size_t>(l2)];
        // The coordinate the two planes share, if any.
        for (int a1 = 0; a1 < 2; ++a1) {
          for (int a2 = 0; a2 < 2; ++a2) {
            if (axes[static_cast<std::size_t>(s1)][static_cast<std::size_t>(a1)]
                != axes[static_cast<std::size_t>(s2)][static_cast<std::size_t>(a2)]) {
              continue;
            }
            NullCandidate cand;
            cand.block = b;
            cand.layer_plus = s1;
            cand.layer_minus = s2;
            cand.shared_coord
                = axes[static_cast<std::size_t>(s1)][static_cast<std::size_t>(a1)];
            cand.d.assign(static_cast<std::size_t>(layout.num_x), 0.0);
            // +alpha * n_shared on plane s1, -alpha * n_shared on plane s2,
            // with alpha = 1. The node's shared coordinate is component a1 (a2)
            // of its 2D position on that plane.
            const auto& v1 = geometry.layers[static_cast<std::size_t>(s1)]
                                 .unique_vertices;
            for (int k = 0; k < static_cast<int>(v1.size()); ++k) {
              cand.d[static_cast<std::size_t>(layout.idxX(s1, k))]
                  = v1[static_cast<std::size_t>(k)](a1);
            }
            const auto& v2 = geometry.layers[static_cast<std::size_t>(s2)]
                                 .unique_vertices;
            for (int k = 0; k < static_cast<int>(v2.size()); ++k) {
              cand.d[static_cast<std::size_t>(layout.idxX(s2, k))]
                  = -v2[static_cast<std::size_t>(k)](a2);
            }
            out.push_back(std::move(cand));
          }
        }
      }
    }
  }
  return out;
}

// Builds a HiGHS model over an arbitrary ReducedLpInput, exactly the way
// initializeHighs() does: same options, same objective normalisation, same
// gauge fixing. Returns the solver and the row bounds.
struct StepModel {
  std::unique_ptr<Highs> highs;
  br::ReducedLpMatrices matrices;
  br::ReducedLpRowLayout rows;
  std::vector<double> row_lower;
  std::vector<double> base_upper;
};

StepModel buildStepModel(const br::ReducedLpInput& input,
                         const baa::BarycentricVarLayout& layout) {
  StepModel model;
  model.matrices = br::assembleReducedLp(input);
  model.rows = model.matrices.rows;
  model.row_lower = model.matrices.row_lower;
  model.base_upper = model.matrices.row_upper;

  model.highs = std::make_unique<Highs>();
  model.highs->setOptionValue("solver", "simplex");
  model.highs->setOptionValue("presolve", "on");
  model.highs->setOptionValue("simplex_strategy", 2);
  model.highs->setOptionValue("primal_feasibility_tolerance", 1e-6);
  model.highs->setOptionValue("dual_feasibility_tolerance", 1e-6);
  model.highs->setOptionValue("small_matrix_value", br::kSmallMatrixValue);
  model.highs->setOptionValue("simplex_iteration_limit", 400000);
  model.highs->setOptionValue("log_to_console", false);
  model.highs->changeObjectiveSense(ObjSense::kMinimize);

  std::vector<double> cost(
      model.matrices.cost.data(),
      model.matrices.cost.data() + model.matrices.cost.size());
  double scale = 0.0;
  for (const double c : cost) {
    scale = std::max(scale, std::abs(c));
  }
  for (double& c : cost) {
    c /= scale;
  }

  std::vector<double> col_lower = model.matrices.col_lower;
  std::vector<double> col_upper = model.matrices.col_upper;
  for (int s = 1; s < baa::kSubsystemCount; ++s) {
    const int col = layout.idxX(s, 0);
    col_lower[static_cast<std::size_t>(col)] = 0.0;
    col_upper[static_cast<std::size_t>(col)] = 0.0;
  }

  model.highs->addCols(model.matrices.cols.num_cols, cost.data(),
                       col_lower.data(), col_upper.data(), 0, nullptr, nullptr,
                       nullptr);
  model.highs->addRows(model.rows.num_rows, model.matrices.row_lower.data(),
                       model.matrices.row_upper.data(),
                       static_cast<int>(model.matrices.value.size()),
                       model.matrices.starts.data(),
                       model.matrices.col_index.data(),
                       model.matrices.value.data());
  return model;
}

// One backward step: install the row bounds for x_next and solve.
std::vector<double> stepOnce(StepModel& model, const br::ReducedLpInput& input,
                             int num_x, const std::vector<double>& x_next) {
  const std::vector<double> upper = br::updateReducedLpRowUpper(
      input, model.rows, model.base_upper, x_next);
  std::vector<int> ids(upper.size());
  for (std::size_t i = 0; i < ids.size(); ++i) {
    ids[i] = static_cast<int>(i);
  }
  model.highs->changeRowsBounds(static_cast<int>(ids.size()), ids.data(),
                                model.row_lower.data(), upper.data());
  model.highs->run();
  const auto& sol = model.highs->getSolution();
  return std::vector<double>(sol.col_value.begin(),
                             sol.col_value.begin() + num_x);
}

double maxAbs(const std::vector<double>& v) {
  double out = 0.0;
  for (const double x : v) {
    out = std::max(out, std::abs(x));
  }
  return out;
}

double integralOfV(const Eigen::VectorXd& w, const std::vector<double>& x) {
  double out = 0.0;
  for (int i = 0; i < w.size(); ++i) {
    out += w(i) * x[static_cast<std::size_t>(i)];
  }
  return out;
}

}  // namespace

// ---------------------------------------------------------------------------
// Hypothesis 2: the gauge fixing leaves a null direction.
//
// Checks the candidate direction against EVERY row of the reduced LP -- every
// phi row of every block vertex, every Psi row of every block region -- and
// against the objective. Not a sample: the whole system.
//
// Then computes the rank of the full row matrix to see whether there are more
// null directions than the ones predicted.
// ---------------------------------------------------------------------------
TEST(common, DISABLED_gauge_nullspace) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  auto approx = makeApproximator(240);
  approx.getIntersectionPoints();
  approx.precomputeMatrices();

  for (int phase = 0; phase < baa::kPhases; ++phase) {
    const auto& geometry = approx.phaseGeometries()[static_cast<std::size_t>(phase)];
    const auto& layout = approx.layouts()[static_cast<std::size_t>(phase)];
    const br::ReducedLpInput& input = approx.reducedLpInput(phase);
    const br::ReducedLpMatrices matrices = br::assembleReducedLp(input);

    std::cout << std::format("\n=== phase {}: num_x = {} ===\n", phase,
                             layout.num_x);

    const std::vector<NullCandidate> candidates
        = buildNullCandidates(geometry, layout, phase);
    std::cout << std::format("predicted null directions: {}\n",
                             candidates.size());

    for (const NullCandidate& cand : candidates) {
      double worst_phi = 0.0;
      double worst_psi = 0.0;
      long long phi_rows = 0;
      long long psi_rows = 0;
      double scale = 0.0;
      for (const auto& block : input.blocks) {
        for (const auto& region : block.regions) {
          for (const auto& psi : region.psi_rows) {
            worst_psi = std::max(worst_psi, std::abs(psi.dot(cand.d)));
            ++psi_rows;
          }
          for (const auto& vertex : region.vertices) {
            worst_phi = std::max(worst_phi, std::abs(vertex.phi.dot(cand.d)));
            ++phi_rows;
            for (std::size_t i = 0; i < vertex.phi.cols.size(); ++i) {
              scale = std::max(
                  scale,
                  std::abs(vertex.phi.vals[i]
                           * cand.d[static_cast<std::size_t>(
                               vertex.phi.cols[i])]));
            }
          }
        }
      }
      double obj = 0.0;
      for (int k = 0; k < layout.num_x; ++k) {
        obj += matrices.cost(k) * cand.d[static_cast<std::size_t>(k)];
      }
      // Does the gauge fixing kill it? The pinned columns are idxX(s, 0),
      // s = 1..4; a direction with a nonzero there is removed by the pin.
      double pinned = 0.0;
      for (int s = 1; s < baa::kSubsystemCount; ++s) {
        pinned = std::max(
            pinned,
            std::abs(cand.d[static_cast<std::size_t>(layout.idxX(s, 0))]));
      }

      std::cout << std::format(
          "  block {} planes {}/{} shared coord {}:\n"
          "    worst |phi^T d| = {:.3e} over {} rows (terms up to {:.3e})\n"
          "    worst |Psi^T d| = {:.3e} over {} rows\n"
          "    |cost^T d|      = {:.3e}   (|cost|inf = {:.3e})\n"
          "    max |d| on a gauge-pinned column = {:.3e}\n",
          cand.block, cand.layer_plus, cand.layer_minus, cand.shared_coord,
          worst_phi, phi_rows, scale, worst_psi, psi_rows, std::abs(obj),
          matrices.cost.cwiseAbs().maxCoeff(), pinned);
    }

    // Rank of the whole z-dependent row system.
    long long total_rows = 0;
    for (const auto& block : input.blocks) {
      for (const auto& region : block.regions) {
        total_rows += static_cast<long long>(region.psi_rows.size());
        total_rows += static_cast<long long>(region.vertices.size());
      }
    }
    Eigen::MatrixXd A
        = Eigen::MatrixXd::Zero(static_cast<int>(total_rows), layout.num_x);
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
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinV);
    const Eigen::VectorXd sv = svd.singularValues();
    const double tol = sv(0) * 1e-10;
    int rank = 0;
    for (int i = 0; i < sv.size(); ++i) {
      if (sv(i) > tol) {
        ++rank;
      }
    }
    std::cout << std::format(
        "  row matrix {}x{}: rank = {}, nullity = {}\n"
        "  smallest five singular values: {:.3e} {:.3e} {:.3e} {:.3e} {:.3e}\n"
        "  largest singular value: {:.3e}\n",
        A.rows(), A.cols(), rank, layout.num_x - rank,
        sv(sv.size() - 1), sv(sv.size() - 2), sv(sv.size() - 3),
        sv(sv.size() - 4), sv(sv.size() - 5), sv(0));

    // The decisive question: does the gauge fixing kill that null space?
    // A direction survives only if it is zero on every pinned column, so the
    // surviving subspace is ker(A) intersect {d : d[pinned] = 0}. Its dimension
    // is nullity minus the rank of the kernel basis restricted to the pinned
    // columns.
    const int nullity = layout.num_x - rank;
    Eigen::MatrixXd K(nullity, baa::kSubsystemCount - 1);
    for (int i = 0; i < nullity; ++i) {
      for (int s = 1; s < baa::kSubsystemCount; ++s) {
        K(i, s - 1) = svd.matrixV()(layout.idxX(s, 0),
                                    layout.num_x - nullity + i);
      }
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> ksvd(K);
    const Eigen::VectorXd ksv = ksvd.singularValues();
    int krank = 0;
    for (int i = 0; i < ksv.size(); ++i) {
      if (ksv(i) > ksv(0) * 1e-10) {
        ++krank;
      }
    }
    std::cout << std::format(
        "  gauge pins {} columns; kernel basis restricted to them has rank {}\n"
        "  singular values of that {}x{} block: ",
        baa::kSubsystemCount - 1, krank, K.rows(), K.cols());
    for (int i = 0; i < ksv.size(); ++i) {
      std::cout << std::format("{:.3e} ", ksv(i));
    }
    std::cout << std::format(
        "\n  >>> null directions surviving the gauge fixing: {}\n",
        nullity - krank);

    // The decisive follow-up. The per-level growth of x is measured on
    // coefficients, which drift is free to inflate. The growth actually
    // reported per level in the run logs is the courier master's objective,
    // w^T z with w = node_weights, and w_j is the integral of basis function j.
    // If w annihilates the surviving directions, that objective is a property
    // of the represented function alone, and no amount of drift can move it --
    // which would mean the gauge cannot be the cause of the level-to-level
    // growth, only of the coefficient blow-up.
    const Eigen::VectorXd& w
        = approx.nodeWeights()[static_cast<std::size_t>(phase)];
    std::cout << "  node-weight objective on the kernel basis (relative):\n";
    for (int i = 0; i < nullity; ++i) {
      Eigen::VectorXd d(layout.num_x);
      for (int k = 0; k < layout.num_x; ++k) {
        d(k) = svd.matrixV()(k, layout.num_x - nullity + i);
      }
      const double num = std::abs(w.head(layout.num_x).dot(d));
      const double den = w.head(layout.num_x).norm() * d.norm();
      std::cout << std::format("    kernel vector {}: |w^T d| = {:.3e}, "
                               "|w^T d| / (|w| |d|) = {:.3e}\n",
                               i, num, den > 0.0 ? num / den : 0.0);
    }
    for (const NullCandidate& cand : candidates) {
      Eigen::VectorXd d(layout.num_x);
      for (int k = 0; k < layout.num_x; ++k) {
        d(k) = cand.d[static_cast<std::size_t>(k)];
      }
      const double num = std::abs(w.head(layout.num_x).dot(d));
      const double den = w.head(layout.num_x).norm() * d.norm();
      std::cout << std::format("    predicted dir (block {} planes {}/{}): "
                               "|w^T d| / (|w| |d|) = {:.3e}\n",
                               cand.block, cand.layer_plus, cand.layer_minus,
                               den > 0.0 ? num / den : 0.0);

      // Split it. Each half alone represents the linear function n_shared on
      // all of Omega -- barycentric interpolation of a linear function is
      // exact -- so each half's node-weight integral must equal the analytic
      // integral of n_shared over Omega = [0,N]^8, which is N^9 / 2. If the
      // halves disagree with that, w and the basis are not indexed the same
      // way, and the master has been maximising something that is not the
      // integral of the represented function.
      Eigen::VectorXd d_plus = Eigen::VectorXd::Zero(layout.num_x);
      Eigen::VectorXd d_minus = Eigen::VectorXd::Zero(layout.num_x);
      for (int k = 0; k < layout.eta_s[static_cast<std::size_t>(
                              cand.layer_plus)]; ++k) {
        const int col = layout.idxX(cand.layer_plus, k);
        d_plus(col) = cand.d[static_cast<std::size_t>(col)];
      }
      for (int k = 0; k < layout.eta_s[static_cast<std::size_t>(
                              cand.layer_minus)]; ++k) {
        const int col = layout.idxX(cand.layer_minus, k);
        d_minus(col) = cand.d[static_cast<std::size_t>(col)];
      }
      const double analytic = std::pow(kN, 9) / 2.0;
      const double ip = w.head(layout.num_x).dot(d_plus);
      const double im = -w.head(layout.num_x).dot(d_minus);
      std::cout << std::format(
          "      integral of n_{} from plane {} = {:.6e}  (ratio to analytic "
          "{:.6f})\n"
          "      integral of n_{} from plane {} = {:.6e}  (ratio to analytic "
          "{:.6f})\n"
          "      analytic N^9/2 = {:.6e}\n",
          cand.shared_coord, cand.layer_plus, ip, ip / analytic,
          cand.shared_coord, cand.layer_minus, im, im / analytic, analytic);
    }
  }
}

// ---------------------------------------------------------------------------
// Hypothesis 1: one backward step multiplies the solution, and the effect is
// or is not a function of dt.
//
// Same geometry, same starting x_next, two models differing only in t_delta.
// Compared at equal physical time: step k at dt against step 2k at dt/2.
// ---------------------------------------------------------------------------
TEST(common, DISABLED_step_amplification) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  auto approx = makeApproximator(240);
  approx.getIntersectionPoints();
  approx.precomputeMatrices();

  constexpr int phase = 0;
  const auto& layout = approx.layouts()[phase];
  const auto& w = approx.nodeWeights()[phase];
  const br::ReducedLpInput& full = approx.reducedLpInput(phase);

  br::ReducedLpInput half = full;
  half.t_delta = full.t_delta / 2.0;

  std::cout << std::format("\ndt = {:.4f} and {:.4f}, num_x = {}\n",
                           full.t_delta, half.t_delta, layout.num_x);

  StepModel model_full = buildStepModel(full, layout);
  StepModel model_half = buildStepModel(half, layout);

  constexpr int kSteps = 8;
  std::vector<double> x_full(static_cast<std::size_t>(layout.num_x), 0.0);
  std::vector<double> x_half(static_cast<std::size_t>(layout.num_x), 0.0);

  std::cout << "\nшаг  физ.время      dt: max|x|      int V        "
               "dt/2: max|x|      int V       отношение int V\n";
  for (int k = 1; k <= kSteps; ++k) {
    x_full = stepOnce(model_full, full, layout.num_x, x_full);
    // Two half steps cover the same physical interval as one full step.
    x_half = stepOnce(model_half, half, layout.num_x, x_half);
    x_half = stepOnce(model_half, half, layout.num_x, x_half);

    const double iv_full = integralOfV(w, x_full);
    const double iv_half = integralOfV(w, x_half);
    std::cout << std::format(
        "{:>3}  {:>8.2f}   {:>12.4e} {:>12.4e}   {:>12.4e} {:>12.4e}   {:>8.3f}\n",
        k, k * full.t_delta, maxAbs(x_full), iv_full, maxAbs(x_half), iv_half,
        iv_full == 0.0 ? 0.0 : iv_half / iv_full);
  }
}

// ---------------------------------------------------------------------------
// Hypothesis 3: does the border LP's output stay under its candidates?
//
// The border condition is V_target <= max_rho Phi_rho pointwise, so the
// gauge-invariant integrals must satisfy
//   int V_target  <=  int max_rho Phi_rho  <=  |P| * max_rho int Phi_rho.
// A ratio far above the candidate count means the border LP is inflating the
// function, not just its coefficients.
// ---------------------------------------------------------------------------
TEST(common, DISABLED_border_vs_candidates) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  auto approx = makeApproximator(240);
  approx.getIntersectionPoints();
  approx.precomputeMatrices();

  constexpr int target_phase = 0;
  constexpr int source_phase = 1;
  const auto& tgt_layout = approx.layouts()[target_phase];
  const auto& src_layout = approx.layouts()[source_phase];
  const auto& w_tgt = approx.nodeWeights()[target_phase];
  const auto& w_src = approx.nodeWeights()[source_phase];

  // Candidates: the source phase marched back from zero, which is what switch
  // count 0 produces in the real run.
  const br::ReducedLpInput& src_input = approx.reducedLpInput(source_phase);
  StepModel src_model = buildStepModel(src_input, src_layout);
  std::vector<std::vector<double>> candidates;
  std::vector<double> x(static_cast<std::size_t>(src_layout.num_x), 0.0);
  for (int k = 0; k < 3; ++k) {
    x = stepOnce(src_model, src_input, src_layout.num_x, x);
    candidates.push_back(x);
  }

  double best_candidate_integral = 0.0;
  for (const auto& c : candidates) {
    best_candidate_integral
        = std::max(best_candidate_integral, std::abs(integralOfV(w_src, c)));
    std::cout << std::format("candidate: max|x| = {:.4e}, int V = {:.4e}\n",
                             maxAbs(c), integralOfV(w_src, c));
  }

  baa::CourierBorderOptions options;
  options.certificate_tol = 1e-2;
  options.max_certificate_cache = 512;
  baa::CourierBorderSolver solver(options);
  solver.prepare(approx.phaseGeometries(), approx.layouts(),
                 approx.nodeWeights(), kN);

  baa::CourierBorderRequest request;
  request.target_phase = target_phase;
  request.source_phase = source_phase;
  request.candidates = std::span<const std::vector<double>>(candidates);
  baa::CourierBorderStats stats;
  const std::vector<double> z
      = solver.solve(request, baa::ApproximationMode::Lower, &stats);

  // The integral comparison is inconclusive on its own: the candidates change
  // sign, so int max_rho Phi_rho can legitimately exceed max_rho int Phi_rho.
  // Check the border condition where it is actually enforced -- at the vertices
  // of the target regions, which is exactly where condition (a) is imposed.
  //
  // Phi_rho at an 8D point needs point location on the SOURCE partition, done
  // here by brute force over each source layer's triangles.
  const auto& src_geom = approx.phaseGeometries()[source_phase];
  const auto& tgt_geom = approx.phaseGeometries()[target_phase];
  auto evalSource = [&](const Eigen::VectorXd& n,
                        const std::vector<double>& x_src, double* out) {
    double acc = 0.0;
    for (int s = 0; s < baa::kSubsystemCount; ++s) {
      const auto& layer = src_geom.layers[static_cast<std::size_t>(s)];
      const Eigen::Vector2d y(n(layer.axes[0]), n(layer.axes[1]));
      bool found = false;
      for (std::size_t t = 0; t < layer.bases.size(); ++t) {
        const Eigen::Vector3d a = layer.bases[t].H * y + layer.bases[t].h;
        if (a(0) >= -1e-7 && a(1) >= -1e-7 && a(2) >= -1e-7) {
          for (int k = 0; k < 3; ++k) {
            acc += a(k) * x_src[static_cast<std::size_t>(src_layout.idxX(
                       s, layer.bases[t].vertex_ids[static_cast<std::size_t>(k)]))];
          }
          found = true;
          break;
        }
      }
      if (!found) {
        return false;
      }
    }
    *out = acc;
    return true;
  };

  double worst_violation = -std::numeric_limits<double>::infinity();
  double worst_at_V = 0.0;
  double worst_at_max = 0.0;
  long long checked = 0;
  long long unlocated = 0;
  const auto& blocks = tgt_geom.blocks;
  const int steps = 12;
  for (int ia = 0; ia < steps; ++ia) {
    for (int ib = 0; ib < steps; ++ib) {
      for (int ic = 0; ic < steps; ++ic) {
        const std::array<int, baa::kBlockCount> js = {
            ia * blocks[0].numRegions() / steps,
            ib * blocks[1].numRegions() / steps,
            ic * blocks[2].numRegions() / steps};
        // One vertex of each block cell -> one vertex of the product region.
        Eigen::VectorXd n = Eigen::VectorXd::Zero(baa::kSpaceDim);
        baa::SparseVec phi;
        bool ok = true;
        for (int b = 0; b < baa::kBlockCount; ++b) {
          const auto& bg = blocks[static_cast<std::size_t>(b)];
          const auto& verts
              = bg.vertices[static_cast<std::size_t>(js[static_cast<std::size_t>(b)])];
          if (verts.empty()) { ok = false; break; }
          const Eigen::VectorXd& v = verts[(ia + ib + ic) % verts.size()];
          for (int c = 0; c < bg.coord_count; ++c) {
            n(bg.coords[static_cast<std::size_t>(c)]) = v(c);
          }
          const baa::SparseVec row = baa::buildPhiRowBlock(
              tgt_geom, tgt_layout, b, js[static_cast<std::size_t>(b)], v,
              baa::kEps);
          for (std::size_t k = 0; k < row.cols.size(); ++k) {
            phi.add(row.cols[k], row.vals[k], baa::kGeomEps);
          }
        }
        if (!ok) { continue; }
        const double v_border = phi.dot(z);
        double best = -std::numeric_limits<double>::infinity();
        bool any = false;
        for (const auto& c : candidates) {
          double val = 0.0;
          if (evalSource(n, c, &val)) { best = std::max(best, val); any = true; }
        }
        if (!any) { ++unlocated; continue; }
        ++checked;
        if (v_border - best > worst_violation) {
          worst_violation = v_border - best;
          worst_at_V = v_border;
          worst_at_max = best;
        }
      }
    }
  }
  std::cout << std::format(
      "\nПоточечная проверка краевого условия V_border(nu) <= max_rho Phi_rho(nu)\n"
      "  проверено вершин: {} (не локализовано: {})\n"
      "  худшее нарушение V - max_rho Phi = {:.6e}\n"
      "  в этой точке: V = {:.6e}, max_rho Phi = {:.6e}\n",
      checked, unlocated, worst_violation, worst_at_V, worst_at_max);

  const double border_integral = std::abs(integralOfV(w_tgt, z));
  std::cout << std::format(
      "\nborder output: max|x| = {:.4e}, int V = {:.4e}\n"
      "iterations = {}, cuts = {}, worst zeta = {:.3e}, hit_box = {}\n"
      "best candidate int V = {:.4e}\n"
      "ratio int V(border) / max_rho int Phi_rho = {:.3f}  "
      "(bound if the border condition holds: {} candidates)\n",
      maxAbs(z), border_integral, stats.iterations, stats.cuts_added,
      stats.worst_zeta, stats.master_hit_box, best_candidate_integral,
      best_candidate_integral == 0.0 ? 0.0
                                     : border_integral / best_candidate_integral,
      candidates.size());
}
