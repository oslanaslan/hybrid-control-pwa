#ifndef HCPWA_COURIER_BORDER_SOLVER_HPP
#define HCPWA_COURIER_BORDER_SOLVER_HPP

#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <memory>
#include <span>
#include <vector>

#include <spdlog/spdlog.h>

#include "barycentric_geometry_types.hpp"

// NOLINTBEGIN(readability-identifier-naming)

namespace barycentric_affine_approximator {

namespace detail {

struct Point2 {
  double x = 0.0;
  double y = 0.0;
};

// Vertices of {triangle} cap {[lo_x, hi_x] x [lo_y, hi_y]}, deduplicated
// within kGeomEps. Empty when the two do not meet.
//
// This is what turns the continuous condition (b) into a finite list of LP
// rows: phi^(s) is affine on each piece {source triangle} cap {rectangle}, so
// requiring the courier to stay under it on a piece is the same as requiring it
// at the piece's vertices.
//
// Exposed for testing. Points exactly on a rectangle edge count as inside, and
// slivers are kept rather than discarded, because dropping a vertex here drops
// an LP row and can break soundness.
std::vector<Point2> clipTriangleToRect(const std::array<Point2, 3>& tri,
                                       double lo_x, double hi_x, double lo_y,
                                       double hi_y);

}  // namespace detail

struct CourierBorderOptions {
  // Benders outer iterations. Hitting the cap is a hard failure: the returned
  // node values would not be certified on every region and therefore would not
  // be a bound at all.
  int max_iterations = 60;

  // Regions whose cut is appended per outer iteration, most violated first.
  // Bounds the per-iteration master growth; Benders does not need every cut at
  // once. The convergence test still sweeps all regions.
  int max_cuts_per_iteration = 512;

  // Certification tolerance on the phase-I objective zeta*. Must sit above the
  // HiGHS feasibility tolerance, otherwise certified regions would re-report as
  // violated forever.
  //
  // It used to be 1e-6, which is exactly kHighsSolutionTol and so did not
  // satisfy its own requirement: a zeta of 1e-6 is indistinguishable from zero
  // to the solver that produced it, the cut built for such a region is not
  // obliged to separate anything, and the separation guard below rejected it.
  // Two orders of margin, matching what kDualIdentityTol needs for the same
  // reason.
  double certificate_tol = 1e-4;

  // Rows of the seeded master relaxation. 0 disables seeding.
  int max_seed_rows = 8192;

  // Couriers kept as screens during a region sweep. A courier that certifies
  // one region very often certifies its neighbours, and re-testing a cached
  // one costs eight lookups against a subproblem solve. 0 disables screening.
  //
  // This is what makes the sweep affordable at all: every region still has to
  // be certified, but only the ones no cached courier covers pay for an LP.
  // Skipping is sound by construction -- a region is skipped only when a
  // concrete feasible courier for it has been exhibited.
  int max_certificate_cache = 64;

  // Safety factor on the data-derived master column box.
  double master_box_scale = 8.0;

  // Preparation memory guard, bytes. prepare() estimates the clipped-vertex
  // pool up front and throws with the estimate rather than exhausting RAM.
  std::size_t max_prepare_bytes = 6UL << 30;

  bool highs_verbose = false;
};

// Per-call diagnostics, so getBorderConditions can log one line without the
// solver knowing anything about the caller's context.
struct CourierBorderStats {
  int iterations = 0;
  long long cuts_added = 0;
  long long subproblems_solved = 0;
  double worst_zeta = 0.0;
  double master_objective = 0.0;
  bool master_hit_box = false;
};

// One call of the border problem. Everything that changes between calls lives
// here; everything that does not lives in the prepared solver.
//
//   target_phase  phase whose node values are being solved for
//   source_phase  phase the candidates Phi_rho are expressed in
//   candidates    x_src(rho), each of size layouts[source_phase].num_x
//
// The span must outlive the solve() call; the solver never copies it.
struct CourierBorderRequest {
  int target_phase = 0;
  int source_phase = 0;
  std::span<const std::vector<double>> candidates;
  // Logging context only; does not affect the result.
  int theta_idx = 0;
  int switch_cnt = 0;
};

// Border condition solver based on courier functions and Benders decomposition
// (companion paper, section "Метод функций-курьеров").
//
// Lower mode certifies  Vt(n) <= max_rho Phi_rho(n)  on all of Omega by
// exhibiting, for every target region j, an affine courier l_j with
//   (a) Vt(nu) <= l_j(nu) at every vertex nu of Omega_target^(j), and
//   (b) l_j <= Phi_{rho(j)} on Omega_target^(j),
// certified plane by plane on the rectangles of the projection lemma. Because
// Vt - l_j is affine on the region, (a) at the vertices gives it on the whole
// region, and (b) then yields Vt <= Phi_{rho(j)} <= max_rho Phi_rho.
//
// Upper mode flips both inequalities and replaces Phi_{rho(j)} by the pointwise
// max over rho, which needs no per-region candidate choice: requiring
// l^(s) >= max_rho phi_rho^(s) on every plane already implies
// sum_s l^(s) >= Phi_rho for every rho.
//
// The set of vertices of Omega_1^(j1) cap Omega_2^(j2) is never built. The
// source partition enters only through evaluating Phi_rho at points and through
// 2D clipping against the projection rectangles.
//
// prepare() is called once, before any worker thread starts. solve() is const
// and reentrant: it allocates its own Highs instances and mutates no solver
// state.
class CourierBorderSolver {
 public:
  explicit CourierBorderSolver(CourierBorderOptions options = {});

  // Geometry ingestion. The referenced objects must outlive the solver and must
  // not be mutated afterwards; only pointers are retained, because
  // region_vertices alone is gigabytes.
  //
  // Call after BarycentricAffineApproximator::getIntersectionPoints(), which is
  // what fills node_weights.
  void prepare(const std::array<PhaseGeometry, kPhases>& geometries,
               const std::array<BarycentricVarLayout, kPhases>& layouts,
               const std::array<Eigen::VectorXd, kPhases>& node_weights,
               double n_max);

  bool prepared() const { return prepared_; }

  // Returns node values of size layouts[target_phase].num_x, certified on every
  // target region. Throws if certification cannot be completed.
  std::vector<double> solve(const CourierBorderRequest& request,
                            ApproximationMode mode,
                            CourierBorderStats* stats = nullptr) const;

  // Independent re-verification of a candidate solution: walks every region,
  // rebuilds its courier from scratch and returns the worst violation, which is
  // <= 0 exactly when z is certified everywhere. Far too slow for production;
  // used by the tests and by an opt-in validation run.
  double worstCertificateResidual(const CourierBorderRequest& request,
                                  ApproximationMode mode,
                                  const std::vector<double>& z) const;

  // Exposed only so the file-local subproblem builder can name it; the type
  // itself stays private to the .cpp.
  struct Impl;

 private:
  // shared_ptr, not unique_ptr, and deliberately so: unique_ptr to an
  // incomplete type forces a user-declared destructor, which would delete the
  // implicit move constructor of this class and, through it, of
  // BarycentricAffineApproximator -- which is returned by value in
  // test/common/barycentric_affine_solver.cpp. shared_ptr type-erases its
  // deleter at construction, so no destructor declaration is needed.
  std::shared_ptr<Impl> impl_;
  CourierBorderOptions options_;
  bool prepared_ = false;
  std::shared_ptr<spdlog::logger> logger_;
};

}  // namespace barycentric_affine_approximator

// NOLINTEND(readability-identifier-naming)

#endif  // HCPWA_COURIER_BORDER_SOLVER_HPP
