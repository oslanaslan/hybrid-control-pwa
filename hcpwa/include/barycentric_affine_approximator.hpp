#ifndef HCPWA_BARYCENTRIC_AFFINE_APPROXIMATOR_HPP
#define HCPWA_BARYCENTRIC_AFFINE_APPROXIMATOR_HPP

#include <Eigen/Core>
#include <algorithm>
#include <stdexcept>
#include <Eigen/Dense>
#include <Highs.h>
#include <algo.hpp>
#include <array>
#include <limits>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <spdlog/spdlog.h>

#include "barycentric_geometry_types.hpp"
#include "courier_border_solver.hpp"
#include "util/band_lp.hpp"
#include "util/block_reduction_lp.hpp"
#include "util/gauge_fix.hpp"
#include "interval_building.hpp"
#include "util/value_function_utils.hpp"

// Keep the same public naming style as the existing GlobalAffineApproximator so
// this new implementation is easy to compare against the production template.
// NOLINTBEGIN(readability-identifier-naming)

namespace barycentric_affine_approximator {

// Axis ids in the paper are one-based. All constants in this implementation are
// zero-based because Eigen vectors and the rest of the C++ code are zero-based.
extern const std::vector<int> kInIds;
extern const std::vector<int> kOutIds;

struct SystemParams {
  double N;
  double F;
  double v;
  double w;
  // beta parameters, with axis ids in comments matching the paper.
  double b51;
  double b57;
  double b84;
  double b86;
  double b31;
  double b36;
  double b24;
  double b27;
  // incoming-flow lower bounds.
  double f2min;
  double f3min;
  double f5min;
  double f8min;
  // incoming-flow upper bounds.
  double f2max;
  double f3max;
  double f5max;
  double f8max;
};

// CTM data of one region restricted to a set of cells, in the order the cells
// were given:
//   drift(n) = a n + f,   g(n) = g_vec^T n + g_scal.
struct CtmRegionData {
  Eigen::MatrixXd a;
  Eigen::VectorXd f;
  Eigen::VectorXd g_vec;
  double g_scal = 0.0;
};

// The disturbance box on a set of cells:
//   c(n) = qc_diag .* n + qc_off,   rho(n) = qr_diag .* n + qr_off.
// Q_c and Q_r are diagonal by construction -- every entry depends on its own
// cell only -- so only the diagonals are kept.
struct BoxRegionData {
  Eigen::VectorXd qc_diag;
  Eigen::VectorXd qc_off;
  Eigen::VectorXd qr_diag;
  Eigen::VectorXd qr_off;
};

// Everything the reduced LP needs about one block-region.
struct BlockSystemMatrices {
  CtmRegionData ctm;
  BoxRegionData box;
};

// How the band LP of one node is solved. The first attempt is `solver`; if it
// does not return an optimal point, the other solver is tried on a fresh
// model, and only then is a primal-feasible non-optimal point accepted -- the
// same tightness-only relaxation the one-step solver used to apply, with the
// residuals re-derived from the formulas afterwards either way. Anything else
// is fatal: the node has no certified point.
struct BandLpOptions {
  // "ipm" (interior point, crossover per run_crossover) or "simplex" (dual).
  // Measured on the N=100 arrangement with nine stages (59 166 rows x 12 852
  // columns): interior point 6.8 s, interior point with crossover 12.3 s, dual
  // simplex 75 s, all three optimal with the same integral. The interior point
  // without crossover returns a strictly feasible point, whose residuals sit
  // below zero rather than at the solver tolerance.
  std::string solver = "ipm";
  bool run_crossover = false;
  int ipm_iteration_limit = 1000;
  // A cycling simplex has to end. One step-4 solve of the one-step LP ran 3.8
  // million iterations with the objective oscillating between exactly two
  // values; the worst honest solve in the same run took 78403, so this cap is
  // five times that.
  int simplex_iteration_limit = 400000;
  double time_limit = std::numeric_limits<double>::infinity();
};

// What one band solve did, for the run log.
struct BandSolveStats {
  int stages = 0;
  int rows = 0;
  int cols = 0;
  long long nnz = 0;
  // The attempt that produced the point: "ipm" or "simplex".
  std::string solver;
  // HighsModelStatus of that attempt, as an int.
  int model_status = 0;
  bool accepted_non_optimal = false;
  long long iterations = 0;
  double seconds = 0.0;
  // s * sum_{k<L} omega_k q^T z_k in the original units.
  double objective = 0.0;
  // dt * sum_{k<=L} omega_k q^T z_k: the integral of the estimate over
  // [t_0, theta] x Omega.
  double integral = 0.0;
  // Largest s * F over all stages, both ends, and its admitted limit.
  double worst_residual = 0.0;
  double residual_limit = 0.0;
};

class BarycentricAffineApproximator {
 private:
  double t_max_;
  int t_split_count_;
  int max_switches_;
  double tau_min_;
  double tau_max_;
  double t_delta_;
  SystemParams system_params_;
  std::vector<double> t_range_;
  std::vector<int> t_index_;

  bool highs_verbose_ = false;
  ApproximationMode approximation_mode_ = ApproximationMode::Upper;
  // Opt-in re-derivation of the courier border conditions. The band LP
  // residual check is not behind this: separability made it cheap enough to
  // run unconditionally, but the courier re-derivation still costs a pass over
  // all M regions.
  bool validate_ = false;

  BandLpOptions band_options_;

  hcpwa::TriangleGeometryOptions geometry_options_;

  ValueFunction value_function_;
  std::vector<Eigen::VectorXd> cube_angle_vertices_;

  std::array<PhaseGeometry, kPhases> phase_geometries_;
  std::array<BarycentricVarLayout, kPhases> layouts_;

  // q = integral over Omega of the hat function of each node (lemma on
  // quadrature). The objective of the band LP and of the courier master.
  // Computed in closed form from triangle areas.
  std::array<Eigen::VectorXd, kPhases> node_weights_;

  // Per phase, per block, per block-region.
  std::array<std::array<std::vector<BlockSystemMatrices>, kBlockCount>, kPhases>
      block_system_;

  // The block data of one stage, one per phase. reduced_inputs_ is what the
  // band assembler, the terminal RHS and the exact worst residual all read, so
  // the three cannot disagree about the row scale.
  std::array<block_reduction::ReducedLpInput, kPhases> reduced_inputs_;
  std::array<block_reduction::ReducedLpColLayout, kPhases> reduced_cols_;
  std::array<block_reduction::ReducedLpRowLayout, kPhases> reduced_rows_;

  // The gauge normalisation of each phase, shared by the band LP and the
  // courier master. Built and verified in precomputeMatrices().
  std::array<GaugeFix, kPhases> gauge_fixes_;

  std::shared_ptr<spdlog::logger> logger_;

  interval_building::ThetaTIndexLists theta_t_index_lists_;
  // Nodes the theta lists carry that no earlier layer can reach; see the
  // constructor and getBorderConditions.
  int unreachable_nodes_ = 0;
  int checked_nodes_ = 0;

  // Border conditions at switching instants. Held by value: CourierBorderSolver
  // keeps its prepared tables behind a shared_ptr precisely so that this class
  // stays implicitly movable.
  CourierBorderOptions courier_options_;
  CourierBorderSolver courier_solver_;

  SparseVec buildPhiRow(int phase, int region, const Eigen::VectorXd& point,
                        double tolerance = kEps) const;

  // buildPhiRow for one coordinate block. point carries only that block's own
  // coordinates, in the block's order, and the columns produced are still the
  // global x columns -- the supports of the three blocks are disjoint, so the
  // three block rows add up to the full phi row of the product region.
  SparseVec buildPhiRowBlock(int phase, int block, int block_region,
                             const Eigen::VectorXd& point,
                             double tolerance = kEps) const;

  // s of step 2.1: the single parameter separating the two approximation
  // directions after the corrections. Every row sign and the objective sign are
  // derived from it.
  double signS() const {
    return approximation_mode_ == ApproximationMode::Upper ? 1.0 : -1.0;
  }

  // Fills node_weights_ from the triangle areas of each projection layer.
  // Called at the end of getIntersectionPoints(), once layouts_ are known.
  void computeNodeWeights();

 public:
  BarycentricAffineApproximator(double t_max, int t_split_count,
                                double tau_min, double tau_max,
                                const SystemParams& system_params,
                                bool highs_verbose = false,
                                ApproximationMode mode
                                = ApproximationMode::Upper);

  ApproximationMode approximationMode() const { return approximation_mode_; }

  // Enables the independent re-derivation of the courier border conditions.
  // The band LP residual check is unconditional and not affected by this.
  void setValidate(bool validate) { validate_ = validate; }

  void setBandLpOptions(const BandLpOptions& options) {
    band_options_ = options;
  }

  // Stops the recursion early, for staged debugging: level r reads only levels
  // below it, so a run capped at r = 2 exercises every part of the machinery
  // -- terminal family, courier border conditions, the band LP with a nonzero
  // terminal value -- at a fraction of the cost. Can only lower the count the
  // constructor derived from the horizon, never raise it.
  void setMaxSwitches(int max_switches) {
    if (max_switches < 0) {
      throw std::invalid_argument("setMaxSwitches: negative level count");
    }
    max_switches_ = std::min(max_switches_, max_switches);
  }
  int maxSwitches() const { return max_switches_; }
  const BandLpOptions& bandLpOptions() const { return band_options_; }

  // Must be set before getIntersectionPoints(). Turning build_8d_vertices off
  // skips the 8D product of the block cells, which nothing but the courier
  // solver still reads.
  void setGeometryOptions(const hcpwa::TriangleGeometryOptions& options) {
    geometry_options_ = options;
  }

  // Must be called before run(); prepare() reads these when it builds the
  // courier tables.
  void setCourierOptions(const CourierBorderOptions& options) {
    courier_options_ = options;
  }

  // The grid Theta of admissible next switching instants after theta: the
  // time-grid points inside [theta + tau_min, theta + tau_max] intersected
  // with [0, T]. Public because the containment of Theta in that window is
  // what makes the transfer step a bound, and a test pins it.
  std::vector<int> admissibleThetaIds(double theta) const;

  double getBetaParamForAxis(int i, int j) const;
  std::pair<double, double> getFMinMaxForAxis(int i) const;

  void getIntersectionPoints();

  std::pair<Eigen::RowVectorXd, Eigen::RowVectorXd> getFIJMinResolution(
      int i, int j, const Eigen::VectorXd& n) const;

  Eigen::VectorXd areaCentroidCoords(int j, int phase) const;

  // The CTM drift and the disturbance box of one region, restricted to a set
  // of cells. Both the full 8D path and the per-block path go through these,
  // so the two cannot drift. n must carry the coordinates of every cell in
  // `cells`; the rest may be anything, including NaN.
  CtmRegionData ctmDataForCells(int phase, const std::vector<int>& cells,
                                const Eigen::VectorXd& n) const;
  BoxRegionData boxDataForCells(const std::vector<int>& cells,
                                const Eigen::VectorXd& n0) const;

  // Centroid of one block-region, in that block's own coordinate order. The
  // region is a product, so concatenating the three block centroids gives the
  // 8D centroid exactly -- branch resolution is the same either way.
  Eigen::VectorXd blockCentroidCoords(int phase, int block,
                                      int block_region) const;

  BlockSystemMatrices getBlockSystemMatrices(int phase, int block,
                                             int block_region) const;

  std::vector<BlockSystemMatrices> precomputeBlockSystemMatrices(int phase,
                                                                 int block);

  // Cross-checks the per-block CTM and box data against the same assemblers run
  // on all eight cells at once. getFIJMinResolution returns rows that depend on
  // n only through which branch is the minimum, so agreement is exact and any
  // difference means the two paths resolved different branches.
  void validateBlockDecomposition(int phase) const;

  // Builds reduced_inputs_[phase] from the block geometry and the block system
  // matrices. Called once, from precomputeMatrices().
  void buildReducedLpInput(int phase);

  std::vector<double> getBorderConditions(int switch_phase, int theta_idx,
                                          double theta, int switch_cnt) const;

  // The band LP of one node: z_L = z_terminal given, num_stages segments of
  // length t_delta before it, the gauge pins of level switch_cnt applied at
  // every stage. Returns z_0, ..., z_{L-1}, z_L. Throws unless every stage is
  // certified by the exact worst residual -- that check is the statement that
  // the result is a bound and does not depend on anything HiGHS asserts.
  // Reentrant: allocates its own solver and mutates no member.
  std::vector<std::vector<double>> solveBandLp(
      int phase, int switch_cnt, int num_stages,
      const std::vector<double>& z_terminal,
      BandSolveStats* stats = nullptr) const;

  // The exact worst residual of every stage of z = (z_0, ..., z_L) over the
  // whole product of regions and vertices, checked against the admitted limit.
  // Returns the worst s * F; throws if any stage has the wrong sign.
  double validateBandResiduals(int phase,
                               const std::vector<std::vector<double>>& z) const;

  void precomputeMatrices();

  // The block data of one stage of one phase, for diagnostics and for driving
  // a band solve from outside run().
  const block_reduction::ReducedLpInput& reducedLpInput(int phase) const {
    return reduced_inputs_[static_cast<std::size_t>(phase)];
  }
  const block_reduction::ReducedLpRowLayout& reducedLpRows(int phase) const {
    return reduced_rows_[static_cast<std::size_t>(phase)];
  }
  const block_reduction::ReducedLpColLayout& reducedLpCols(int phase) const {
    return reduced_cols_[static_cast<std::size_t>(phase)];
  }
  const std::array<GaugeFix, kPhases>& gaugeFixes() const {
    return gauge_fixes_;
  }

  // Read-only views of the ingested geometry, for diagnostics and for driving
  // the border solver outside run().
  const std::array<PhaseGeometry, kPhases>& phaseGeometries() const {
    return phase_geometries_;
  }
  const std::array<BarycentricVarLayout, kPhases>& layouts() const {
    return layouts_;
  }
  const std::array<Eigen::VectorXd, kPhases>& nodeWeights() const {
    return node_weights_;
  }
  double tDelta() const { return t_delta_; }

  void run(const std::string& output_folder_path, int n_threads = 2);

  void dumpInitParamsToJson(const std::string& filepath) const;
};

}  // namespace barycentric_affine_approximator

// NOLINTEND(readability-identifier-naming)

#endif  // HCPWA_BARYCENTRIC_AFFINE_APPROXIMATOR_HPP
