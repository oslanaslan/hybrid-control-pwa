#ifndef HCPWA_BARYCENTRIC_AFFINE_APPROXIMATOR_HPP
#define HCPWA_BARYCENTRIC_AFFINE_APPROXIMATOR_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Highs.h>
#include <algo.hpp>
#include <array>
#include <memory>
#include <mutex>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <spdlog/spdlog.h>

#include "barycentric_geometry_types.hpp"
#include "courier_border_solver.hpp"
#include "util/block_reduction_lp.hpp"
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
  // Opt-in re-derivation of the courier border conditions. The main LP residual
  // check is no longer behind this: separability made it cheap enough to run
  // unconditionally, but the courier re-derivation still costs a pass over all
  // M regions.
  bool validate_ = false;

  // Objective weights of the reduced LP. The default reproduces the product
  // objective exactly and keeps the LP bounded; see ObjectiveWeights.
  block_reduction::ObjectiveWeights objective_weights_
      = block_reduction::ObjectiveWeights::ProductCount;
  // Weak epsilon * ||z||_1 regularization for a reproducible tie-break.
  double tie_break_eps_ = 0.0;

  hcpwa::TriangleGeometryOptions geometry_options_;

  ValueFunction value_function_;
  std::vector<Eigen::VectorXd> cube_angle_vertices_;

  std::array<PhaseGeometry, kPhases> phase_geometries_;
  std::array<BarycentricVarLayout, kPhases> layouts_;

  // w = integral over Omega of phi^(phase)(n) dn, the objective of the border LP
  // (step 2.2, section 7). Computed in closed form from triangle areas.
  std::array<Eigen::VectorXd, kPhases> node_weights_;

  // Per phase, per block, per block-region.
  std::array<std::array<std::vector<BlockSystemMatrices>, kBlockCount>, kPhases>
      block_system_;

  // The reduced LP, one per phase. reduced_inputs_ is the block data the row
  // builder, the per-step RHS update and the exact worst residual all read, so
  // the three cannot disagree about the row scale.
  std::array<block_reduction::ReducedLpInput, kPhases> reduced_inputs_;
  std::array<block_reduction::ReducedLpColLayout, kPhases> reduced_cols_;
  std::array<block_reduction::ReducedLpRowLayout, kPhases> reduced_rows_;

  std::shared_ptr<spdlog::logger> logger_;
  std::vector<std::unique_ptr<Highs>> highs_solvers_;
  std::vector<std::vector<double>> row_lowers_;
  std::vector<std::vector<double>> row_uppers_;
  std::vector<std::unique_ptr<std::mutex>> solver_mutexes_;

  interval_building::ThetaTIndexLists theta_t_index_lists_;

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

  std::vector<int> locateRegions(int phase, const Eigen::VectorXd& point,
                                 double tolerance = kEps) const;

  double evaluateBarycentricValue(int phase, const std::vector<double>& x,
                                  const Eigen::VectorXd& point,
                                  double tolerance = kEps) const;

  std::vector<int> admissibleThetaIds(double theta) const;

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
  // The main LP residual check is unconditional and not affected by this.
  void setValidate(bool validate) { validate_ = validate; }

  void setObjectiveWeights(block_reduction::ObjectiveWeights weights) {
    objective_weights_ = weights;
  }
  void setTieBreakEps(double eps) { tie_break_eps_ = eps; }

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

  double getBetaParamForAxis(int i, int j) const;
  std::pair<double, double> getFMinMaxForAxis(int i) const;

  void getIntersectionPoints();

  std::pair<Eigen::RowVectorXd, Eigen::RowVectorXd> getFIJMinResolution(
      int i, int j, const Eigen::VectorXd& n) const;

  Eigen::VectorXd areaCentroidCoords(int j, int phase) const;

  std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd, double>
  getAMatrFVecGVecAndGScalJ(int j, int phase) const;

  std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd,
             Eigen::VectorXd>
  getQQForArea(int j, int phase) const;

  // The two assemblers above, restricted to a set of cells. Both the full 8D
  // path and the per-block path go through these, so the two cannot drift.
  // n must carry the coordinates of every cell in `cells`; the rest may be
  // anything, including NaN.
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

  // Assembles the reduced LP of one phase. Pure function of the precomputed
  // block data: precomputeMatrices() fills reduced_inputs_ single-threaded,
  // before any solver exists.
  block_reduction::ReducedLpMatrices prepareLpMatrices(int phase) const;

  // Builds reduced_inputs_[phase] from the block geometry and the block system
  // matrices. Called once, from precomputeMatrices().
  void buildReducedLpInput(int phase);

  std::vector<double> getBorderConditions(int switch_phase, int theta_idx,
                                          double theta, int switch_cnt) const;

  std::tuple<std::unique_ptr<Highs>, std::vector<double>, std::vector<double>>
  initializeHighs(int phase);

  void updateHighsRhsUpperBounds(int phase, int solver_index,
                                 const std::vector<double>& x_next);

  std::vector<double> solveLp(int phase, int solver_index);

  std::vector<double> solveMainLpStep(int phase, int solver_index,
                                      const std::vector<double>& x_next);

  void precomputeMatrices();

  // Recomputes the exact worst residual of one segment over the whole product
  // of regions and vertices and checks its sign. Separability turns what used
  // to be a pass over prod_b R_b pairs into sum_b R_b, so this runs on every
  // step: it is the check that the constructed function is a bound.
  void validateStepResiduals(int phase, const std::vector<double>& x_next,
                             const std::vector<double>& z) const;

  // The assembled reduced LP of one phase, for diagnostics and for driving a
  // single step from outside run().
  const block_reduction::ReducedLpInput& reducedLpInput(int phase) const {
    return reduced_inputs_[static_cast<std::size_t>(phase)];
  }
  const block_reduction::ReducedLpRowLayout& reducedLpRows(int phase) const {
    return reduced_rows_[static_cast<std::size_t>(phase)];
  }
  const block_reduction::ReducedLpColLayout& reducedLpCols(int phase) const {
    return reduced_cols_[static_cast<std::size_t>(phase)];
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

  void run(const std::string& output_folder_path, int n_threads = 2);

  void dumpInitParamsToJson(const std::string& filepath) const;
};

}  // namespace barycentric_affine_approximator

// NOLINTEND(readability-identifier-naming)

#endif  // HCPWA_BARYCENTRIC_AFFINE_APPROXIMATOR_HPP
