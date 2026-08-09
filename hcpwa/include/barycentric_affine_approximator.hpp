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

constexpr int kPhases = 2;
constexpr int kSpaceDim = 8;
constexpr int kSubsystemCount = 5;

// Geometry and LP tolerances are intentionally separated by name in comments in
// the source, but kept numerically equal for this first implementation.
constexpr double kEps = 1e-5;
constexpr double kGeomEps = 1e-8;

// Tolerance for the optional end-to-end residual check. It has to sit above the
// HiGHS feasibility tolerance, otherwise the check would fire on solutions that
// the solver legitimately reports as optimal.
constexpr double kResidualValidationTol = 1e-4;

using ValueFunction = hcpwa::util::ValueFunction;

// Direction of the approximation. It enters every LP row and the objective only
// through the sign s of step 2.1:
//   s = +1 (Upper): the residual must satisfy max_omega F <= 0;
//   s = -1 (Lower): the residual must satisfy min_omega F >= 0.
// The name deliberately mirrors global_affine_approximator::ApproximationMode.
// The two namespaces are independent and no translation unit includes both
// headers, so the repeated name cannot create an ambiguity.
enum class ApproximationMode { Upper, Lower };

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

// Sparse vector used for formula-level LP coefficients such as phi_{j,nu} and
// a_{j,nu}. The representation is deliberately simple: the barycentric LP is
// difficult to debug, so we keep sparse operations explicit and local.
struct SparseVec {
  std::vector<int> cols;
  std::vector<double> vals;

  void add(int col, double value, double eps = kEps);
  double dot(const std::vector<double>& x) const;
};

// Sparse representation of Psi_j in the formula grad V = Psi_j x. There are
// exactly 8 rows because the state dimension is fixed at m = 8 in the paper and
// in this codebase.
struct SparsePsi {
  std::array<SparseVec, kSpaceDim> rows;
};

// Which endpoint of the time segment [t_{k-1}, t_k] a residual row controls.
// Both endpoints are required: the slope is constant on the segment, but the
// value x(t) varies, and the residual depends on the value through
// q_j = Psi_j x(t). Controlling one endpoint only leaves the condition violated
// on the rest of the segment (step 2.1, section 2.2).
enum class RhsKind {
  // Endpoint with the unknown value z = x_{k-1}. The modulus |Psi_j z| is taken
  // at the unknown, hence the y-lift.
  Left,
  // Endpoint with the known value x_next. The modulus is a number there, so the
  // row carries no y block and the whole RHS is arithmetic.
  Right,
};

// One residual row together with everything needed to recompute its RHS after
// x_next changes. Left rows need only phi; Right rows also need the per-(j,nu)
// data entering beta = q^T m + s rho^T |q| + g with q = Psi_j x_next.
struct ResidualRhsTerm {
  int row_id = 0;
  RhsKind kind = RhsKind::Left;
  SparseVec phi;

  // Used by RhsKind::Right only.
  int region = -1;
  Eigen::VectorXd m;    // A_j nu + f_j + c_j(nu)
  Eigen::VectorXd rho;  // radius of the disturbance box at nu
  double g = 0.0;       // g_i(nu)
};

// Barycentric coordinates on one 2D simplex:
//   alpha(z) = H z + h,
// where z is the projected 2D point P_s n and alpha has three components.
struct TriangleBasis {
  std::array<int, 3> vertex_ids{};
  Eigen::Matrix<double, 3, 2> H = Eigen::Matrix<double, 3, 2>::Zero();
  Eigen::Vector3d h = Eigen::Vector3d::Zero();
};

// One projection layer s. It owns the unique 2D vertices for that layer and the
// barycentric map for every triangle in that layer.
struct ProjectionLayer {
  std::array<int, 2> axes{};
  std::vector<hcpwa::TriangleWithUniqueVertices> triangles;
  std::vector<Eigen::Vector2d> unique_vertices;
  std::vector<TriangleBasis> bases;
};

// Geometry for one phase. A full 8D region stores both its vertices and the
// five triangle ids used to select the local barycentric charts.
struct PhaseGeometry {
  std::array<ProjectionLayer, kSubsystemCount> layers;
  std::vector<std::vector<Eigen::VectorXd>> region_vertices;
  std::vector<std::array<int, kSubsystemCount>> region_triangle_ids;
};

// One nonempty cell of the common refinement Omega_0^(j0) cap Omega_1^(j1).
// The upper border LP needs only the vertex set, but the lower one selects one
// family member per cell (step 2.2, section 6.1), so cell membership must
// survive the vertex deduplication.
struct RefinementCell {
  int phase0_area_id = 0;
  int phase1_area_id = 0;
  // Indices into common_refinement_vertices_.
  std::vector<int> vertex_ids;
};

// Layout of the LP variable vector:
//   Y = [x; y_0; y_1; ...; y_{M-1}],
// where x contains all unique 2D barycentric vertex values and y_j is the
// 8-vector auxiliary for |Psi_j x| in full region j.
struct BarycentricVarLayout {
  std::array<int, kSubsystemCount> eta_s{};
  std::array<int, kSubsystemCount> offset_s{};
  int num_x = 0;
  int num_regions = 0;
  int num_cols = 0;

  int idxX(int subsystem, int vertex_id) const;
  int idxY(int region, int dim) const;
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
  // Off by default: the check walks every (region, vertex) pair on every time
  // step, which is far too expensive for a production run.
  bool validate_ = false;

  ValueFunction value_function_;
  std::vector<Eigen::VectorXd> cube_angle_vertices_;
  std::vector<Eigen::VectorXd> common_refinement_vertices_;
  std::vector<RefinementCell> refinement_cells_;

  std::array<PhaseGeometry, kPhases> phase_geometries_;
  std::array<BarycentricVarLayout, kPhases> layouts_;

  // w = integral over Omega of phi^(phase)(n) dn, the objective of the border LP
  // (step 2.2, section 7). Computed in closed form from triangle areas.
  std::array<Eigen::VectorXd, kPhases> node_weights_;

  // phi^(phase)(g) for every g in common_refinement_vertices_. The border LP is
  // solved once per (level, theta, phase), so locating every g and rebuilding
  // its phi row on each call would repeat the same work thousands of times.
  std::array<std::vector<SparseVec>, kPhases> phi_at_refinement_;

  std::vector<std::vector<Eigen::MatrixXd>> A_j_matrs_;
  std::vector<std::vector<Eigen::VectorXd>> f_j_vecs_;
  std::vector<std::vector<Eigen::MatrixXd>> Q_c_j_matrs_;
  std::vector<std::vector<Eigen::VectorXd>> q_c_j_vecs_;
  std::vector<std::vector<Eigen::MatrixXd>> Q_r_j_matrs_;
  std::vector<std::vector<Eigen::VectorXd>> q_r_j_vecs_;
  std::vector<std::vector<Eigen::VectorXd>> g_j_vecs_;
  std::vector<std::vector<double>> g_j_scals_;

  std::array<std::vector<ResidualRhsTerm>, kPhases> rhs_terms_;
  std::array<std::vector<SparsePsi>, kPhases> psi_by_region_;

  std::shared_ptr<spdlog::logger> logger_;
  std::vector<std::unique_ptr<Highs>> highs_solvers_;
  std::vector<std::vector<double>> row_lowers_;
  std::vector<std::vector<double>> row_uppers_;
  std::vector<std::unique_ptr<std::mutex>> solver_mutexes_;

  interval_building::ThetaTIndexLists theta_t_index_lists_;

  SparseVec buildPhiRow(int phase, int region, const Eigen::VectorXd& point,
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

  // Recomputes both endpoint residuals of one segment directly from the
  // formulas and checks their sign. Guarded by validate_ because it costs a full
  // pass over all (region, vertex) pairs.
  void validateStepResiduals(int phase, const std::vector<double>& x_next,
                             const std::vector<double>& z) const;

 public:
  BarycentricAffineApproximator(double t_max, int t_split_count,
                                double tau_min, double tau_max,
                                const SystemParams& system_params,
                                bool highs_verbose = false,
                                ApproximationMode mode
                                = ApproximationMode::Upper);

  ApproximationMode approximationMode() const { return approximation_mode_; }

  // Enables the per-step residual recomputation of validateStepResiduals().
  // Intended for debugging runs on a reduced geometry.
  void setValidate(bool validate) { validate_ = validate; }

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

  std::tuple<std::vector<Eigen::MatrixXd>, std::vector<Eigen::VectorXd>,
             std::vector<Eigen::MatrixXd>, std::vector<Eigen::VectorXd>,
             std::vector<Eigen::MatrixXd>, std::vector<Eigen::VectorXd>,
             std::vector<Eigen::VectorXd>, std::vector<double>>
  precomputeSystemMatrices(int phase);

  std::tuple<std::vector<int>, std::vector<int>, std::vector<double>,
             std::vector<double>, std::vector<double>, Eigen::RowVectorXd>
  prepareLpMatrices(int phase);

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

  void run(const std::string& output_folder_path, int n_threads = 2);

  void dumpInitParamsToJson(const std::string& filepath) const;
};

}  // namespace barycentric_affine_approximator

// NOLINTEND(readability-identifier-naming)

#endif  // HCPWA_BARYCENTRIC_AFFINE_APPROXIMATOR_HPP
