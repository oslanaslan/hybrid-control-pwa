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

using ValueFunction = hcpwa::util::ValueFunction;

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

// Stores the row id and sparse phi_{j,nu} vector needed to update the dynamic
// RHS term:
//   row_upper = -g(nu) - phi_{j,nu}^T x_next / dt.
struct ResidualRhsTerm {
  int row_id = 0;
  SparseVec phi;
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

  ValueFunction value_function_;
  std::vector<Eigen::VectorXd> cube_angle_vertices_;
  std::vector<Eigen::VectorXd> common_refinement_vertices_;

  std::array<PhaseGeometry, kPhases> phase_geometries_;
  std::array<BarycentricVarLayout, kPhases> layouts_;

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

  int locateRegion(int phase, const Eigen::VectorXd& point,
                   double tolerance = kEps) const;

  double evaluateBarycentricValue(int phase, const std::vector<double>& x,
                                  const Eigen::VectorXd& point,
                                  double tolerance = kEps) const;

  std::vector<int> admissibleThetaIds(double theta) const;

 public:
  BarycentricAffineApproximator(double t_max, int t_split_count,
                                double tau_min, double tau_max,
                                const SystemParams& system_params,
                                bool highs_verbose = false);

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
