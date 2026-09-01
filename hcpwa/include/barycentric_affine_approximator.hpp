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
#include "util/block_reduction_lp.hpp"
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

// Weight of the tie-break term epsilon * ||z||_1 in the main-step objective.
// Small enough that it only chooses among points the LP already considers
// optimal, large enough to survive the solver's own feasibility tolerance of
// 1e-6.
constexpr double kTieBreakEpsilon = 1e-8;

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

// The 8-row SparsePsi and the per-full-region ResidualRhsTerm that used to be
// declared here are gone. Psi is now stored per block region with coord_count
// rows (hcpwa::util::block_lp::BlockPsi), and the per-row right-hand-side data
// lives in the assembled reduced LP as block_lp::RhsTerm.

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

// Number of coordinate blocks per phase: A, B, C.
constexpr int kBlockCount = hcpwa::util::block_lp::kBlockCount;

// Geometry of one coordinate block of one phase.
//
// The five projection planes of a phase split into three groups sharing no
// coordinate, so an 8D region is the Cartesian product of three low-dimensional
// block polytopes (step 7, Lemma 1). Carrying the three block vertex sets
// instead of materialising the product is what makes the correct geometry
// affordable: the product has 108-192 vertices per region across 1.36M
// regions, the blocks have about 509 regions in total.
struct BlockGeometry {
  std::array<int, 3> coords{};     // state coordinates, ascending
  int coord_count = 0;             // 3 for A and B, 2 for C
  std::array<int, 2> layer_ids{};  // indices into PhaseGeometry::layers
  int layer_count = 0;             // 2 for A and B, 1 for C

  // local_axis[l][k] is the position within `coords` of the k-th axis of this
  // block's l-th layer, where l runs over 0..layer_count-1.
  //
  // Two of the six entries are NOT contiguous: block B of phase 0 has layer 3
  // on axes (1,6) inside coords (1,3,6), giving (0,2), and block B' of phase 1
  // has layer 2 on axes (3,7) inside coords (3,5,7), also giving (0,2).
  // Assuming (0,1) then (1,2) uniformly is wrong in exactly those two places
  // and produces a silently wrong phi rather than a crash, so this table is
  // derived by searching projectionAxesForPhase() and then asserted. It is
  // never written down as a constant in production code.
  std::array<std::array<int, 2>, 2> local_axis{};

  int num_regions = 0;
  std::vector<std::array<int, 2>> triangle_ids;
  // Per block region, all vertices, each of length coord_count.
  std::vector<std::vector<Eigen::VectorXd>> vertices;

  // Per block region, filled by precomputeSystemMatrices(). All restricted to
  // the block's own coordinates, so the matrices are coord_count square rather
  // than 8x8.
  std::vector<hcpwa::util::block_lp::BlockPsi> psi;
  std::vector<Eigen::MatrixXd> a_matr;
  std::vector<Eigen::VectorXd> f_vec;
  std::vector<Eigen::MatrixXd> q_c_matr;
  std::vector<Eigen::VectorXd> q_c_vec;
  std::vector<Eigen::MatrixXd> q_r_matr;
  std::vector<Eigen::VectorXd> q_r_vec;
  std::vector<Eigen::VectorXd> g_vec;
  std::vector<double> g_scal;
};

// Geometry for one phase.
//
// region_vertices is no longer populated: the 8D vertex lists it used to hold
// were truncated to 12 of an area's 108-192 vertices, and this path now works
// from `blocks` instead. region_triangle_ids survives because the border path
// indexes areas by their five simplex ids.
struct PhaseGeometry {
  std::array<ProjectionLayer, kSubsystemCount> layers;
  std::vector<std::array<int, kSubsystemCount>> region_triangle_ids;
  std::array<BlockGeometry, kBlockCount> blocks;

  // Number of 8D regions, i.e. M_A * M_B * M_C. Region ids run in the order
  //   region = (j_A * M_B + j_B) * M_C + j_C
  // which is the order the geometry pipeline assembles them in.
  int num_regions = 0;
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

// Layout of the x block of the LP variable vector: all unique 2D barycentric
// vertex values, grouped by projection layer.
//
// The y and mu columns are NOT described here. They belong to the reduced LP
// and are laid out by hcpwa::util::block_lp::ReducedLp, which owns
//   [ x | y_A | y_B | y_C | muR | muL ].
// The old num_regions and idxY(region, dim) fields, which indexed one
// 8-vector per full 8D region, are deliberately deleted rather than
// repurposed: leaving them in place would let a missed call site keep
// compiling while indexing a layout that no longer exists.
struct BarycentricVarLayout {
  std::array<int, kSubsystemCount> eta_s{};
  std::array<int, kSubsystemCount> offset_s{};
  int num_x = 0;

  int idxX(int subsystem, int vertex_id) const;
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

  // The per-region CTM matrices used to live here as eight parallel vectors
  // indexed by full region. They are now per block region and live inside
  // PhaseGeometry::blocks, alongside the geometry they are derived from.

  // The assembled reduced LP per phase. Owns the constraint matrix, the
  // objective, the column layout and the per-row data needed to refresh the
  // right-hand sides when x_next changes.
  std::array<hcpwa::util::block_lp::ReducedLp, kPhases> reduced_lps_;

  std::shared_ptr<spdlog::logger> logger_;
  std::vector<std::unique_ptr<Highs>> highs_solvers_;
  std::vector<std::vector<double>> row_lowers_;
  std::vector<std::unique_ptr<std::mutex>> solver_mutexes_;

  interval_building::ThetaTIndexLists theta_t_index_lists_;

  SparseVec buildPhiRow(int phase, int region, const Eigen::VectorXd& point,
                        double tolerance = kEps) const;

  // phi restricted to one block, evaluated at a block vertex.
  //
  // Rows are built blockwise but columns stay global: the entries land in the
  // same x columns as buildPhiRow() would put them, via idxX(layer_id, ...).
  // Mixing block-local and global column order is the classic way to get this
  // wrong, so the block index never reaches a column.
  SparseVec buildPhiRowBlock(int phase, int block_id, int j_block,
                             const Eigen::VectorXd& nu_block,
                             double tolerance = kEps) const;

  // Checks the block geometry invariants that do not need an LP: vertex counts,
  // every block vertex lying inside its own selected simplices, and each
  // block's phi summing to its layer count with no negative entries. Called
  // from getIntersectionPoints(); costs O(sum of block vertex counts), about
  // 3.5e3 rows on the production geometry.
  void validateBlockGeometry(int phase) const;

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

  // Tier-3 check: is the constructed function actually a bound?
  //
  // Samples random triples (a,b,c) from R_A x R_B x R_C and verifies the
  // ORIGINAL condition s*F(a,b,c) <= 0 at both endpoints of the segment, using
  // F = F_A + F_B + F_C (step 7, Lemma 2). This is the property the whole
  // construction exists to provide, and checking it directly is the only
  // verification available: the production product LP is unsolvable, so there
  // is no second implementation to diff against.
  //
  // It is also the check that would have caught the vertex truncation at once.
  // Under truncation only 12 of an area's 108-192 vertices were constrained,
  // so almost any sampled triple lands outside the enforced set and shows
  // s*F > 0.
  //
  // Cost is O(1) per triple with no LP and no 8D vertex, so the sample can be
  // large. Guarded by validate_.
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

  // Expands a block point into an 8-vector, padding the coordinates outside
  // the block with NaN. See the definition for why NaN and not zero.
  Eigen::VectorXd blockPointToFullState(int phase, int block_id,
                                        const Eigen::VectorXd& nu_block) const;

  // Centroid of one block region, computed from its own COMPLETE vertex set.
  //
  // This is the representative point at which each min branch is resolved. The
  // old whole-region centroid averaged a truncated vertex set whose hull was
  // segment x segment x triangle, a 4-dimensional slice of an 8-dimensional
  // region; midpoints of two vertices sharing a facet land on the region
  // boundary, where the min is tied, so branch resolution was being decided in
  // the kEps tie window or throwing outright. Complete per-block vertex sets
  // fix that as a side effect.
  Eigen::VectorXd blockCentroidCoords(int phase, int block_id,
                                      int j_block) const;

  // A, f, g restricted to one block, resolved at that block region's centroid.
  std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd, double>
  getBlockAMatrFVecGVecAndGScalJ(int phase, int block_id, int j_block) const;

  // Disturbance-box centre and radius maps restricted to one block.
  std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd,
             Eigen::VectorXd>
  getBlockQQ(int phase, int block_id, int j_block) const;

  // Fills the per-block-region Psi and CTM data in phase_geometries_[phase].
  void precomputeSystemMatrices(int phase);

  // Assembles the reduced LP of step 7 section 10 for one phase: three block
  // loops emitting (R), (L) and (Y) rows, then the two (Sigma) coupling rows.
  hcpwa::util::block_lp::ReducedLp prepareLpMatrices(int phase);

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
