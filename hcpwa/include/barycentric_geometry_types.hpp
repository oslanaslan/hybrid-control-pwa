#ifndef HCPWA_BARYCENTRIC_GEOMETRY_TYPES_HPP
#define HCPWA_BARYCENTRIC_GEOMETRY_TYPES_HPP

#include <Eigen/Core>
#include <algo.hpp>
#include <array>
#include <stdexcept>
#include <vector>

#include "util/value_function_utils.hpp"

// Geometry and layout types shared by BarycentricAffineApproximator and
// CourierBorderSolver. Extracted verbatim from barycentric_affine_approximator
// .hpp so the border solver can depend on the geometry without depending on the
// approximator, which holds a solver by value and would otherwise close an
// include cycle.
//
// Everything here keeps the naming style of the approximator header.
// NOLINTBEGIN(readability-identifier-naming)

namespace barycentric_affine_approximator {

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

// Sparse vector used for formula-level LP coefficients such as phi_{j,nu} and
// a_{j,nu}. The representation is deliberately simple: the barycentric LP is
// difficult to debug, so we keep sparse operations explicit and local.
struct SparseVec {
  std::vector<int> cols;
  std::vector<double> vals;

  void add(int col, double value, double eps = kEps);
  double dot(const std::vector<double>& x) const;
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

// The two axes each projection plane is built on, zero-based.
//
// This table is a contract shared by three places and must not drift: the prism
// tuple order produced by compute_intersection_points() in user_algo.cpp, the
// layer order of PhaseGeometry, and the projection lemma the courier border
// solver relies on -- every plane of one phase must draw its two axes from two
// different coordinate groups of the other phase.
inline std::array<std::array<int, 2>, kSubsystemCount> projectionAxesForPhase(
    int phase) {
  if (phase == 0) {
    // Phase 0 tuple order: [31, 36, 24, 27, 58].
    return {{{0, 2}, {2, 5}, {1, 3}, {1, 6}, {4, 7}}};
  }
  if (phase == 1) {
    // Phase 1 tuple order: [51, 57, 84, 86, 23].
    return {{{0, 4}, {4, 6}, {3, 7}, {5, 7}, {1, 2}}};
  }
  throw std::invalid_argument("projectionAxesForPhase: invalid phase");
}

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

}  // namespace barycentric_affine_approximator

// NOLINTEND(readability-identifier-naming)

#endif  // HCPWA_BARYCENTRIC_GEOMETRY_TYPES_HPP
