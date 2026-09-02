#ifndef HCPWA_BARYCENTRIC_GEOMETRY_TYPES_HPP
#define HCPWA_BARYCENTRIC_GEOMETRY_TYPES_HPP

#include <Eigen/Core>
#include <algo.hpp>
#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
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

// Tolerance for the end-to-end residual check, as an absolute floor plus a
// term relative to the size of what is being measured.
//
// F is a sum of terms -- the slope, (Psi z)^T m, rho^T |Psi z|, g -- that very
// nearly cancel. On the N=160 geometry those terms reach 1e2, and the LP
// satisfies its rows only to the solver's own 1e-6, so F carries an absolute
// error of about 1e-4 before anything is wrong. A bare 1e-4 threshold sits
// exactly on that floor and fires on rounding: an observed worst of 1.20e-4
// against a value function of order 1e2 is 1e-6 relative.
//
// The check still asserts s * F <= 0 at every vertex of every region -- it is
// the statement that the result is a bound, and nothing about that weakens.
// Only the threshold is expressed in units the quantity actually has.
constexpr double kResidualValidationTol = 1e-4;
constexpr double kResidualValidationRelTol = 1e-6;

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

// The coordinate blocks of a phase, mirrored from hcpwa::BlockRegions.
constexpr int kBlockCount = hcpwa::kBlockCount;
constexpr int kMaxBlockCoords = hcpwa::kMaxBlockCoords;
constexpr int kMaxBlockLayers = hcpwa::kMaxBlockLayers;

// One coordinate block of one phase. The 8D region is the direct product of the
// three blocks, and every projection plane lies entirely inside one of them, so
// each block carries a self-contained slice of the geometry.
struct BlockGeometry {
  // Global state coordinates of this block, ascending.
  std::array<int, kMaxBlockCoords> coords{};
  int coord_count = 0;
  // Which of the five projection layers belong to this block.
  std::array<int, kMaxBlockLayers> layer_ids{};
  int layer_count = 0;
  // Where each layer's two axes sit inside coords. Two of the six entries are
  // non-adjacent -- phase 0 layer 3 is (0,2) and phase 1 layer 2 is (0,2) --
  // which is why this is derived by search rather than written down.
  std::array<std::array<int, 2>, kMaxBlockLayers> local_axis{};

  // Per block-region, one triangle id per layer of this block.
  std::vector<std::array<int, kMaxBlockLayers>> triangle_ids;
  // Per block-region, every vertex, in this block's own coordinate order.
  std::vector<std::vector<Eigen::VectorXd>> vertices;
  // Per block-region, the bounding box of those vertices.
  std::vector<Eigen::VectorXd> aabb_lower;
  std::vector<Eigen::VectorXd> aabb_upper;

  int numRegions() const { return static_cast<int>(vertices.size()); }
};

// Geometry for one phase. A full 8D region stores both its vertices and the
// five triangle ids used to select the local barycentric charts.
struct PhaseGeometry {
  std::array<ProjectionLayer, kSubsystemCount> layers;
  // The 8D product of the block cells. Empty in production: the LP and the
  // border solver are both assembled block by block, and materialising the
  // product costs prod_b |V_b| vertices per region. Tests that compare the two
  // paths ask the geometry layer for it explicitly.
  std::vector<std::vector<Eigen::VectorXd>> region_vertices;
  std::vector<std::array<int, kSubsystemCount>> region_triangle_ids;
  std::array<BlockGeometry, kBlockCount> blocks;
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

// The four CTM flows of a phase, as zero-based (from, to) cell pairs. Row `to`
// of A gains the flow, row `from` loses it, and g is their sum; that single
// rule reproduces the whole assembly of A, f, g and g0.
//
// Every flow reads only its own two cells and both of them lie in the same
// coordinate block, which is what lets the CTM data be built one block at a
// time.
inline std::array<std::array<int, 2>, 4> phaseFlows(int phase) {
  if (phase == 0) {
    // f31, f36, f24, f27 in the paper's one-based cell numbering.
    return {{{2, 0}, {2, 5}, {1, 3}, {1, 6}}};
  }
  if (phase == 1) {
    // f51, f57, f84, f86.
    return {{{4, 0}, {4, 6}, {7, 3}, {7, 5}}};
  }
  throw std::invalid_argument("phaseFlows: invalid phase");
}

// Recovers the coordinate blocks of a phase from projectionAxesForPhase alone:
// two state coordinates belong to the same block exactly when a chain of
// projection planes links them, so the blocks are the connected components of
// the plane incidence graph. Kept next to the axis table because it is the
// definition the geometry code has to agree with, and the assertions in the
// approximator compare the two.
//
// Components come out ordered by their smallest coordinate, which is not the
// order the geometry code emits blocks in (phase 1 emits {0,4,6}, {3,5,7},
// {1,2}). Compare the two as unordered collections.
inline std::array<std::vector<int>, kBlockCount> coordinateBlocksForPhase(
    int phase) {
  const auto axes = projectionAxesForPhase(phase);
  std::array<int, kSpaceDim> parent{};
  for (int i = 0; i < kSpaceDim; ++i) {
    parent[static_cast<std::size_t>(i)] = i;
  }
  const std::function<int(int)> find = [&parent](int i) {
    while (parent[static_cast<std::size_t>(i)] != i) {
      i = parent[static_cast<std::size_t>(i)];
    }
    return i;
  };
  for (const auto& pair : axes) {
    const int ra = find(pair[0]);
    const int rb = find(pair[1]);
    if (ra != rb) {
      parent[static_cast<std::size_t>(ra)] = rb;
    }
  }

  std::vector<int> roots;
  std::array<std::vector<int>, kBlockCount> blocks;
  for (int i = 0; i < kSpaceDim; ++i) {
    const int root = find(i);
    auto it = std::find(roots.begin(), roots.end(), root);
    if (it == roots.end()) {
      if (static_cast<int>(roots.size()) == kBlockCount) {
        throw std::runtime_error(
            "coordinateBlocksForPhase: more components than blocks");
      }
      roots.push_back(root);
      it = roots.end() - 1;
    }
    blocks[static_cast<std::size_t>(it - roots.begin())].push_back(i);
  }
  if (static_cast<int>(roots.size()) != kBlockCount) {
    throw std::runtime_error(
        "coordinateBlocksForPhase: fewer components than blocks");
  }
  return blocks;
}

// Layout of the x block of the LP variable vector: all unique 2D barycentric
// vertex values, one contiguous run per projection layer.
//
// The auxiliary y block is no longer described here. It used to be one
// 8-vector per product region; it is now one vector per block-region, and its
// columns -- together with the epigraph and tie-break columns -- are laid out
// by block_reduction::ReducedLpColLayout. Removing idxY and num_regions rather
// than repurposing them makes every missed call site a compile error.
struct BarycentricVarLayout {
  std::array<int, kSubsystemCount> eta_s{};
  std::array<int, kSubsystemCount> offset_s{};
  int num_x = 0;

  int idxX(int subsystem, int vertex_id) const;
};

// The value selector phi_{j,nu}: V(nu) = phi^T x. Free functions rather than
// members so that the block path and the product path can be compared on a
// hand-built PhaseGeometry, without running the arrangement.
//
// Both write global x columns. The block form reads only its own block's
// coordinates, in the block's order, and projects through local_axis; because
// the three blocks own disjoint sets of planes, summing the three block rows
// reproduces the product row exactly.
SparseVec buildPhiRow(const PhaseGeometry& geometry,
                      const BarycentricVarLayout& layout, int region,
                      const Eigen::VectorXd& point, double tolerance = kEps);

SparseVec buildPhiRowBlock(const PhaseGeometry& geometry,
                           const BarycentricVarLayout& layout, int block,
                           int block_region, const Eigen::VectorXd& point,
                           double tolerance = kEps);

// Converts the geometry layer's block factorisation into the form the LP uses,
// deriving local_axis and checking every structural invariant on the way:
// the blocks partition the eight coordinates and the five planes, and they are
// the connected components of the plane incidence graph.
std::array<BlockGeometry, kBlockCount> blockGeometryFromRegions(
    int phase, const std::array<hcpwa::BlockRegions, kBlockCount>& src);

// Where each of a block's planes puts its two axes inside the block's
// coordinate list. Throws if a plane straddles two blocks, which is the
// property the whole reduction rests on.
std::array<std::array<int, 2>, kMaxBlockLayers> localAxesForBlock(
    int phase, const std::array<int, kMaxBlockCoords>& coords, int coord_count,
    const std::array<int, kMaxBlockLayers>& layer_ids, int layer_count);

}  // namespace barycentric_affine_approximator

// NOLINTEND(readability-identifier-naming)

#endif  // HCPWA_BARYCENTRIC_GEOMETRY_TYPES_HPP
