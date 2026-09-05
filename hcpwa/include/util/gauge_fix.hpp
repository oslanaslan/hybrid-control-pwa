#ifndef HCPWA_UTIL_GAUGE_FIX_HPP
#define HCPWA_UTIL_GAUGE_FIX_HPP

#include <vector>

#include "barycentric_geometry_types.hpp"
#include "util/block_reduction_lp.hpp"

// The gauge normalisation of the node values (paper, remark on non-uniqueness
// of the representation, condition gauge_fix).
//
// The map from node values z to the function V has a kernel: within a block,
// two projection planes share a coordinate c, and any piecewise-linear g(c)
// with breaks only on lines {c = const} that are edges of *both* planes'
// triangulations can be added to one plane and subtracted from the other;
// and, in the coupled form of the band LP, a constant can be moved from a
// plane of block A or B to the plane of block C. None of V, its gradient, the
// residuals or the integral see these directions, so an LP over z is
// degenerate along them unless they are pinned.
//
// Pins, per phase:
//   line_pins      on the first plane of each shared-coordinate pair, the
//                  nodes on the boundary line {other coordinate = 0} at the
//                  shared-coordinate values Xi where a line is an edge of
//                  both triangulations -- one pin per kernel function g;
//   constant_pins  one node of the second plane of each pair, against the
//                  constant exchange with block C (coupled form only);
//   group_c_columns every node of the block-C plane, which at level r = 0 is
//                  identically zero (remark on group C).
//
// pinsForLevel() says which apply at a level; the same set must be used in
// every LP that shares the node values of that level (band LP and courier
// master alike). buildGaugeFix() derives everything from the geometry, and
// verifyGaugeFix() checks it against the numerically computed kernel.
//
// NOLINTBEGIN(readability-identifier-naming)

namespace barycentric_affine_approximator {

struct GaugeFix {
  // One shared-coordinate pair of planes.
  struct Pair {
    int block = -1;
    int layer_plus = -1;   // carries the line pins
    int layer_minus = -1;  // carries the constant pin
    int shared_coord = -1;
    // Values of the shared coordinate whose line is an edge of both
    // triangulations, ascending. Always contains both ends of the square.
    std::vector<double> xi;
  };
  std::vector<Pair> pairs;
  int group_c_block = -1;
  int group_c_layer = -1;

  std::vector<int> line_pins;
  std::vector<int> constant_pins;
  std::vector<int> group_c_columns;

  // Sorted, unique. r = 0: line pins and the whole block-C plane. r >= 1:
  // line pins and the constant pins. The constant pins must NOT be added at
  // r = 0: with the C plane at zero the constants of the other planes are
  // meaningful and pinning one would narrow the class.
  std::vector<int> pinsForLevel(int switch_cnt) const;
};

GaugeFix buildGaugeFix(const PhaseGeometry& geometry,
                       const BarycentricVarLayout& layout, int phase);

struct GaugeFixReport {
  // Dimension of the kernel of [phi rows; Psi rows] over the x columns.
  int kernel_dim = 0;
  // Rank of that kernel restricted to the line pins; must equal kernel_dim.
  int line_rank = 0;
  // Rank of the kernel extended by the two constant exchanges, restricted to
  // line and constant pins; must equal kernel_dim + 2.
  int full_rank = 0;
  double sigma_max = 0.0;
  double sigma_min_kept = 0.0;
  double sigma_max_dropped = 0.0;
};

// Throws unless the pins meet the kernel transversally: the line pins alone
// span the block-form kernel, the line and constant pins together span it
// extended by the constant exchanges, the kernel does not touch the C plane,
// and its dimension is the sum of |Xi| over the pairs, which is what the
// remark predicts. Cost: one SVD of a (sum_b R_b + Psi rows) x num_x matrix.
GaugeFixReport verifyGaugeFix(const GaugeFix& gauge,
                              const block_reduction::ReducedLpInput& input,
                              const BarycentricVarLayout& layout);

}  // namespace barycentric_affine_approximator

// NOLINTEND(readability-identifier-naming)

#endif  // HCPWA_UTIL_GAUGE_FIX_HPP
