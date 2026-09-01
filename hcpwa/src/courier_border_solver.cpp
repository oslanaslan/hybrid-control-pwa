#include "courier_border_solver.hpp"

#include <Highs.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <format>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>

#include <spdlog/sinks/stdout_color_sinks.h>

// NOLINTBEGIN(readability-identifier-naming)

namespace barycentric_affine_approximator {

namespace {

// Same option block as the rest of the barycentric code, kept numerically
// identical so the border LP and the in-phase LPs agree on what "optimal" is.
constexpr double kHighsSolutionTol = 1e-6;
constexpr double kHighsSmallMatrixValue = 1e-9;

// Multipliers must satisfy the stationarity identities exactly up to rounding;
// this is a correctness gate, not a convergence knob.
constexpr double kDualIdentityTol = 1e-7;

void applyHighsOptions(Highs& highs, bool verbose, bool presolve) {
  highs.setOptionValue("solver", "simplex");
  highs.setOptionValue("presolve", presolve ? "on" : "off");
  highs.setOptionValue("simplex_strategy", 2);
  highs.setOptionValue("kkt_tolerance", kHighsSolutionTol);
  highs.setOptionValue("primal_feasibility_tolerance", kHighsSolutionTol);
  highs.setOptionValue("dual_feasibility_tolerance", kHighsSolutionTol);
  highs.setOptionValue("primal_residual_tolerance", kHighsSolutionTol);
  highs.setOptionValue("dual_residual_tolerance", kHighsSolutionTol);
  highs.setOptionValue("optimality_tolerance", kHighsSolutionTol);
  highs.setOptionValue("small_matrix_value", kHighsSmallMatrixValue);
  highs.setOptionValue("log_to_console", verbose);
  highs.changeObjectiveSense(ObjSense::kMinimize);
}

using detail::Point2;

// Sutherland-Hodgman clip of a convex polygon against one half-plane of an
// axis-aligned rectangle. `inside` is evaluated with a kGeomEps band and points
// exactly on the boundary count as inside, so an edge-incident vertex is never
// dropped -- dropping one would drop an LP row and could break soundness.
template <typename InsideFn, typename IntersectFn>
std::vector<Point2> clipHalfPlane(const std::vector<Point2>& poly,
                                  InsideFn inside, IntersectFn intersect) {
  std::vector<Point2> out;
  if (poly.empty()) {
    return out;
  }
  out.reserve(poly.size() + 1);
  for (std::size_t i = 0; i < poly.size(); ++i) {
    const Point2& cur = poly[i];
    const Point2& nxt = poly[(i + 1) % poly.size()];
    const bool cur_in = inside(cur);
    const bool nxt_in = inside(nxt);
    if (cur_in) {
      out.push_back(cur);
    }
    if (cur_in != nxt_in) {
      out.push_back(intersect(cur, nxt));
    }
  }
  return out;
}

}  // namespace

namespace detail {

std::vector<Point2> clipTriangleToRect(const std::array<Point2, 3>& tri,
                                       double lo_x, double hi_x, double lo_y,
                                       double hi_y) {
  std::vector<Point2> poly = {tri[0], tri[1], tri[2]};

  auto lerp = [](const Point2& a, const Point2& b, double t) {
    return Point2{a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)};
  };

  // x >= lo_x
  poly = clipHalfPlane(
      poly, [&](const Point2& p) { return p.x >= lo_x - kGeomEps; },
      [&](const Point2& a, const Point2& b) {
        const double d = b.x - a.x;
        return std::abs(d) < kGeomEps ? a : lerp(a, b, (lo_x - a.x) / d);
      });
  // x <= hi_x
  poly = clipHalfPlane(
      poly, [&](const Point2& p) { return p.x <= hi_x + kGeomEps; },
      [&](const Point2& a, const Point2& b) {
        const double d = b.x - a.x;
        return std::abs(d) < kGeomEps ? a : lerp(a, b, (hi_x - a.x) / d);
      });
  // y >= lo_y
  poly = clipHalfPlane(
      poly, [&](const Point2& p) { return p.y >= lo_y - kGeomEps; },
      [&](const Point2& a, const Point2& b) {
        const double d = b.y - a.y;
        return std::abs(d) < kGeomEps ? a : lerp(a, b, (lo_y - a.y) / d);
      });
  // y <= hi_y
  poly = clipHalfPlane(
      poly, [&](const Point2& p) { return p.y <= hi_y + kGeomEps; },
      [&](const Point2& a, const Point2& b) {
        const double d = b.y - a.y;
        return std::abs(d) < kGeomEps ? a : lerp(a, b, (hi_y - a.y) / d);
      });

  std::vector<Point2> out;
  for (const Point2& p : poly) {
    bool seen = false;
    for (const Point2& q : out) {
      if (std::abs(p.x - q.x) <= kGeomEps && std::abs(p.y - q.y) <= kGeomEps) {
        seen = true;
        break;
      }
    }
    if (!seen) {
      out.push_back(p);
    }
  }
  return out;
}

}  // namespace detail

namespace {

using detail::clipTriangleToRect;

// Uniform bucket grid over one projection layer. locateRegions() in the
// approximator is O(#regions) per query, which is unusable at the call counts
// here; this is O(triangles per bucket).
struct LayerIndex {
  int nx = 0;
  int ny = 0;
  double cw = 1.0;
  double ch = 1.0;
  double n_max = 0.0;
  std::vector<int> bucket_start;
  std::vector<int> bucket_tri;
  const ProjectionLayer* layer = nullptr;

  int bucketOf(double x, double y) const {
    int ix = static_cast<int>(x / cw);
    int iy = static_cast<int>(y / ch);
    ix = std::clamp(ix, 0, nx - 1);
    iy = std::clamp(iy, 0, ny - 1);
    return iy * nx + ix;
  }

  void build(const ProjectionLayer& l, double n) {
    layer = &l;
    n_max = n;
    const auto count = static_cast<int>(l.triangles.size());
    const int side = std::clamp(
        static_cast<int>(std::ceil(std::sqrt(std::max(1, count)))), 4, 256);
    nx = side;
    ny = side;
    cw = n_max / static_cast<double>(nx);
    ch = n_max / static_cast<double>(ny);

    std::vector<std::vector<int>> buckets(
        static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny));
    for (int t = 0; t < count; ++t) {
      const auto& tri = l.triangles[t];
      const double xs[3] = {static_cast<double>(tri.a[0]),
                            static_cast<double>(tri.b[0]),
                            static_cast<double>(tri.c[0])};
      const double ys[3] = {static_cast<double>(tri.a[1]),
                            static_cast<double>(tri.b[1]),
                            static_cast<double>(tri.c[1])};
      const double lo_x = *std::min_element(xs, xs + 3) - kGeomEps;
      const double hi_x = *std::max_element(xs, xs + 3) + kGeomEps;
      const double lo_y = *std::min_element(ys, ys + 3) - kGeomEps;
      const double hi_y = *std::max_element(ys, ys + 3) + kGeomEps;
      const int ix0 = std::clamp(static_cast<int>(lo_x / cw), 0, nx - 1);
      const int ix1 = std::clamp(static_cast<int>(hi_x / cw), 0, nx - 1);
      const int iy0 = std::clamp(static_cast<int>(lo_y / ch), 0, ny - 1);
      const int iy1 = std::clamp(static_cast<int>(hi_y / ch), 0, ny - 1);
      for (int iy = iy0; iy <= iy1; ++iy) {
        for (int ix = ix0; ix <= ix1; ++ix) {
          buckets[static_cast<std::size_t>(iy) * nx + ix].push_back(t);
        }
      }
    }

    bucket_start.assign(buckets.size() + 1, 0);
    for (std::size_t b = 0; b < buckets.size(); ++b) {
      bucket_start[b + 1]
          = bucket_start[b] + static_cast<int>(buckets[b].size());
    }
    bucket_tri.reserve(static_cast<std::size_t>(bucket_start.back()));
    for (const auto& b : buckets) {
      bucket_tri.insert(bucket_tri.end(), b.begin(), b.end());
    }
  }

  // Triangle ids whose bounding box meets the rectangle. Deduplicated and
  // sorted, so the clip pool order is deterministic.
  std::vector<int> overlapping(double lo_x, double hi_x, double lo_y,
                               double hi_y) const {
    const int ix0 = std::clamp(static_cast<int>(lo_x / cw), 0, nx - 1);
    const int ix1 = std::clamp(static_cast<int>(hi_x / cw), 0, nx - 1);
    const int iy0 = std::clamp(static_cast<int>(lo_y / ch), 0, ny - 1);
    const int iy1 = std::clamp(static_cast<int>(hi_y / ch), 0, ny - 1);
    std::vector<int> out;
    for (int iy = iy0; iy <= iy1; ++iy) {
      for (int ix = ix0; ix <= ix1; ++ix) {
        const int b = iy * nx + ix;
        out.insert(out.end(), bucket_tri.begin() + bucket_start[b],
                   bucket_tri.begin() + bucket_start[b + 1]);
      }
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
  }

  // Triangle containing (x, y) together with its barycentric coordinates.
  // Triangulations tile the square, so a miss is a hard error. A point on a
  // shared edge may resolve to either incident triangle; both give the same
  // value because neighbours share node ids.
  bool locate(double x, double y, int* tri_id, double alpha[3]) const {
    const int b = bucketOf(x, y);
    for (int k = bucket_start[b]; k < bucket_start[b + 1]; ++k) {
      const int t = bucket_tri[k];
      const TriangleBasis& basis = layer->bases[t];
      const Eigen::Vector3d a
          = basis.H * Eigen::Vector2d(x, y) + basis.h;
      if (a(0) >= -1e-7 && a(1) >= -1e-7 && a(2) >= -1e-7) {
        *tri_id = t;
        alpha[0] = a(0);
        alpha[1] = a(1);
        alpha[2] = a(2);
        return true;
      }
    }
    return false;
  }
};

// One vertex of a clipped polygon, carrying everything needed to evaluate any
// candidate's plane-s piece there: the 2D point and the three source-layer node
// columns with their barycentric weights.
struct ClipVertex {
  double gx = 0.0;
  double gy = 0.0;
  std::array<int, 3> col{};
  std::array<double, 3> alpha{};

  double eval(const std::vector<double>& x_src) const {
    return alpha[0] * x_src[static_cast<std::size_t>(col[0])]
           + alpha[1] * x_src[static_cast<std::size_t>(col[1])]
           + alpha[2] * x_src[static_cast<std::size_t>(col[2])];
  }
};

// Where one source plane's two axes land among the target blocks. The
// projection lemma says they always land in two *different* blocks, which is
// what makes the projection of a target region onto that plane an exact
// rectangle -- the product of one interval from each block -- rather than a
// bounding box of something larger.
struct SourcePlaneSplit {
  std::array<int, 2> block{};
  std::array<int, 2> local{};
};

// Prepared data for one target phase.
struct TargetData {
  // Per source plane s: the shared clip-vertex pool and the [begin, end) range
  // of each distinct rectangle in it.
  std::array<std::vector<ClipVertex>, kSubsystemCount> pool;
  std::array<std::vector<int>, kSubsystemCount> rect_begin;
  std::array<std::vector<int>, kSubsystemCount> rect_end;
  // Per source plane s, the interned rectangle of each block-cell pair, indexed
  // j_{b0} * M_{b1} + j_{b1}. This replaces one entry per product region: on
  // the N=100 arrangement that is 88 thousand pairs instead of 1.6 million
  // regions, and it is exact rather than a shared upper bound.
  std::array<std::vector<int>, kSubsystemCount> pair_rect;
  std::array<SourcePlaneSplit, kSubsystemCount> split{};

  // Deduplicated vertices and centroid of every block cell of the target
  // phase. Deduplicating once per block cell instead of once per region is the
  // same saving again: duplicate rows make the subproblem duals arbitrary.
  std::array<std::vector<std::vector<Eigen::VectorXd>>, kBlockCount>
      block_vertices;
  std::array<std::vector<Eigen::VectorXd>, kBlockCount> block_centroid;
  std::array<int, kBlockCount> block_regions{};
};

// Column layout of the region subproblem:
//   [ c_0 (2) d_0 | ... | c_4 (2) d_4 | zeta | t_A t_B t_C ]
constexpr int kSubZetaCol = 3 * kSubsystemCount;
constexpr int kSubCols = kSubZetaCol + 1 + kBlockCount;

int subTCol(int block) { return kSubZetaCol + 1 + block; }

// Phase-I result of one region, with the multipliers the Benders cut needs.
struct SubproblemResult {
  double zeta = 0.0;
  double kappa = 0.0;
  // Per block, one entry per unique vertex of that block's cell.
  std::array<std::vector<double>, kBlockCount> lambda;
  // Per emitted group-2 row, in order, together with its (plane, pool index).
  std::vector<double> mu;
  std::vector<std::pair<int, int>> group2;
};

}  // namespace

struct CourierBorderSolver::Impl {
  const std::array<PhaseGeometry, kPhases>* geometries = nullptr;
  const std::array<BarycentricVarLayout, kPhases>* layouts = nullptr;
  const std::array<Eigen::VectorXd, kPhases>* node_weights = nullptr;
  double n_max = 0.0;

  std::array<std::array<LayerIndex, kSubsystemCount>, kPhases> layer_index;
  std::array<TargetData, kPhases> target;

  // Number of product regions of a phase, and its decomposition into block
  // cell indices. The geometry nests its loops as
  //   for j_A { for j_B { for j_C } },
  // so j = (j_A * M_B + j_B) * M_C + j_C.
  int numRegions(int phase) const {
    const auto& blocks = (*geometries)[static_cast<std::size_t>(phase)].blocks;
    return blocks[0].numRegions() * blocks[1].numRegions()
           * blocks[2].numRegions();
  }

  std::array<int, kBlockCount> decodeRegion(int phase, int region) const {
    const auto& blocks = (*geometries)[static_cast<std::size_t>(phase)].blocks;
    std::array<int, kBlockCount> out{};
    for (int b = kBlockCount - 1; b >= 0; --b) {
      out[static_cast<std::size_t>(b)]
          = region % blocks[static_cast<std::size_t>(b)].numRegions();
      region /= blocks[static_cast<std::size_t>(b)].numRegions();
    }
    return out;
  }

  // Value of the target-phase barycentric function restricted to one block, at
  // one vertex of one block cell. Recomputed on the fly: caching phi rows would
  // still be large, and this is now a two-plane evaluation instead of five.
  double valueAtBlock(int phase, int block, int block_region,
                      const Eigen::VectorXd& nu,
                      const std::vector<double>& z) const {
    const PhaseGeometry& geom = (*geometries)[static_cast<std::size_t>(phase)];
    const BarycentricVarLayout& lay
        = (*layouts)[static_cast<std::size_t>(phase)];
    const BlockGeometry& bg = geom.blocks[static_cast<std::size_t>(block)];
    const auto& tri_ids
        = bg.triangle_ids[static_cast<std::size_t>(block_region)];
    double acc = 0.0;
    for (int l = 0; l < bg.layer_count; ++l) {
      const int s = bg.layer_ids[static_cast<std::size_t>(l)];
      const ProjectionLayer& layer = geom.layers[static_cast<std::size_t>(s)];
      const TriangleBasis& basis = layer.bases[static_cast<std::size_t>(
          tri_ids[static_cast<std::size_t>(l)])];
      const auto& local = bg.local_axis[static_cast<std::size_t>(l)];
      const Eigen::Vector3d a
          = basis.H * Eigen::Vector2d(nu(local[0]), nu(local[1])) + basis.h;
      for (int k = 0; k < 3; ++k) {
        acc += a(k)
               * z[static_cast<std::size_t>(lay.idxX(
                   s, basis.vertex_ids[static_cast<std::size_t>(k)]))];
      }
    }
    return acc;
  }

  // Phi_rho at an arbitrary point of Omega, evaluated on the source partition.
  bool evalSource(int phase, const Eigen::VectorXd& n,
                  const std::vector<double>& x_src, double* out) const {
    const BarycentricVarLayout& lay = (*layouts)[phase];
    double acc = 0.0;
    for (int s = 0; s < kSubsystemCount; ++s) {
      const LayerIndex& idx = layer_index[static_cast<std::size_t>(phase)][s];
      const ProjectionLayer& layer = (*geometries)[phase].layers[s];
      int tri = -1;
      double alpha[3];
      if (!idx.locate(n(layer.axes[0]), n(layer.axes[1]), &tri, alpha)) {
        return false;
      }
      const TriangleBasis& basis = layer.bases[static_cast<std::size_t>(tri)];
      for (int k = 0; k < 3; ++k) {
        acc += alpha[k]
               * x_src[static_cast<std::size_t>(
                   lay.idxX(s, basis.vertex_ids[k]))];
      }
    }
    *out = acc;
    return true;
  }
};

namespace {

// Builds and solves the phase-I subproblem of one target region, whose block
// cells are js. Condition (a) of the courier is
//   Vt(nu) - l_j(nu) <= 0  at every vertex of the region,
// and both sides split over the target blocks: Vt does because every target
// plane lies inside one block, and l_j does because each of its terms
// c_s[k] * nu_axis is a single coordinate, which lies in one block. So
//   max over the product of vertices  =  sum over blocks of the block maximum,
// and bounding each block maximum by t_b with sum_b t_b <= sum_s d_s + zeta is
// the same condition, checked at |V_A| + |V_B| + |V_C| vertices instead of
// their product. No vertex is dropped.
SubproblemResult solveRegionSubproblem(
    const CourierBorderSolver::Impl& impl, const TargetData& data, int tgt,
    const std::array<int, kBlockCount>& js, double sigma,
    const std::vector<double>& z,
    const std::function<double(int, int)>& target_value, bool want_duals,
    bool verbose) {
  const double kInf = std::numeric_limits<double>::infinity();

  std::vector<int> sstart = {0};
  std::vector<int> sindex;
  std::vector<double> svalue;
  std::vector<double> slower;
  std::vector<double> supper;

  // Group 1, block by block:
  //   sum_{(s,k) whose axis is in block b} c_s[k] nu_local + t_b
  //       >= sigma * V_b(nu; z).
  std::array<int, kBlockCount> g1_offset{};
  SubproblemResult out;
  for (int b = 0; b < kBlockCount; ++b) {
    g1_offset[static_cast<std::size_t>(b)]
        = static_cast<int>(slower.size());
    const auto& verts = data.block_vertices[static_cast<std::size_t>(b)][
        static_cast<std::size_t>(js[static_cast<std::size_t>(b)])];
    for (const Eigen::VectorXd& nu : verts) {
      for (int s = 0; s < kSubsystemCount; ++s) {
        const SourcePlaneSplit& split = data.split[static_cast<std::size_t>(s)];
        for (int k = 0; k < 2; ++k) {
          if (split.block[static_cast<std::size_t>(k)] != b) {
            continue;
          }
          sindex.push_back(3 * s + k);
          svalue.push_back(nu(split.local[static_cast<std::size_t>(k)]));
        }
      }
      sindex.push_back(subTCol(b));
      svalue.push_back(1.0);
      sstart.push_back(static_cast<int>(sindex.size()));
      slower.push_back(sigma
                       * impl.valueAtBlock(tgt, b,
                                           js[static_cast<std::size_t>(b)], nu,
                                           z));
      supper.push_back(kInf);
    }
  }

  // Coupling: sum_b t_b - sum_s d_s - zeta <= 0.
  const int coup_row = static_cast<int>(slower.size());
  for (int b = 0; b < kBlockCount; ++b) {
    sindex.push_back(subTCol(b));
    svalue.push_back(1.0);
  }
  for (int s = 0; s < kSubsystemCount; ++s) {
    sindex.push_back(3 * s + 2);
    svalue.push_back(-1.0);
  }
  sindex.push_back(kSubZetaCol);
  svalue.push_back(-1.0);
  sstart.push_back(static_cast<int>(sindex.size()));
  slower.push_back(-kInf);
  supper.push_back(0.0);

  // Group 2: c_s . g + d_s <= sigma * T_j^(s)(g). Unchanged by the reduction:
  // the rectangles it runs over are now indexed by block-cell pairs, but they
  // are the same rectangles.
  const int g2_offset = static_cast<int>(slower.size());
  for (int s = 0; s < kSubsystemCount; ++s) {
    const SourcePlaneSplit& split = data.split[static_cast<std::size_t>(s)];
    const int m1 = data.block_regions[static_cast<std::size_t>(split.block[1])];
    const std::size_t pair
        = static_cast<std::size_t>(js[static_cast<std::size_t>(split.block[0])])
              * m1
          + static_cast<std::size_t>(js[static_cast<std::size_t>(
              split.block[1])]);
    const int rect = data.pair_rect[s][pair];
    for (int e = data.rect_begin[s][static_cast<std::size_t>(rect)];
         e < data.rect_end[s][static_cast<std::size_t>(rect)]; ++e) {
      const ClipVertex& cv = data.pool[s][static_cast<std::size_t>(e)];
      sindex.push_back(3 * s + 0);
      svalue.push_back(cv.gx);
      sindex.push_back(3 * s + 1);
      svalue.push_back(cv.gy);
      sindex.push_back(3 * s + 2);
      svalue.push_back(1.0);
      sstart.push_back(static_cast<int>(sindex.size()));
      slower.push_back(-kInf);
      supper.push_back(sigma * target_value(s, e));
      out.group2.emplace_back(s, e);
    }
  }

  std::vector<double> sobj(static_cast<std::size_t>(kSubCols), 0.0);
  sobj[static_cast<std::size_t>(kSubZetaCol)] = 1.0;
  std::vector<double> sclo(static_cast<std::size_t>(kSubCols), -kInf);
  std::vector<double> schi(static_cast<std::size_t>(kSubCols), kInf);
  sclo[static_cast<std::size_t>(kSubZetaCol)] = 0.0;

  Highs sub;
  // presolve off: 19 columns make it pure overhead, and it removes any question
  // about postsolve dual recovery, which the cut depends on.
  applyHighsOptions(sub, verbose, /*presolve=*/false);
  if (sub.addCols(kSubCols, sobj.data(), sclo.data(), schi.data(), 0, nullptr,
                  nullptr, nullptr)
      != HighsStatus::kOk) {
    throw std::runtime_error("CourierBorderSolver: sub addCols failed");
  }
  if (sub.addRows(static_cast<int>(slower.size()), slower.data(),
                  supper.data(), static_cast<int>(svalue.size()),
                  sstart.data(), sindex.data(), svalue.data())
      != HighsStatus::kOk) {
    throw std::runtime_error("CourierBorderSolver: sub addRows failed");
  }
  if (sub.run() != HighsStatus::kOk
      || sub.getModelStatus() != HighsModelStatus::kOptimal) {
    throw std::runtime_error(std::format(
        "CourierBorderSolver: subproblem of block cells ({}, {}, {}) not "
        "optimal, status {}",
        js[0], js[1], js[2], static_cast<int>(sub.getModelStatus())));
  }

  const auto& ssol = sub.getSolution();
  out.zeta = ssol.col_value[static_cast<std::size_t>(kSubZetaCol)];
  if (!want_duals) {
    return out;
  }
  if (!ssol.dual_valid || ssol.row_dual.size() != slower.size()) {
    throw std::runtime_error(
        "CourierBorderSolver: no duals for the region subproblem");
  }

  // lambda on the >= rows, kappa and mu on the <= rows. Rather than trust the
  // sign convention, the caller verifies the stationarity identities that free
  // c_s, d_s, t_b and a basic zeta force. A flipped sign fails those
  // immediately instead of silently producing an invalid cut.
  out.kappa = -ssol.row_dual[static_cast<std::size_t>(coup_row)];
  for (int b = 0; b < kBlockCount; ++b) {
    const std::size_t count
        = data.block_vertices[static_cast<std::size_t>(b)][
              static_cast<std::size_t>(js[static_cast<std::size_t>(b)])].size();
    out.lambda[static_cast<std::size_t>(b)].resize(count);
    for (std::size_t i = 0; i < count; ++i) {
      out.lambda[static_cast<std::size_t>(b)][i] = ssol.row_dual[
          static_cast<std::size_t>(g1_offset[static_cast<std::size_t>(b)]) + i];
    }
  }
  out.mu.resize(out.group2.size());
  for (std::size_t i = 0; i < out.group2.size(); ++i) {
    out.mu[i] = -ssol.row_dual[static_cast<std::size_t>(g2_offset) + i];
  }
  return out;
}

}  // namespace

CourierBorderSolver::CourierBorderSolver(CourierBorderOptions options)
    : options_(options) {
  logger_ = spdlog::get("courier_border_solver");
  if (!logger_) {
    logger_ = spdlog::stdout_color_mt("courier_border_solver");
  }
  logger_->set_level(spdlog::level::info);
}

void CourierBorderSolver::prepare(
    const std::array<PhaseGeometry, kPhases>& geometries,
    const std::array<BarycentricVarLayout, kPhases>& layouts,
    const std::array<Eigen::VectorXd, kPhases>& node_weights, double n_max) {
  if (n_max <= 0.0) {
    throw std::invalid_argument("CourierBorderSolver::prepare: n_max must be positive");
  }
  auto impl = std::make_shared<Impl>();
  impl->geometries = &geometries;
  impl->layouts = &layouts;
  impl->node_weights = &node_weights;
  impl->n_max = n_max;

  for (int phase = 0; phase < kPhases; ++phase) {
    for (int s = 0; s < kSubsystemCount; ++s) {
      impl->layer_index[static_cast<std::size_t>(phase)][s].build(
          geometries[static_cast<std::size_t>(phase)].layers[s], n_max);
    }
  }

  // Deduplicated block-cell vertices and their centroids. A product region is
  // the direct product of three block cells, so its centroid is the
  // concatenation of theirs and its vertex set is their product; neither needs
  // to be materialised.
  for (int phase = 0; phase < kPhases; ++phase) {
    TargetData& data = impl->target[static_cast<std::size_t>(phase)];
    const auto& blocks = geometries[static_cast<std::size_t>(phase)].blocks;
    for (int b = 0; b < kBlockCount; ++b) {
      const BlockGeometry& bg = blocks[static_cast<std::size_t>(b)];
      if (bg.numRegions() == 0) {
        throw std::runtime_error(
            "CourierBorderSolver::prepare: target phase has no block cells");
      }
      data.block_regions[static_cast<std::size_t>(b)] = bg.numRegions();
      auto& verts = data.block_vertices[static_cast<std::size_t>(b)];
      auto& centroids = data.block_centroid[static_cast<std::size_t>(b)];
      verts.resize(static_cast<std::size_t>(bg.numRegions()));
      centroids.reserve(static_cast<std::size_t>(bg.numRegions()));
      for (int j = 0; j < bg.numRegions(); ++j) {
        const auto& raw = bg.vertices[static_cast<std::size_t>(j)];
        if (raw.empty()) {
          throw std::runtime_error(
              "CourierBorderSolver::prepare: block cell has no vertices");
        }
        // Duplicate rows in the subproblem make its multipliers arbitrary, so
        // they are removed once here rather than once per region.
        auto& unique = verts[static_cast<std::size_t>(j)];
        for (const auto& v : raw) {
          bool dup = false;
          for (const auto& u : unique) {
            if ((u - v).lpNorm<Eigen::Infinity>() <= kGeomEps) {
              dup = true;
              break;
            }
          }
          if (!dup) {
            unique.push_back(v);
          }
        }
        Eigen::VectorXd centroid = Eigen::VectorXd::Zero(bg.coord_count);
        for (const auto& v : raw) {
          centroid += v;
        }
        centroid /= static_cast<double>(raw.size());
        centroids.push_back(std::move(centroid));
      }
    }
  }

  std::size_t pool_bytes = 0;
  for (int phase = 0; phase < kPhases; ++phase) {
    const int src = 1 - phase;
    const PhaseGeometry& tgt_geom = geometries[static_cast<std::size_t>(phase)];
    const PhaseGeometry& src_geom = geometries[static_cast<std::size_t>(src)];
    const BarycentricVarLayout& src_lay = layouts[static_cast<std::size_t>(src)];
    TargetData& data = impl->target[static_cast<std::size_t>(phase)];

    for (int s = 0; s < kSubsystemCount; ++s) {
      const ProjectionLayer& src_layer = src_geom.layers[s];
      const LayerIndex& idx = impl->layer_index[static_cast<std::size_t>(src)][s];

      // Locate the plane's two axes among the target blocks. They must land in
      // different blocks: that is the projection lemma, and everything below
      // depends on it.
      SourcePlaneSplit& split = data.split[static_cast<std::size_t>(s)];
      for (int k = 0; k < 2; ++k) {
        const int axis = src_layer.axes[static_cast<std::size_t>(k)];
        split.block[static_cast<std::size_t>(k)] = -1;
        for (int b = 0; b < kBlockCount; ++b) {
          const BlockGeometry& bg
              = tgt_geom.blocks[static_cast<std::size_t>(b)];
          for (int c = 0; c < bg.coord_count; ++c) {
            if (bg.coords[static_cast<std::size_t>(c)] == axis) {
              split.block[static_cast<std::size_t>(k)] = b;
              split.local[static_cast<std::size_t>(k)] = c;
            }
          }
        }
        if (split.block[static_cast<std::size_t>(k)] < 0) {
          throw std::runtime_error(
              "CourierBorderSolver::prepare: source axis is in no target "
              "block");
        }
      }
      if (split.block[0] == split.block[1]) {
        throw std::runtime_error(std::format(
            "CourierBorderSolver::prepare: source plane {} of phase {} draws "
            "both axes from target block {}; the region projection is then not "
            "a rectangle and the method does not apply",
            s, src, split.block[0]));
      }

      const int b0 = split.block[0];
      const int b1 = split.block[1];
      const BlockGeometry& bg0 = tgt_geom.blocks[static_cast<std::size_t>(b0)];
      const BlockGeometry& bg1 = tgt_geom.blocks[static_cast<std::size_t>(b1)];
      const int m0 = bg0.numRegions();
      const int m1 = bg1.numRegions();

      data.pair_rect[s].assign(static_cast<std::size_t>(m0) * m1, -1);
      // Ordered map, not a hash: the pool order feeds LP row order, and HiGHS
      // picks different vertices of a degenerate optimum for different row
      // orders. Reproducibility matters more here than the lookup constant.
      std::map<std::array<std::int64_t, 4>, int> interned;

      for (int j0 = 0; j0 < m0; ++j0) {
        for (int j1 = 0; j1 < m1; ++j1) {
          // Exact, not a bound: the region is a product, so its projection on
          // this plane is the product of one interval from each block.
          double lo_x = bg0.aabb_lower[static_cast<std::size_t>(j0)](
              split.local[0]);
          double hi_x = bg0.aabb_upper[static_cast<std::size_t>(j0)](
              split.local[0]);
          double lo_y = bg1.aabb_lower[static_cast<std::size_t>(j1)](
              split.local[1]);
          double hi_y = bg1.aabb_upper[static_cast<std::size_t>(j1)](
              split.local[1]);
          // Inflating enlarges the rectangle, which only strengthens condition
          // (b); the certificate stays sound and merely gets more conservative.
          lo_x = std::max(0.0, lo_x - kGeomEps);
          lo_y = std::max(0.0, lo_y - kGeomEps);
          hi_x = std::min(n_max, hi_x + kGeomEps);
          hi_y = std::min(n_max, hi_y + kGeomEps);

          auto q = [](double v) {
            return static_cast<std::int64_t>(std::llround(v / kGeomEps));
          };
          const std::array<std::int64_t, 4> key = {q(lo_x), q(hi_x), q(lo_y),
                                                   q(hi_y)};
          const std::size_t pair = static_cast<std::size_t>(j0) * m1 + j1;
          auto it = interned.find(key);
          if (it != interned.end()) {
            data.pair_rect[s][pair] = it->second;
            continue;
          }

          const int rect_id = static_cast<int>(data.rect_begin[s].size());
          const std::size_t pool_before = data.pool[s].size();
          data.rect_begin[s].push_back(static_cast<int>(data.pool[s].size()));

          for (int t : idx.overlapping(lo_x, hi_x, lo_y, hi_y)) {
            const auto& tri = src_layer.triangles[static_cast<std::size_t>(t)];
            const std::array<Point2, 3> corners
                = {Point2{static_cast<double>(tri.a[0]),
                          static_cast<double>(tri.a[1])},
                   Point2{static_cast<double>(tri.b[0]),
                          static_cast<double>(tri.b[1])},
                   Point2{static_cast<double>(tri.c[0]),
                          static_cast<double>(tri.c[1])}};
            const std::vector<Point2> cut
                = clipTriangleToRect(corners, lo_x, hi_x, lo_y, hi_y);
            if (cut.empty()) {
              continue;
            }
            const TriangleBasis& basis
                = src_layer.bases[static_cast<std::size_t>(t)];
            for (const Point2& g : cut) {
              ClipVertex cv;
              cv.gx = g.x;
              cv.gy = g.y;
              const Eigen::Vector3d a
                  = basis.H * Eigen::Vector2d(g.x, g.y) + basis.h;
              for (int k = 0; k < 3; ++k) {
                cv.col[static_cast<std::size_t>(k)]
                    = src_lay.idxX(s, basis.vertex_ids[static_cast<std::size_t>(k)]);
                cv.alpha[static_cast<std::size_t>(k)] = a(k);
              }
              // Same point reached from two triangles gives the same value, the
              // barycentric function being continuous, so one entry is enough.
              bool dup = false;
              for (int e = data.rect_begin[s].back();
                   e < static_cast<int>(data.pool[s].size()); ++e) {
                const ClipVertex& q2 = data.pool[s][static_cast<std::size_t>(e)];
                if (std::abs(q2.gx - g.x) <= kGeomEps
                    && std::abs(q2.gy - g.y) <= kGeomEps) {
                  dup = true;
                  break;
                }
              }
              if (!dup) {
                data.pool[s].push_back(cv);
              }
            }
          }

          data.rect_end[s].push_back(static_cast<int>(data.pool[s].size()));
          if (data.rect_end[s].back() == data.rect_begin[s].back()) {
            throw std::runtime_error(std::format(
                "CourierBorderSolver::prepare: projection rectangle of block "
                "cells ({}, {}) on source plane {} meets no source triangle",
                j0, j1, s));
          }
          interned.emplace(key, rect_id);
          data.pair_rect[s][pair] = rect_id;

          pool_bytes += static_cast<std::size_t>(data.pool[s].size()
                                                 - pool_before) * sizeof(ClipVertex);
          if (pool_bytes > options_.max_prepare_bytes) {
            throw std::runtime_error(
                "CourierBorderSolver::prepare: clipped-vertex pool exceeds "
                "max_prepare_bytes at "
                + std::to_string(pool_bytes) + " bytes");
          }
        }
      }

      logger_->info(
          "courier prepare: target_phase={} source_plane={} block cells={}x{} "
          "pairs={} rectangles={} clip_vertices={}",
          phase, s, m0, m1, data.pair_rect[s].size(),
          data.rect_begin[s].size(), data.pool[s].size());
    }
  }

  impl_ = std::move(impl);
  prepared_ = true;
}

std::vector<double> CourierBorderSolver::solve(
    const CourierBorderRequest& request, ApproximationMode mode,
    CourierBorderStats* stats) const {
  if (!prepared_) {
    throw std::runtime_error("CourierBorderSolver::solve: prepare() not called");
  }
  if (request.target_phase < 0 || request.target_phase >= kPhases
      || request.source_phase < 0 || request.source_phase >= kPhases
      || request.target_phase == request.source_phase) {
    throw std::invalid_argument("CourierBorderSolver::solve: invalid phases");
  }
  if (request.candidates.empty()) {
    throw std::invalid_argument("CourierBorderSolver::solve: no candidates");
  }

  const Impl& impl = *impl_;
  const int tgt = request.target_phase;
  const int src = request.source_phase;
  const PhaseGeometry& tgt_geom = (*impl.geometries)[static_cast<std::size_t>(tgt)];
  const PhaseGeometry& src_geom = (*impl.geometries)[static_cast<std::size_t>(src)];
  const BarycentricVarLayout& tgt_lay = (*impl.layouts)[static_cast<std::size_t>(tgt)];
  const BarycentricVarLayout& src_lay = (*impl.layouts)[static_cast<std::size_t>(src)];
  const TargetData& data = impl.target[static_cast<std::size_t>(tgt)];
  const int n_regions = impl.numRegions(tgt);
  const auto n_cand = static_cast<int>(request.candidates.size());
  const int n_cols = tgt_lay.num_x;

  for (const auto& c : request.candidates) {
    if (static_cast<int>(c.size()) != src_lay.num_x) {
      throw std::invalid_argument(
          "CourierBorderSolver::solve: candidate has wrong size");
    }
  }

  // sigma = +1 for Lower, -1 for Upper. It enters exactly four places: the
  // master objective, the master row bounds, and the two subproblem RHS groups.
  const bool is_upper = mode == ApproximationMode::Upper;
  const double sigma = is_upper ? -1.0 : 1.0;

  // Upper needs the pointwise max over rho per clip vertex, which is
  // region-independent and computed once. Requiring l^(s) >= max_rho phi_rho^(s)
  // on every plane already implies sum_s l^(s) >= Phi_rho for every rho, so the
  // conjunction over rho collapses and Upper costs the same as Lower.
  std::array<std::vector<double>, kSubsystemCount> upper_T;
  if (is_upper) {
    for (int s = 0; s < kSubsystemCount; ++s) {
      const auto n_pool = data.pool[s].size();
      upper_T[s].assign(n_pool, -std::numeric_limits<double>::infinity());
      for (std::size_t e = 0; e < n_pool; ++e) {
        for (int c = 0; c < n_cand; ++c) {
          upper_T[s][e]
              = std::max(upper_T[s][e], data.pool[s][e].eval(request.candidates[c]));
        }
      }
    }
  }

  // Lower picks one candidate per region. Any choice is sound, so this is
  // purely about tightness; the centroid rule costs O(|P|) per region, whereas
  // "argmax over rho of the min over vertices" would cost O(|V_j| * |P| * 15).
  // The centroid of a region is the concatenation of its block centroids: the
  // region is a product, so the mean of each coordinate is the mean over that
  // block's own vertices. Nothing 8-dimensional is stored for it.
  auto regionCentroid = [&](const std::array<int, kBlockCount>& js) {
    Eigen::VectorXd centroid = Eigen::VectorXd::Zero(kSpaceDim);
    for (int b = 0; b < kBlockCount; ++b) {
      const BlockGeometry& bg = tgt_geom.blocks[static_cast<std::size_t>(b)];
      const Eigen::VectorXd& block_centroid
          = data.block_centroid[static_cast<std::size_t>(b)][
              static_cast<std::size_t>(js[static_cast<std::size_t>(b)])];
      for (int c = 0; c < bg.coord_count; ++c) {
        centroid(bg.coords[static_cast<std::size_t>(c)]) = block_centroid(c);
      }
    }
    return centroid;
  };

  std::vector<int> rho_of_region;
  if (!is_upper) {
    rho_of_region.assign(static_cast<std::size_t>(n_regions), 0);
    for (int j = 0; j < n_regions; ++j) {
      const Eigen::VectorXd centroid
          = regionCentroid(impl.decodeRegion(tgt, j));
      int best = 0;
      double best_val = -std::numeric_limits<double>::infinity();
      for (int c = 0; c < n_cand; ++c) {
        double val = 0.0;
        if (!impl.evalSource(src, centroid, request.candidates[c], &val)) {
          continue;
        }
        if (val > best_val) {
          best_val = val;
          best = c;
        }
      }
      rho_of_region[static_cast<std::size_t>(j)] = best;
    }
  }

  // Target value T_j^(s)(g) for one region and one source plane.
  auto targetValue = [&](int j, int s, int pool_idx) {
    if (is_upper) {
      return upper_T[s][static_cast<std::size_t>(pool_idx)];
    }
    return data.pool[s][static_cast<std::size_t>(pool_idx)].eval(
        request.candidates[static_cast<std::size_t>(
            rho_of_region[static_cast<std::size_t>(j)])]);
  };

  // One region vertex, assembled from one vertex of each block cell.
  auto regionVertex = [&](const std::array<int, kBlockCount>& js,
                          const std::array<std::size_t, kBlockCount>& ks) {
    Eigen::VectorXd nu = Eigen::VectorXd::Zero(kSpaceDim);
    for (int b = 0; b < kBlockCount; ++b) {
      const BlockGeometry& bg = tgt_geom.blocks[static_cast<std::size_t>(b)];
      const Eigen::VectorXd& v
          = data.block_vertices[static_cast<std::size_t>(b)][
              static_cast<std::size_t>(js[static_cast<std::size_t>(b)])][
              ks[static_cast<std::size_t>(b)]];
      for (int c = 0; c < bg.coord_count; ++c) {
        nu(bg.coords[static_cast<std::size_t>(c)]) = v(c);
      }
    }
    return nu;
  };

  // phi row of one region vertex, as the sum of the three block rows. Their
  // supports are disjoint -- each block owns its own projection planes -- so the
  // sum is the product row exactly.
  auto regionPhiRow = [&](const std::array<int, kBlockCount>& js,
                          const std::array<std::size_t, kBlockCount>& ks) {
    SparseVec phi;
    for (int b = 0; b < kBlockCount; ++b) {
      const SparseVec block_row = buildPhiRowBlock(
          tgt_geom, tgt_lay, b, js[static_cast<std::size_t>(b)],
          data.block_vertices[static_cast<std::size_t>(b)][
              static_cast<std::size_t>(js[static_cast<std::size_t>(b)])][
              ks[static_cast<std::size_t>(b)]],
          kEps);
      for (std::size_t k = 0; k < block_row.cols.size(); ++k) {
        phi.add(block_row.cols[k], block_row.vals[k], kGeomEps);
      }
    }
    return phi;
  };

  // ---- master ----------------------------------------------------------
  // Data-derived column box, so the very first solve is bounded before any cut
  // exists. It only restricts the feasible set, so it can cost tightness but
  // never soundness.
  double v_hi = 0.0;
  double v_lo = 0.0;
  for (int s = 0; s < kSubsystemCount; ++s) {
    double hi = -std::numeric_limits<double>::infinity();
    double lo = std::numeric_limits<double>::infinity();
    for (int k = 0; k < src_lay.eta_s[static_cast<std::size_t>(s)]; ++k) {
      for (int c = 0; c < n_cand; ++c) {
        const double v = request.candidates[static_cast<std::size_t>(c)]
                                           [static_cast<std::size_t>(
                                               src_lay.idxX(s, k))];
        hi = std::max(hi, v);
        lo = std::min(lo, v);
      }
    }
    v_hi += hi;
    v_lo += lo;
  }
  const double box = options_.master_box_scale * (v_hi - v_lo + 1.0);

  std::vector<double> col_lower(static_cast<std::size_t>(n_cols), -box);
  std::vector<double> col_upper(static_cast<std::size_t>(n_cols), box);
  // Same gauge fixing as initializeHighs(): layer 0 free, first node of every
  // later layer pinned to zero. Removes the additive barycentric nullspace
  // without changing the represented function. Any new LP over this basis must
  // use the identical convention or its output is not comparable.
  for (int s = 1; s < kSubsystemCount; ++s) {
    if (tgt_lay.eta_s[static_cast<std::size_t>(s)] == 0) {
      throw std::runtime_error(
          "CourierBorderSolver::solve: empty target projection layer");
    }
    const int col = tgt_lay.idxX(s, 0);
    col_lower[static_cast<std::size_t>(col)] = 0.0;
    col_upper[static_cast<std::size_t>(col)] = 0.0;
  }

  Eigen::VectorXd obj
      = -sigma * (*impl.node_weights)[static_cast<std::size_t>(tgt)];

  std::vector<int> row_start = {0};
  std::vector<int> row_index;
  std::vector<double> row_value;
  std::vector<double> row_lower;
  std::vector<double> row_upper;
  const double kInf = std::numeric_limits<double>::infinity();

  auto appendRow = [&](const SparseVec& row, double bound) {
    for (std::size_t k = 0; k < row.cols.size(); ++k) {
      row_index.push_back(row.cols[k]);
      row_value.push_back(row.vals[k]);
    }
    row_start.push_back(static_cast<int>(row_index.size()));
    if (is_upper) {
      row_lower.push_back(bound);
      row_upper.push_back(kInf);
    } else {
      row_lower.push_back(-kInf);
      row_upper.push_back(bound);
    }
  };

  // Seed rows: Vt(g) <= max_rho Phi_rho(g) at a stride sample of region
  // vertices. This is a necessary condition on the true solution in both modes,
  // so the rows are valid; they only speed up convergence.
  if (options_.max_seed_rows > 0) {
    const int stride = std::max(
        1, n_regions / std::max(1, options_.max_seed_rows / 8));
    int seeded = 0;
    for (int j = 0; j < n_regions && seeded < options_.max_seed_rows;
         j += stride) {
      const std::array<int, kBlockCount> js = impl.decodeRegion(tgt, j);
      std::array<std::size_t, kBlockCount> counts{};
      for (int b = 0; b < kBlockCount; ++b) {
        counts[static_cast<std::size_t>(b)]
            = data.block_vertices[static_cast<std::size_t>(b)][
                static_cast<std::size_t>(js[static_cast<std::size_t>(b)])]
                  .size();
      }
      for (std::size_t ka = 0; ka < counts[0] && seeded < options_.max_seed_rows;
           ++ka) {
        for (std::size_t kb = 0;
             kb < counts[1] && seeded < options_.max_seed_rows; ++kb) {
          for (std::size_t kc = 0;
               kc < counts[2] && seeded < options_.max_seed_rows; ++kc) {
            const std::array<std::size_t, kBlockCount> ks = {ka, kb, kc};
            const Eigen::VectorXd nu = regionVertex(js, ks);
            double best = -std::numeric_limits<double>::infinity();
            bool ok = false;
            for (int c = 0; c < n_cand; ++c) {
              double val = 0.0;
              if (impl.evalSource(
                      src, nu,
                      request.candidates[static_cast<std::size_t>(c)], &val)) {
                best = std::max(best, val);
                ok = true;
              }
            }
            if (!ok) {
              continue;
            }
            appendRow(regionPhiRow(js, ks), best);
            ++seeded;
          }
        }
      }
    }
  }

  CourierBorderStats local_stats;
  std::vector<double> z(static_cast<std::size_t>(n_cols), 0.0);

  for (int iter = 0; iter < options_.max_iterations; ++iter) {
    local_stats.iterations = iter + 1;

    Highs master;
    applyHighsOptions(master, options_.highs_verbose, /*presolve=*/true);
    if (master.addCols(n_cols, obj.data(), col_lower.data(), col_upper.data(), 0,
                       nullptr, nullptr, nullptr)
        != HighsStatus::kOk) {
      throw std::runtime_error("CourierBorderSolver::solve: master addCols failed");
    }
    const auto n_rows = static_cast<int>(row_lower.size());
    if (n_rows > 0
        && master.addRows(n_rows, row_lower.data(), row_upper.data(),
                          static_cast<int>(row_value.size()), row_start.data(),
                          row_index.data(), row_value.data())
               != HighsStatus::kOk) {
      throw std::runtime_error("CourierBorderSolver::solve: master addRows failed");
    }
    if (master.run() != HighsStatus::kOk
        || master.getModelStatus() != HighsModelStatus::kOptimal) {
      throw std::runtime_error(
          "CourierBorderSolver::solve: master LP not optimal, status "
          + std::to_string(static_cast<int>(master.getModelStatus()))
          + " at iteration " + std::to_string(iter));
    }
    const auto& sol = master.getSolution();
    for (int i = 0; i < n_cols; ++i) {
      z[static_cast<std::size_t>(i)] = sol.col_value[static_cast<std::size_t>(i)];
    }
    local_stats.master_objective = master.getObjectiveValue();
    local_stats.master_hit_box = false;
    for (int i = 0; i < n_cols; ++i) {
      if (std::abs(std::abs(z[static_cast<std::size_t>(i)]) - box) < 1e-6) {
        local_stats.master_hit_box = true;
        break;
      }
    }

    // ---- sweep every region ------------------------------------------
    struct Violation {
      double zeta;
      SparseVec row;
      double rhs;
    };
    std::vector<Violation> violations;
    double worst = 0.0;

    for (int j = 0; j < n_regions; ++j) {
      const std::array<int, kBlockCount> js = impl.decodeRegion(tgt, j);
      const SubproblemResult sub = solveRegionSubproblem(
          impl, data, tgt, js, sigma, z,
          [&](int s, int e) { return targetValue(j, s, e); },
          /*want_duals=*/true, /*verbose=*/false);
      ++local_stats.subproblems_solved;

      if (sub.zeta <= options_.certificate_tol) {
        continue;
      }
      worst = std::max(worst, sub.zeta);

      // Stationarity identities forced by free c_s, d_s, t_b and a basic zeta.
      // There are five families now instead of two, and they replace the old
      // pair rather than adding to it: kappa = 1 from zeta, one lambda mass per
      // block from t_b, one mu mass per plane from d_s, and ten barycentre
      // identities from the ten free c_s[k]. They are the only guard against a
      // sign mistake in the new "<=" coupling row against the ">=" group-1
      // rows.
      if (std::abs(sub.kappa - 1.0) > kDualIdentityTol) {
        throw std::runtime_error(
            "CourierBorderSolver::solve: coupling multiplier is "
            + std::to_string(sub.kappa) + ", expected 1");
      }
      for (int b = 0; b < kBlockCount; ++b) {
        double mass = 0.0;
        for (double v : sub.lambda[static_cast<std::size_t>(b)]) {
          if (v < -kDualIdentityTol) {
            throw std::runtime_error(
                "CourierBorderSolver::solve: negative lambda multiplier");
          }
          mass += v;
        }
        if (std::abs(mass - 1.0) > kDualIdentityTol) {
          throw std::runtime_error(
              "CourierBorderSolver::solve: lambda mass of block "
              + std::to_string(b) + " is " + std::to_string(mass)
              + ", expected 1");
        }
      }
      std::array<double, kSubsystemCount> mu_sum{};
      for (std::size_t i = 0; i < sub.group2.size(); ++i) {
        if (sub.mu[i] < -kDualIdentityTol) {
          throw std::runtime_error(
              "CourierBorderSolver::solve: negative mu multiplier");
        }
        mu_sum[static_cast<std::size_t>(sub.group2[i].first)] += sub.mu[i];
      }
      for (int s = 0; s < kSubsystemCount; ++s) {
        if (std::abs(mu_sum[static_cast<std::size_t>(s)] - 1.0)
            > kDualIdentityTol) {
          throw std::runtime_error(
              "CourierBorderSolver::solve: mu mass of plane "
              + std::to_string(s) + " is "
              + std::to_string(mu_sum[static_cast<std::size_t>(s)])
              + ", expected 1");
        }
      }
      // Ten barycentre identities, one per free c_s[k]: the lambda-weighted
      // mean of that axis over the owning block's vertices must equal the
      // mu-weighted mean of the same coordinate over the plane's clip vertices.
      for (int s = 0; s < kSubsystemCount; ++s) {
        const SourcePlaneSplit& split = data.split[static_cast<std::size_t>(s)];
        for (int k = 0; k < 2; ++k) {
          const int b = split.block[static_cast<std::size_t>(k)];
          const auto& verts = data.block_vertices[static_cast<std::size_t>(b)][
              static_cast<std::size_t>(js[static_cast<std::size_t>(b)])];
          double left = 0.0;
          for (std::size_t i = 0; i < verts.size(); ++i) {
            left += sub.lambda[static_cast<std::size_t>(b)][i]
                    * verts[i](split.local[static_cast<std::size_t>(k)]);
          }
          double right = 0.0;
          for (std::size_t i = 0; i < sub.group2.size(); ++i) {
            if (sub.group2[i].first != s) {
              continue;
            }
            const ClipVertex& cv = data.pool[s][static_cast<std::size_t>(
                sub.group2[i].second)];
            right += sub.mu[i] * (k == 0 ? cv.gx : cv.gy);
          }
          if (std::abs(left - right)
              > kDualIdentityTol * std::max(1.0, impl.n_max)) {
            throw std::runtime_error(std::format(
                "CourierBorderSolver::solve: barycentre identity of plane {} "
                "axis {} fails, {} vs {}", s, k, left, right));
          }
        }
      }

      // Cut: sigma * (sum_b sum_i lambda_{b,i} phi_b(nu_b^i))^T z
      //          <= sigma * sum mu T.
      // Not weaker than the product cut: it coincides with it for the product
      // multiplier lambda_A (x) lambda_B (x) lambda_C, and it cuts off exactly
      // the same zeta*.
      SparseVec row;
      for (int b = 0; b < kBlockCount; ++b) {
        const auto& verts = data.block_vertices[static_cast<std::size_t>(b)][
            static_cast<std::size_t>(js[static_cast<std::size_t>(b)])];
        for (std::size_t i = 0; i < verts.size(); ++i) {
          const double weight = sub.lambda[static_cast<std::size_t>(b)][i];
          if (weight == 0.0) {
            continue;
          }
          const SparseVec phi = buildPhiRowBlock(
              tgt_geom, tgt_lay, b, js[static_cast<std::size_t>(b)], verts[i],
              kEps);
          for (std::size_t k = 0; k < phi.cols.size(); ++k) {
            row.add(phi.cols[k], weight * phi.vals[k], kGeomEps);
          }
        }
      }
      // Each block contributes one unit of barycentric mass per plane it owns,
      // and lambda is a distribution on each block, so the aggregate row still
      // sums to the five planes.
      double row_mass = 0.0;
      for (double v : row.vals) {
        row_mass += v;
      }
      if (std::abs(row_mass - static_cast<double>(kSubsystemCount))
          > kDualIdentityTol * kSubsystemCount) {
        throw std::runtime_error(
            "CourierBorderSolver::solve: cut row mass is "
            + std::to_string(row_mass) + ", expected "
            + std::to_string(kSubsystemCount));
      }

      double rhs = 0.0;
      for (std::size_t i = 0; i < sub.group2.size(); ++i) {
        rhs += sub.mu[i]
               * targetValue(j, sub.group2[i].first, sub.group2[i].second);
      }

      // A cut that does not separate the current point means the multipliers
      // are wrong and Benders would stall silently.
      const double sep = sigma * (row.dot(z) - rhs);
      if (sep < 0.5 * sub.zeta) {
        throw std::runtime_error(
            "CourierBorderSolver::solve: generated cut does not separate z*, "
            "separation " + std::to_string(sep) + " vs zeta "
            + std::to_string(sub.zeta));
      }
      violations.push_back(Violation{sub.zeta, std::move(row), rhs});
    }

    local_stats.worst_zeta = worst;

    if (violations.empty()) {
      if (stats != nullptr) {
        *stats = local_stats;
      }
      return z;
    }

    std::sort(violations.begin(), violations.end(),
              [](const Violation& a, const Violation& b) {
                return a.zeta > b.zeta;
              });
    const int take = std::min(static_cast<int>(violations.size()),
                              options_.max_cuts_per_iteration);
    for (int i = 0; i < take; ++i) {
      appendRow(violations[static_cast<std::size_t>(i)].row,
                violations[static_cast<std::size_t>(i)].rhs);
      ++local_stats.cuts_added;
    }
  }

  if (stats != nullptr) {
    *stats = local_stats;
  }
  // Deliberately no partial result and no fallback: run() writes whatever comes
  // back into the value function and marches on it for the rest of the backward
  // sweep, so an uncertified vector would poison everything downstream.
  throw std::runtime_error(
      "CourierBorderSolver::solve: not certified after "
      + std::to_string(options_.max_iterations)
      + " iterations; worst zeta = " + std::to_string(local_stats.worst_zeta)
      + ", target_phase = " + std::to_string(request.target_phase)
      + ", theta_idx = " + std::to_string(request.theta_idx)
      + ", switch_cnt = " + std::to_string(request.switch_cnt));
}

double CourierBorderSolver::worstCertificateResidual(
    const CourierBorderRequest& request, ApproximationMode mode,
    const std::vector<double>& z) const {
  if (!prepared_) {
    throw std::runtime_error(
        "CourierBorderSolver::worstCertificateResidual: prepare() not called");
  }
  // Rebuild each region's courier from scratch against the given z and report
  // the worst phase-I objective. Zero means every region is certified.
  const Impl& impl = *impl_;
  const int tgt = request.target_phase;
  const PhaseGeometry& tgt_geom
      = (*impl.geometries)[static_cast<std::size_t>(tgt)];
  const bool is_upper = mode == ApproximationMode::Upper;
  const double sigma = is_upper ? -1.0 : 1.0;
  const TargetData& data = impl.target[static_cast<std::size_t>(tgt)];
  const auto n_cand = static_cast<int>(request.candidates.size());
  const int n_regions = impl.numRegions(tgt);

  std::vector<int> rho_of_region;
  if (!is_upper) {
    rho_of_region.assign(static_cast<std::size_t>(n_regions), 0);
    for (int j = 0; j < n_regions; ++j) {
      const std::array<int, kBlockCount> js = impl.decodeRegion(tgt, j);
      Eigen::VectorXd centroid = Eigen::VectorXd::Zero(kSpaceDim);
      for (int b = 0; b < kBlockCount; ++b) {
        const BlockGeometry& bg = tgt_geom.blocks[static_cast<std::size_t>(b)];
        const Eigen::VectorXd& block_centroid
            = data.block_centroid[static_cast<std::size_t>(b)][
                static_cast<std::size_t>(js[static_cast<std::size_t>(b)])];
        for (int c = 0; c < bg.coord_count; ++c) {
          centroid(bg.coords[static_cast<std::size_t>(c)]) = block_centroid(c);
        }
      }
      double best_val = -std::numeric_limits<double>::infinity();
      for (int c = 0; c < n_cand; ++c) {
        double val = 0.0;
        if (impl.evalSource(request.source_phase, centroid,
                            request.candidates[static_cast<std::size_t>(c)],
                            &val)
            && val > best_val) {
          best_val = val;
          rho_of_region[static_cast<std::size_t>(j)] = c;
        }
      }
    }
  }

  double worst = 0.0;
  for (int j = 0; j < n_regions; ++j) {
    const std::array<int, kBlockCount> js = impl.decodeRegion(tgt, j);
    // Same candidate solve() would pick. Soundness only needs a courier for
    // SOME rho -- if one exists for rho(j) then Vt <= Phi_{rho(j)} <=
    // max_rho Phi_rho -- so re-verifying against this rho is a valid check,
    // whereas the min over all candidates would demand something the method
    // never claimed.
    auto target_value = [&](int s, int e) {
      const ClipVertex& cv = data.pool[s][static_cast<std::size_t>(e)];
      if (!is_upper) {
        return cv.eval(request.candidates[static_cast<std::size_t>(
            rho_of_region[static_cast<std::size_t>(j)])]);
      }
      double best = -std::numeric_limits<double>::infinity();
      for (int c = 0; c < n_cand; ++c) {
        best = std::max(
            best, cv.eval(request.candidates[static_cast<std::size_t>(c)]));
      }
      return best;
    };

    const SubproblemResult sub = solveRegionSubproblem(
        impl, data, tgt, js, sigma, z, target_value, /*want_duals=*/false,
        /*verbose=*/false);
    worst = std::max(worst, sub.zeta);
  }
  return worst;
}

}  // namespace barycentric_affine_approximator

// NOLINTEND(readability-identifier-naming)
