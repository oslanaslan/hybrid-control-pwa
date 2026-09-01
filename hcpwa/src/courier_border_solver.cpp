#include "courier_border_solver.hpp"

#include <Highs.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
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

// Prepared data for one target phase.
struct TargetData {
  // Per source plane s: the shared clip-vertex pool, the [begin, end) range of
  // each distinct rectangle in it, and the rectangle each region maps to.
  std::array<std::vector<ClipVertex>, kSubsystemCount> pool;
  std::array<std::vector<int>, kSubsystemCount> rect_begin;
  std::array<std::vector<int>, kSubsystemCount> rect_end;
  std::array<std::vector<int>, kSubsystemCount> region_rect;
};

}  // namespace

struct CourierBorderSolver::Impl {
  const std::array<PhaseGeometry, kPhases>* geometries = nullptr;
  const std::array<BarycentricVarLayout, kPhases>* layouts = nullptr;
  const std::array<Eigen::VectorXd, kPhases>* node_weights = nullptr;
  double n_max = 0.0;

  std::array<std::array<LayerIndex, kSubsystemCount>, kPhases> layer_index;
  std::array<TargetData, kPhases> target;

  // Value of the target-phase barycentric function at one region vertex.
  // Recomputed on the fly: caching phi rows for every (region, vertex) pair
  // would be tens of gigabytes.
  void phiRow(int phase, int region, const Eigen::VectorXd& nu,
              SparseVec* out) const {
    const PhaseGeometry& geom = (*geometries)[phase];
    const BarycentricVarLayout& lay = (*layouts)[phase];
    const auto& tri_ids = geom.region_triangle_ids[static_cast<std::size_t>(region)];
    for (int s = 0; s < kSubsystemCount; ++s) {
      const ProjectionLayer& layer = geom.layers[s];
      const TriangleBasis& basis = layer.bases[static_cast<std::size_t>(tri_ids[s])];
      const Eigen::Vector3d a
          = basis.H * Eigen::Vector2d(nu(layer.axes[0]), nu(layer.axes[1]))
            + basis.h;
      for (int k = 0; k < 3; ++k) {
        // kGeomEps, never kEps: the default tolerance of SparseVec::add is
        // large enough to delete real barycentric mass near a simplex edge.
        out->add(lay.idxX(s, basis.vertex_ids[k]), a(k), kGeomEps);
      }
    }
  }

  double valueAt(int phase, int region, const Eigen::VectorXd& nu,
                 const std::vector<double>& z) const {
    const PhaseGeometry& geom = (*geometries)[phase];
    const BarycentricVarLayout& lay = (*layouts)[phase];
    const auto& tri_ids = geom.region_triangle_ids[static_cast<std::size_t>(region)];
    double acc = 0.0;
    for (int s = 0; s < kSubsystemCount; ++s) {
      const ProjectionLayer& layer = geom.layers[s];
      const TriangleBasis& basis = layer.bases[static_cast<std::size_t>(tri_ids[s])];
      const Eigen::Vector3d a
          = basis.H * Eigen::Vector2d(nu(layer.axes[0]), nu(layer.axes[1]))
            + basis.h;
      for (int k = 0; k < 3; ++k) {
        acc += a(k)
               * z[static_cast<std::size_t>(lay.idxX(s, basis.vertex_ids[k]))];
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

  std::size_t pool_bytes = 0;
  for (int phase = 0; phase < kPhases; ++phase) {
    const int src = 1 - phase;
    const PhaseGeometry& tgt_geom = geometries[static_cast<std::size_t>(phase)];
    const PhaseGeometry& src_geom = geometries[static_cast<std::size_t>(src)];
    const BarycentricVarLayout& src_lay = layouts[static_cast<std::size_t>(src)];
    TargetData& data = impl->target[static_cast<std::size_t>(phase)];
    const auto n_regions = static_cast<int>(tgt_geom.region_vertices.size());
    if (n_regions == 0) {
      throw std::runtime_error(
          "CourierBorderSolver::prepare: target phase has no regions");
    }

    for (int s = 0; s < kSubsystemCount; ++s) {
      const ProjectionLayer& src_layer = src_geom.layers[s];
      const LayerIndex& idx = impl->layer_index[static_cast<std::size_t>(src)][s];
      const int ax0 = src_layer.axes[0];
      const int ax1 = src_layer.axes[1];

      data.region_rect[s].assign(static_cast<std::size_t>(n_regions), -1);
      // Ordered map, not a hash: the pool order feeds LP row order, and HiGHS
      // picks different vertices of a degenerate optimum for different row
      // orders. Reproducibility matters more here than the lookup constant.
      std::map<std::array<std::int64_t, 4>, int> interned;

      for (int j = 0; j < n_regions; ++j) {
        const auto& verts = tgt_geom.region_vertices[static_cast<std::size_t>(j)];
        if (verts.empty()) {
          throw std::runtime_error(
              "CourierBorderSolver::prepare: region has no vertices");
        }
        double lo_x = verts[0](ax0);
        double hi_x = lo_x;
        double lo_y = verts[0](ax1);
        double hi_y = lo_y;
        for (const auto& v : verts) {
          lo_x = std::min(lo_x, v(ax0));
          hi_x = std::max(hi_x, v(ax0));
          lo_y = std::min(lo_y, v(ax1));
          hi_y = std::max(hi_y, v(ax1));
        }
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
        auto it = interned.find(key);
        if (it != interned.end()) {
          data.region_rect[s][static_cast<std::size_t>(j)] = it->second;
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
          throw std::runtime_error(
              "CourierBorderSolver::prepare: projection rectangle of region "
              + std::to_string(j) + " on source plane " + std::to_string(s)
              + " meets no source triangle");
        }
        interned.emplace(key, rect_id);
        data.region_rect[s][static_cast<std::size_t>(j)] = rect_id;

        pool_bytes += static_cast<std::size_t>(data.pool[s].size()
                                               - pool_before) * sizeof(ClipVertex);
        if (pool_bytes > options_.max_prepare_bytes) {
          throw std::runtime_error(
              "CourierBorderSolver::prepare: clipped-vertex pool exceeds "
              "max_prepare_bytes at "
              + std::to_string(pool_bytes) + " bytes");
        }
      }

      logger_->info(
          "courier prepare: target_phase={} source_plane={} regions={} "
          "rectangles={} clip_vertices={}",
          phase, s, n_regions, data.rect_begin[s].size(), data.pool[s].size());
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
  const auto n_regions = static_cast<int>(tgt_geom.region_vertices.size());
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
  std::vector<int> rho_of_region;
  if (!is_upper) {
    rho_of_region.assign(static_cast<std::size_t>(n_regions), 0);
    for (int j = 0; j < n_regions; ++j) {
      const auto& verts = tgt_geom.region_vertices[static_cast<std::size_t>(j)];
      Eigen::VectorXd centroid = Eigen::VectorXd::Zero(kSpaceDim);
      for (const auto& v : verts) {
        centroid += v;
      }
      centroid /= static_cast<double>(verts.size());
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
      const auto& verts = tgt_geom.region_vertices[static_cast<std::size_t>(j)];
      for (const auto& nu : verts) {
        if (seeded >= options_.max_seed_rows) {
          break;
        }
        double best = -std::numeric_limits<double>::infinity();
        bool ok = false;
        for (int c = 0; c < n_cand; ++c) {
          double val = 0.0;
          if (impl.evalSource(src, nu, request.candidates[static_cast<std::size_t>(c)],
                              &val)) {
            best = std::max(best, val);
            ok = true;
          }
        }
        if (!ok) {
          continue;
        }
        SparseVec phi;
        impl.phiRow(tgt, j, nu, &phi);
        appendRow(phi, best);
        ++seeded;
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
      const auto& verts = tgt_geom.region_vertices[static_cast<std::size_t>(j)];

      // Deduplicate first: duplicate rows create dual degeneracy and make the
      // multipliers arbitrary.
      std::vector<const Eigen::VectorXd*> uniq;
      for (const auto& v : verts) {
        bool dup = false;
        for (const auto* u : uniq) {
          if ((*u - v).lpNorm<Eigen::Infinity>() <= kGeomEps) {
            dup = true;
            break;
          }
        }
        if (!dup) {
          uniq.push_back(&v);
        }
      }

      const auto n_v = static_cast<int>(uniq.size());
      const int n_sub_cols = 3 * kSubsystemCount + 1;
      const int zeta_col = n_sub_cols - 1;

      std::vector<int> sstart = {0};
      std::vector<int> sindex;
      std::vector<double> svalue;
      std::vector<double> slower;
      std::vector<double> supper;

      // Group 1: sum_s (c_s . Pr_s(nu) + d_s) + zeta >= sigma * Vt(nu; z).
      for (int i = 0; i < n_v; ++i) {
        const Eigen::VectorXd& nu = *uniq[static_cast<std::size_t>(i)];
        for (int s = 0; s < kSubsystemCount; ++s) {
          const ProjectionLayer& sl = src_geom.layers[s];
          sindex.push_back(3 * s + 0);
          svalue.push_back(nu(sl.axes[0]));
          sindex.push_back(3 * s + 1);
          svalue.push_back(nu(sl.axes[1]));
          sindex.push_back(3 * s + 2);
          svalue.push_back(1.0);
        }
        sindex.push_back(zeta_col);
        svalue.push_back(1.0);
        sstart.push_back(static_cast<int>(sindex.size()));
        slower.push_back(sigma * impl.valueAt(tgt, j, nu, z));
        supper.push_back(kInf);
      }

      // Group 2: c_s . g + d_s <= sigma * T_j^(s)(g).
      std::vector<std::pair<int, int>> group2;  // (plane, pool index)
      for (int s = 0; s < kSubsystemCount; ++s) {
        const int rect = data.region_rect[s][static_cast<std::size_t>(j)];
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
          supper.push_back(sigma * targetValue(j, s, e));
          group2.emplace_back(s, e);
        }
      }

      std::vector<double> sobj(static_cast<std::size_t>(n_sub_cols), 0.0);
      sobj[static_cast<std::size_t>(zeta_col)] = 1.0;
      std::vector<double> sclo(static_cast<std::size_t>(n_sub_cols), -kInf);
      std::vector<double> schi(static_cast<std::size_t>(n_sub_cols), kInf);
      sclo[static_cast<std::size_t>(zeta_col)] = 0.0;

      Highs sub;
      // presolve off: 16 columns make it pure overhead, and it removes any
      // question about postsolve dual recovery, which this cut depends on.
      applyHighsOptions(sub, /*verbose=*/false, /*presolve=*/false);
      if (sub.addCols(n_sub_cols, sobj.data(), sclo.data(), schi.data(), 0,
                      nullptr, nullptr, nullptr)
          != HighsStatus::kOk) {
        throw std::runtime_error("CourierBorderSolver::solve: sub addCols failed");
      }
      if (sub.addRows(static_cast<int>(slower.size()), slower.data(),
                      supper.data(), static_cast<int>(svalue.size()),
                      sstart.data(), sindex.data(), svalue.data())
          != HighsStatus::kOk) {
        throw std::runtime_error("CourierBorderSolver::solve: sub addRows failed");
      }
      if (sub.run() != HighsStatus::kOk
          || sub.getModelStatus() != HighsModelStatus::kOptimal) {
        throw std::runtime_error(
            "CourierBorderSolver::solve: subproblem of region "
            + std::to_string(j) + " not optimal, status "
            + std::to_string(static_cast<int>(sub.getModelStatus())));
      }
      ++local_stats.subproblems_solved;

      const double zeta = sub.getSolution().col_value[static_cast<std::size_t>(zeta_col)];
      if (zeta <= options_.certificate_tol) {
        continue;
      }
      worst = std::max(worst, zeta);

      const auto& ssol = sub.getSolution();
      if (!ssol.dual_valid
          || ssol.row_dual.size() != slower.size()) {
        throw std::runtime_error(
            "CourierBorderSolver::solve: no duals for the region subproblem");
      }

      // lambda on the >= rows, mu on the <= rows. Rather than trust the sign
      // convention, verify the stationarity identities that free c_s, d_s and a
      // basic zeta force: lambda and every mu_s are probability distributions
      // whose barycentres agree in each source plane. A flipped sign fails this
      // immediately instead of silently producing an invalid cut.
      std::vector<double> lambda(static_cast<std::size_t>(n_v));
      for (int i = 0; i < n_v; ++i) {
        lambda[static_cast<std::size_t>(i)] = ssol.row_dual[static_cast<std::size_t>(i)];
      }
      std::vector<double> mu(group2.size());
      for (std::size_t i = 0; i < group2.size(); ++i) {
        mu[i] = -ssol.row_dual[static_cast<std::size_t>(n_v) + i];
      }

      double lam_sum = 0.0;
      for (double v : lambda) {
        if (v < -kDualIdentityTol) {
          throw std::runtime_error(
              "CourierBorderSolver::solve: negative lambda multiplier");
        }
        lam_sum += v;
      }
      if (std::abs(lam_sum - 1.0) > kDualIdentityTol) {
        throw std::runtime_error(
            "CourierBorderSolver::solve: lambda does not sum to one, got "
            + std::to_string(lam_sum));
      }
      std::array<double, kSubsystemCount> mu_sum{};
      for (std::size_t i = 0; i < group2.size(); ++i) {
        if (mu[i] < -kDualIdentityTol) {
          throw std::runtime_error(
              "CourierBorderSolver::solve: negative mu multiplier");
        }
        mu_sum[static_cast<std::size_t>(group2[i].first)] += mu[i];
      }
      for (int s = 0; s < kSubsystemCount; ++s) {
        if (std::abs(mu_sum[static_cast<std::size_t>(s)] - lam_sum)
            > kDualIdentityTol) {
          throw std::runtime_error(
              "CourierBorderSolver::solve: mu of plane " + std::to_string(s)
              + " does not match lambda mass");
        }
      }

      // Cut: sigma * (sum_nu lambda_nu phi_nu)^T z <= sigma * sum mu T.
      SparseVec row;
      for (int i = 0; i < n_v; ++i) {
        if (lambda[static_cast<std::size_t>(i)] == 0.0) {
          continue;
        }
        SparseVec phi;
        impl.phiRow(tgt, j, *uniq[static_cast<std::size_t>(i)], &phi);
        for (std::size_t k = 0; k < phi.cols.size(); ++k) {
          row.add(phi.cols[k],
                  lambda[static_cast<std::size_t>(i)] * phi.vals[k], kGeomEps);
        }
      }
      double rhs = 0.0;
      for (std::size_t i = 0; i < group2.size(); ++i) {
        rhs += mu[i] * targetValue(j, group2[i].first, group2[i].second);
      }

      // A cut that does not separate the current point means the multipliers
      // are wrong and Benders would stall silently.
      const double sep = sigma * (row.dot(z) - rhs);
      if (sep < 0.5 * zeta) {
        throw std::runtime_error(
            "CourierBorderSolver::solve: generated cut does not separate z*, "
            "separation " + std::to_string(sep) + " vs zeta "
            + std::to_string(zeta));
      }
      violations.push_back(Violation{zeta, std::move(row), rhs});
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
  const PhaseGeometry& tgt_geom = (*impl.geometries)[static_cast<std::size_t>(tgt)];
  const bool is_upper = mode == ApproximationMode::Upper;
  const double sigma = is_upper ? -1.0 : 1.0;
  const TargetData& data = impl.target[static_cast<std::size_t>(tgt)];
  const PhaseGeometry& src_geom
      = (*impl.geometries)[static_cast<std::size_t>(request.source_phase)];
  const auto n_cand = static_cast<int>(request.candidates.size());
  const double kInf = std::numeric_limits<double>::infinity();

  std::vector<int> rho_of_region;
  if (!is_upper) {
    rho_of_region.assign(tgt_geom.region_vertices.size(), 0);
    for (std::size_t j = 0; j < tgt_geom.region_vertices.size(); ++j) {
      const auto& verts = tgt_geom.region_vertices[j];
      Eigen::VectorXd centroid = Eigen::VectorXd::Zero(kSpaceDim);
      for (const auto& v : verts) {
        centroid += v;
      }
      centroid /= static_cast<double>(verts.size());
      double best_val = -std::numeric_limits<double>::infinity();
      for (int c = 0; c < n_cand; ++c) {
        double val = 0.0;
        if (impl.evalSource(request.source_phase, centroid,
                            request.candidates[static_cast<std::size_t>(c)],
                            &val)
            && val > best_val) {
          best_val = val;
          rho_of_region[j] = c;
        }
      }
    }
  }

  double worst = 0.0;
  for (int j = 0; j < static_cast<int>(tgt_geom.region_vertices.size()); ++j) {
    const auto& verts = tgt_geom.region_vertices[static_cast<std::size_t>(j)];
    const int n_sub_cols = 3 * kSubsystemCount + 1;
    const int zeta_col = n_sub_cols - 1;

    std::vector<int> sstart = {0};
    std::vector<int> sindex;
    std::vector<double> svalue;
    std::vector<double> slower;
    std::vector<double> supper;

    for (const auto& nu : verts) {
      for (int s = 0; s < kSubsystemCount; ++s) {
        const ProjectionLayer& sl = src_geom.layers[s];
        sindex.push_back(3 * s + 0);
        svalue.push_back(nu(sl.axes[0]));
        sindex.push_back(3 * s + 1);
        svalue.push_back(nu(sl.axes[1]));
        sindex.push_back(3 * s + 2);
        svalue.push_back(1.0);
      }
      sindex.push_back(zeta_col);
      svalue.push_back(1.0);
      sstart.push_back(static_cast<int>(sindex.size()));
      slower.push_back(sigma * impl.valueAt(tgt, j, nu, z));
      supper.push_back(kInf);
    }

    for (int s = 0; s < kSubsystemCount; ++s) {
      const int rect = data.region_rect[s][static_cast<std::size_t>(j)];
      for (int e = data.rect_begin[s][static_cast<std::size_t>(rect)];
           e < data.rect_end[s][static_cast<std::size_t>(rect)]; ++e) {
        const ClipVertex& cv = data.pool[s][static_cast<std::size_t>(e)];
        double t_val;
        if (is_upper) {
          t_val = -kInf;
          for (int c = 0; c < n_cand; ++c) {
            t_val = std::max(t_val, cv.eval(request.candidates[static_cast<std::size_t>(c)]));
          }
        } else {
          // Same candidate solve() would pick. Soundness only needs a courier
          // for SOME rho -- if one exists for rho(j) then
          // Vt <= Phi_{rho(j)} <= max_rho Phi_rho -- so re-verifying against
          // this rho is a valid check, whereas the min over all candidates
          // would demand something the method never claimed.
          t_val = cv.eval(request.candidates[static_cast<std::size_t>(rho_of_region[
              static_cast<std::size_t>(j)])]);
        }
        sindex.push_back(3 * s + 0);
        svalue.push_back(cv.gx);
        sindex.push_back(3 * s + 1);
        svalue.push_back(cv.gy);
        sindex.push_back(3 * s + 2);
        svalue.push_back(1.0);
        sstart.push_back(static_cast<int>(sindex.size()));
        slower.push_back(-kInf);
        supper.push_back(sigma * t_val);
      }
    }

    std::vector<double> sobj(static_cast<std::size_t>(n_sub_cols), 0.0);
    sobj[static_cast<std::size_t>(zeta_col)] = 1.0;
    std::vector<double> sclo(static_cast<std::size_t>(n_sub_cols), -kInf);
    std::vector<double> schi(static_cast<std::size_t>(n_sub_cols), kInf);
    sclo[static_cast<std::size_t>(zeta_col)] = 0.0;

    Highs sub;
    applyHighsOptions(sub, /*verbose=*/false, /*presolve=*/false);
    sub.addCols(n_sub_cols, sobj.data(), sclo.data(), schi.data(), 0, nullptr,
                nullptr, nullptr);
    sub.addRows(static_cast<int>(slower.size()), slower.data(), supper.data(),
                static_cast<int>(svalue.size()), sstart.data(), sindex.data(),
                svalue.data());
    if (sub.run() != HighsStatus::kOk
        || sub.getModelStatus() != HighsModelStatus::kOptimal) {
      throw std::runtime_error(
          "CourierBorderSolver::worstCertificateResidual: subproblem of region "
          + std::to_string(j) + " not optimal");
    }
    worst = std::max(
        worst, sub.getSolution().col_value[static_cast<std::size_t>(zeta_col)]);
  }
  return worst;
}

}  // namespace barycentric_affine_approximator

// NOLINTEND(readability-identifier-naming)
