#include <algo.hpp>
#include <algorithm>
#include <cmath>
#include <vector>
#include "cddwrap/lineareq.hpp"
#include "morph.hpp"
#include "types.hpp"
#include <Eigen/Dense>
#include <Highs.h>

namespace hcpwa {

template <typename T>
static void InplaceUnique(std::vector<T>& data) {
  std::vector<T> found;
  for (const auto& item : data) {
    if (!std::ranges::contains(found, item)) {
      found.push_back(item);
    }
  }
  data.swap(found);
}

std::vector<PolygonResolution> SplitAABBWithLines(AABB<2> aabb,
                                                  const LineSet<2>& lines) {
  std::vector<PolygonResolution> result;
  const std::uint64_t options = 1 << lines.size();
  cddwrap::matrix<double> inequalities(3, lines.size() + 4);
  {
    inequalities[lines.size()][0] = 1;
    inequalities[lines.size()][1] = 0;
    inequalities[lines.size()][2] = (double)aabb.second[0];

    inequalities[lines.size() + 1][0] = 0;
    inequalities[lines.size() + 1][1] = 1;
    inequalities[lines.size() + 1][2] = (double)aabb.second[1];

    inequalities[lines.size() + 2][0] = -1;
    inequalities[lines.size() + 2][1] = 0;
    inequalities[lines.size() + 2][2] = (double)aabb.first[0];

    inequalities[lines.size() + 3][0] = 0;
    inequalities[lines.size() + 3][1] = -1;
    inequalities[lines.size() + 3][2] = (double)aabb.first[1];
  }
  for (std::uint64_t option_ind = 0; option_ind < options; option_ind++) {
    for (std::uint64_t i = 0; i < lines.size(); i++) {
      const bool is_gt = (option_ind & (1 << i)) != 0;
      const Line<2>& line = is_gt ? -lines[i] : lines[i];
      inequalities[i][0] = (double)line[0];
      inequalities[i][1] = (double)line[1];
      inequalities[i][2] = -(double)line[2];
    }

    auto points = GetHullPoints(inequalities);
    if (points.rows() < 3) {
      continue;
    }
    PolygonResolution poly;
    poly.is_gt = option_ind;
    for (std::size_t i = 0; i < points.rows(); i++) {
      poly.polygon.emplace_back(points[i, 0], points[i, 1]);
    }
    result.push_back(std::move(poly));
  }
  return result;
}

void NormalizeVertices(std::vector<hcpwa::PolygonResolution>& data,
                       hcpwa::UniquePool<hcpwa::Vec<2>>& pool) {
  for (auto& poly : data) {
    for (auto& point : poly.polygon) {
      point = pool.Unique(point);
    }
    InplaceUnique(poly.polygon);
  }
  for (std::size_t i = 0; i < data.size();) {
    if (data[i].polygon.size() < 3) {
      std::swap(data[i], data.back());
      data.pop_back();
    } else {
      i++;
    }
  }
}

std::vector<std::pair<hcpwa::Triangle, std::size_t>> Triangulate(
    const std::span<hcpwa::PolygonResolution>& data) {
  using IndexedTriangle = std::pair<Triangle, std::size_t>;
  std::vector<IndexedTriangle> result;

  for (std::size_t i = 0; i < data.size(); i++) {
    const auto& poly = data[i].polygon;
    if (poly.size() < 3) {
      continue;
    }
    // The vertices arrive from LinesToPoints, which reads them out of cddlib's
    // extreme-point enumeration -- a set, with no boundary order. A fan needs
    // one. Without it the fan's triangles overlap each other and leave the rest
    // of the cell uncovered, and because the caller sums |area| the two errors
    // partly cancel and the total looks plausible. Measured before this fix, on
    // the N=160 production geometry: every projection layer covered between
    // 80.6% and 98.1% of [0,N]^2 with 1.9% to 9.3% covered twice, while the
    // summed |area| sat within 1% of N^2 on three layers out of ten and was
    // exact on two of them.
    //
    // These cells are intersections of half-planes, hence convex, so sorting
    // the vertices by angle about the centroid is exactly the boundary order
    // and the fan from any one of them is then a true triangulation.
    std::vector<hcpwa::Vec<2>> ordered(poly.begin(), poly.end());
    double cx = 0.0;
    double cy = 0.0;
    for (const auto& p : ordered) {
      cx += static_cast<double>(p[0]);
      cy += static_cast<double>(p[1]);
    }
    cx /= static_cast<double>(ordered.size());
    cy /= static_cast<double>(ordered.size());
    std::sort(ordered.begin(), ordered.end(),
              [cx, cy](const hcpwa::Vec<2>& p, const hcpwa::Vec<2>& q) {
                return std::atan2(static_cast<double>(p[1]) - cy,
                                  static_cast<double>(p[0]) - cx)
                       < std::atan2(static_cast<double>(q[1]) - cy,
                                    static_cast<double>(q[0]) - cx);
              });
    const auto& a = ordered[0];
    for (std::size_t j = 2; j < ordered.size(); j++) {
      result.emplace_back(
          Triangle{.a = a, .b = ordered[j - 1], .c = ordered[j]}, i);
    }
  }
  return result;
}

hcpwa::LineSet<8> CalcPrism(const Triangle& triangle,
                            const std::array<int, 2>& dims) {
  std::vector<hcpwa::Line<8>> result;
  {
    hcpwa::Line<2> line
        = cross(hcpwa::Vec<3>{triangle.a, 1}, hcpwa::Vec<3>{triangle.b, 1});
    if (dot(line, hcpwa::Vec<3>{triangle.c, 1}) > 0) {
      line = -line;
    }
    result.push_back(hcpwa::Promote<9>(line, {dims[0], dims[1], 8}));
  }
  {
    hcpwa::Line<2> line
        = cross(hcpwa::Vec<3>{triangle.b, 1}, hcpwa::Vec<3>{triangle.c, 1});
    if (dot(line, hcpwa::Vec<3>{triangle.a, 1}) > 0) {
      line = -line;
    }
    result.push_back(hcpwa::Promote<9>(line, {dims[0], dims[1], 8}));
  }
  {
    hcpwa::Line<2> line
        = cross(hcpwa::Vec<3>{triangle.a, 1}, hcpwa::Vec<3>{triangle.c, 1});
    if (dot(line, hcpwa::Vec<3>{triangle.b, 1}) > 0) {
      line = -line;
    }
    result.push_back(hcpwa::Promote<9>(line, {dims[0], dims[1], 8}));
  }
  return result;
}

hcpwa::LineSet<8> CalcPrism(const Polygon& polygon,
                            const std::array<int, 2>& dims) {
  std::vector<hcpwa::Line<8>> result;
  for (int i = 0; i < polygon.size(); i++) {
    auto a = polygon[i];
    auto b = polygon[(i + 1) % polygon.size()];
    auto c = polygon[(i + 2) % polygon.size()];
    hcpwa::Line<2> line = cross(hcpwa::Vec<3>{a, 1}, hcpwa::Vec<3>{b, 1});
    if (dot(line, hcpwa::Vec<3>{c, 1}) > 0) {
      line = -line;
    }
    result.push_back(hcpwa::Promote<9>(line, {dims[0], dims[1], 8}));
  }
  return result;
}

template <int N>
std::vector<hcpwa::Vec<N>> LinesToPoints(const hcpwa::LineSet<N>& data) {
  cddwrap::matrix<double> inequalities(N + 1, data.size());
  for (int i = 0; i < data.size(); i++) {
    for (int j = 0; j < N; j++) {
      inequalities[i][j] = (double)data[i][j];
    }
    inequalities[i][N] = -(double)data[i][N];
  }

  auto points = GetHullPoints(inequalities);
  if (points.rows() < N) {
    return {};
  }
  std::vector<hcpwa::Vec<N>> result;
  for (std::size_t i = 0; i < points.rows(); i++) {
    hcpwa::Vec<N> point = kZeroVec;
    for (int n = 0; n < N; n++) {
      point[n] = points[i, n];
    }
    // result.emplace_back(points[i, 0], points[i, 1], points[i, 2],
    // points[i, 3],
    //                     points[i, 4], points[i, 5], points[i, 6],
    //                     points[i, 7]);
    result.push_back(point);
  }
  return result;
}

template std::vector<hcpwa::Vec<8>> LinesToPoints(
    const hcpwa::LineSet<8>& data);
template std::vector<hcpwa::Vec<3>> LinesToPoints(
    const hcpwa::LineSet<3>& data);
template std::vector<hcpwa::Vec<2>> LinesToPoints(
    const hcpwa::LineSet<2>& data);

std::vector<hcpwa::TriangleWithUniqueVertices> GetTrianglesWithUniqueVertices(
    const AABB<2>& aabb, std::vector<hcpwa::PolygonResolution>& polygons) {
  using Comp = hcpwa::FixedPrecisionComparator<hcpwa::Vec<2>, 0.01f>;
  hcpwa::UniquePool<hcpwa::Vec<2>> pool(Comp{});
  pool.Unique(aabb.first[0, 1]);
  pool.Unique(aabb.second[0, 1]);
  pool.Unique({aabb.first[0], aabb.second[1]});
  pool.Unique({aabb.second[0], aabb.first[1]});
  hcpwa::NormalizeVertices(polygons, pool);
  auto triangles = Triangulate(polygons);
  std::vector<hcpwa::TriangleWithUniqueVertices> result;
  for (std::size_t i = 0; i < triangles.size(); i++) {
    result.push_back(hcpwa::TriangleWithUniqueVertices{
        triangles[i].first, pool.Index(triangles[i].first.a),
        pool.Index(triangles[i].first.b), pool.Index(triangles[i].first.c),
        triangles[i].second});
  }
  return result;
}

}  // namespace hcpwa