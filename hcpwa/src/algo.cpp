#include <algo.hpp>
#include <algorithm>
#include "lineareq.hpp"
#include "matrix.hpp"
#include "morph.hpp"

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
    inequalities[lines.size()][2] = (double)aabb.second.x;

    inequalities[lines.size() + 1][0] = 0;
    inequalities[lines.size() + 1][1] = 1;
    inequalities[lines.size() + 1][2] = (double)aabb.second.y;

    inequalities[lines.size() + 2][0] = -1;
    inequalities[lines.size() + 2][1] = 0;
    inequalities[lines.size() + 2][2] = (double)aabb.first.x;

    inequalities[lines.size() + 3][0] = 0;
    inequalities[lines.size() + 3][1] = -1;
    inequalities[lines.size() + 3][2] = (double)aabb.first.y;
  }
  for (std::uint64_t option_ind = 0; option_ind < options; option_ind++) {
    for (std::uint64_t i = 0; i < lines.size(); i++) {
      const bool is_gt = (option_ind & (1 << i)) != 0;
      const Line<2>& line = is_gt ? -lines[i] : lines[i];
      inequalities[i][0] = (double)line.x;
      inequalities[i][1] = (double)line.y;
      inequalities[i][2] = -(double)line.z;
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
    const auto& a = poly[0];
    for (std::size_t j = 2; j < poly.size(); j++) {
      result.emplace_back(Triangle{.a = a, .b = poly[j - 1], .c = poly[j]}, i);
    }
  }
  return result;
}

std::vector<hcpwa::Line<8>> CalcPrism(const Triangle& triangle,
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

}  // namespace hcpwa