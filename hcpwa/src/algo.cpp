#include <algo.hpp>
#include <algorithm>
#include "lineareq.hpp"
#include "matrix.hpp"

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
  matrix<double> inequalities(3, lines.size() + 4);
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
      inequalities[i][2] = -(double)(line.z);
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
}  // namespace hcpwa