#include <algo.hpp>
#include "lineareq.hpp"
#include "matrix.hpp"

namespace hcpwa {
std::vector<PolygonResolution> SplitAABBWithLines(AABB<2> aabb,
                                                  const LineSet<2>& lines) {
  std::vector<PolygonResolution> result;
  const std::uint64_t options = 1 << lines.size();
  matrix<double> inequalities(3, lines.size() + 4);
  {
    inequalities[lines.size()][0] = 1;
    inequalities[lines.size()][1] = 0;
    inequalities[lines.size()][2] = aabb.second.x;

    inequalities[lines.size() + 1][0] = 0;
    inequalities[lines.size() + 1][1] = 1;
    inequalities[lines.size() + 1][2] = aabb.second.y;

    inequalities[lines.size() + 2][0] = -1;
    inequalities[lines.size() + 2][1] = 0;
    inequalities[lines.size() + 2][2] = aabb.first.x;

    inequalities[lines.size() + 3][0] = 0;
    inequalities[lines.size() + 3][1] = -1;
    inequalities[lines.size() + 3][2] = aabb.first.y;
  }
  for (std::uint64_t option_ind = 0; option_ind < options; option_ind++) {
    for (std::uint64_t i = 0; i < lines.size(); i++) {
      const bool is_gt = (option_ind & (1 << i)) != 0;
      const Line<2>& line = is_gt ? -lines[i] : lines[i];
      inequalities[i][0] = line.x;
      inequalities[i][1] = line.y;
      inequalities[i][2] = -line.z;
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
}  // namespace hcpwa