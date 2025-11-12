#include <algo.hpp>
#include <symbolic.hpp>
#include <morph.hpp>
#include "cddwrap/cdd.hpp"
#include "utility.hpp"

namespace hcpwa {

// NOLINTNEXTLINE
using namespace hcpwa::symbols;

AreasVerticesResult compute_areas_vertices(int phase, float N, float F, float v, float w,
                                           float b51, float b57, float b84, float b86,
                                           float b31, float b36, float b24, float b27,
                                           float f2min, float f3min, float f5min, float f8min,
                                           float f2max, float f3max, float f5max, float f8max) {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;

  hcpwa::AABB<8> aabb = {{0, 0, 0, 0, 0, 0, 0, 0}, {N, N, N, N, N, N, N, N}};
  hcpwa::AABB<2> aabb2d = {{0, 0}, {N, N}};

  constexpr auto n1 = X<0>{};
  constexpr auto n2 = X<1>{};
  constexpr auto n3 = X<2>{};
  constexpr auto n4 = X<3>{};
  constexpr auto n5 = X<4>{};
  constexpr auto n6 = X<5>{};
  constexpr auto n7 = X<6>{};
  constexpr auto n8 = X<7>{};

  auto f31 = SymMin(b31 * F, v * b31 * n3, w * (N - n1));
  auto f36 = SymMin(b36 * F, v * b36 * n3, w * (N - n6));
  auto f24 = SymMin(b24 * F, v * b24 * n2, w * (N - n4));
  auto f27 = SymMin(b27 * F, v * b27 * n2, w * (N - n7));
  auto f51 = SymMin(b51 * F, v * b51 * n5, w * (N - n1));
  auto f57 = SymMin(b57 * F, v * b57 * n5, w * (N - n7));
  auto f84 = SymMin(b84 * F, v * b84 * n8, w * (N - n4));
  auto f86 = SymMin(b86 * F, v * b86 * n8, w * (N - n6));
  // f1out(t) = min{vn1(τ),F}
  auto f1out = SymMin(v * n1, F);
  // f4out(t) = min{vn4(τ),F}
  auto f4out = SymMin(v * n4, F);
  // f6out(t) = min{vn6(τ),F}
  auto f6out = SymMin(v * n6, F);
  // f7out(t) = min{vn7(τ),F}
  auto f7out = SymMin(v * n7, F);
  // f_i_in(t) = min{f_i_min,w(N−n_i(τ))} for i = 2, 3, 5, 8
  auto f2in = SymMin(f2min, w * (N - n2));
  auto f3in = SymMin(f3min, w * (N - n3));
  auto f5in = SymMin(f5min, w * (N - n5));
  auto f8in = SymMin(f8min, w * (N - n8));

  // Plane (5, 1)
  hcpwa::LineSet<8> lines51;
  {
    if (phase == 0) {
      auto resolutions = MinResolutions(f31);
      constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
      auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
      lines51.append_range(lines);
    } else {
      auto resolutions = MinResolutions(f51);
      constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
      auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
      lines51.append_range(lines);
    }
  }
  {
    auto resolutions = MinResolutions(f1out);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines51.append_range(lines);
  }
  {
    auto resolutions = MinResolutions(f5in);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines51.append_range(lines);
  }
  // Plane (5, 7)
  hcpwa::LineSet<8> lines57;
  {
    if (phase == 0) {
      auto resolutions = MinResolutions(f36);
      constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
      auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
      lines57.append_range(lines);
    } else {
      auto resolutions = MinResolutions(f57);
      constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
      auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
      lines57.append_range(lines);
    }
  }
  {
    auto resolutions = MinResolutions(f5in);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines57.append_range(lines);
  }
  {
    auto resolutions = MinResolutions(f7out);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines57.append_range(lines);
  }
  // Plane (8, 4)
  hcpwa::LineSet<8> lines84;
  {
    if (phase == 0) {
      auto resolutions = MinResolutions(f24);
      constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
      auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
      lines84.append_range(lines);
    } else {
      auto resolutions = MinResolutions(f84);
      constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
      auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
      lines84.append_range(lines);
    }
  }
  {
    auto resolutions = MinResolutions(f4out);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines84.append_range(lines);
  }
  {
    auto resolutions = MinResolutions(f8in);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines84.append_range(lines);
  }
  // Plane (8, 6)
  hcpwa::LineSet<8> lines86;
  {
    if (phase == 0) {
      auto resolutions = MinResolutions(f27);
      constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
      auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
      lines86.append_range(lines);
    } else {
      auto resolutions = MinResolutions(f86);
      constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
      auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
      lines86.append_range(lines);
    }
  }
  {
    auto resolutions = MinResolutions(f6out);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines86.append_range(lines);
  }
  {
    auto resolutions = MinResolutions(f8in);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines86.append_range(lines);
  }

  auto polygon51 = hcpwa::SplitAABBWithLines(aabb2d, hcpwa::DimensionCast<2, 8>(lines51, {1 - 1, 5 - 1}));
  auto polygon57 = hcpwa::SplitAABBWithLines(aabb2d, hcpwa::DimensionCast<2, 8>(lines57, {5 - 1, 7 - 1}));
  auto polygon84 = hcpwa::SplitAABBWithLines(aabb2d, hcpwa::DimensionCast<2, 8>(lines84, {4 - 1, 8 - 1}));
  auto polygon86 = hcpwa::SplitAABBWithLines(aabb2d, hcpwa::DimensionCast<2, 8>(lines86, {6 - 1, 8 - 1}));

  // Triangulated polygons with unique vertices
  auto triangles51 = hcpwa::GetTrianglesWithUniqueVertices(aabb2d, polygon51);
  auto triangles57 = hcpwa::GetTrianglesWithUniqueVertices(aabb2d, polygon57);
  auto triangles84 = hcpwa::GetTrianglesWithUniqueVertices(aabb2d, polygon84);
  auto triangles86 = hcpwa::GetTrianglesWithUniqueVertices(aabb2d, polygon86);

  // Calculate prisms
  std::vector<hcpwa::LineSet<8>> prisms51;
  std::vector<hcpwa::LineSet<8>> prisms57;
  std::vector<hcpwa::LineSet<8>> prisms84;
  std::vector<hcpwa::LineSet<8>> prisms86;
  for (auto& triangle : triangles51) {
    auto prism = hcpwa::CalcPrism(triangle, {1 - 1, 5 - 1});
    prisms51.push_back(prism);
  }
  for (auto& triangle : triangles57) {
    auto prism = hcpwa::CalcPrism(triangle, {5 - 1, 7 - 1});
    prisms57.push_back(prism);
  }
  for (auto& triangle : triangles84) {
    auto prism = hcpwa::CalcPrism(triangle, {4 - 1, 8 - 1});
    prisms84.push_back(prism);
  }
  for (auto& triangle : triangles86) {
    auto prism = hcpwa::CalcPrism(triangle, {6 - 1, 8 - 1});
    prisms86.push_back(prism);
  }

  // Intersect all prisms
  std::vector<std::vector<size_t>> intersection_prism_indices;
  std::vector<std::vector<hcpwa::Vec<8>>> intersection_points;
  for (size_t prism51_idx = 0; prism51_idx < prisms51.size(); prism51_idx++) {
    for (size_t prism57_idx = 0; prism57_idx < prisms57.size(); prism57_idx++) {
      for (size_t prism84_idx = 0; prism84_idx < prisms84.size(); prism84_idx++) {
        for (size_t prism86_idx = 0; prism86_idx < prisms86.size(); prism86_idx++) {
          hcpwa::LineSet<8> concatenated_prisms;
          auto &prism51 = prisms51[prism51_idx];
          auto &prism57 = prisms57[prism57_idx];
          auto &prism84 = prisms84[prism84_idx];
          auto &prism86 = prisms86[prism86_idx];
          concatenated_prisms.append_range(prism51);
          concatenated_prisms.append_range(prism57);
          concatenated_prisms.append_range(prism84);
          concatenated_prisms.append_range(prism86);
          auto intersection = hcpwa::LinesToPoints(concatenated_prisms);
          if (intersection.size() > 0) {
            intersection_points.push_back(intersection);
            // Store the indices of the prisms that form the intersection
            std::vector<size_t> prism_indices = {prism51_idx, prism57_idx, prism84_idx, prism86_idx};
            intersection_prism_indices.push_back(prism_indices);
          }
        }
      }
    }
  }
  
  AreasVerticesResult result;
  result.triangles51 = triangles51;
  result.triangles57 = triangles57;
  result.triangles84 = triangles84;
  result.triangles86 = triangles86;
  result.intersection_points = intersection_points;
  result.intersection_prism_indices = intersection_prism_indices;
  return result;
}

}  // namespace hcpwa

