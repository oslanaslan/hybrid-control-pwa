#include <algo.hpp>
#include <symbolic.hpp>
#include <morph.hpp>
#include "cddwrap/cdd.hpp"
#include "utility.hpp"

namespace hcpwa {

// NOLINTNEXTLINE
using namespace hcpwa::symbols;

AreasVerticesResult compute_areas_vertices(float N, float F, float v, float w,
                                           float b51, float b57, float b84, float b86,
                                           float b31, float b36, float b24, float b27,
                                           float f2min, float f3min, float f5min, float f8min,
                                           float f2max, float f3max, float f5max, float f8max) {
  // cddwrap::global_init();
  // defer _ = &cddwrap::global_free;

  hcpwa::AABB<8> aabb = {{0, 0, 0, 0, 0, 0, 0, 0}, {N, N, N, N, N, N, N, N}};
  hcpwa::AABB<2> aabb2d = {{0, 0}, {N, N}};
  const auto aabb_bounds = hcpwa::AABBBounds(aabb);

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
  // f_out Lower bounds
  auto f1out_lower = SymMin(v * n1, F);
  auto f4out_lower = SymMin(v * n4, F);
  auto f6out_lower = SymMin(v * n6, F);
  auto f7out_lower = SymMin(v * n7, F);
  // f_in Lower bounds
  auto f2in_lower = SymMin(f2min, w * (N - n2));
  auto f3in_lower = SymMin(f3min, w * (N - n3));
  auto f5in_lower = SymMin(f5min, w * (N - n5));
  auto f8in_lower = SymMin(f8min, w * (N - n8));
  // f_in Upper bounds
  auto f2in_upper = SymMin(f2max, w * (N - n2));
  auto f3in_upper = SymMin(f3max, w * (N - n3));
  auto f5in_upper = SymMin(f5max, w * (N - n5));
  auto f8in_upper = SymMin(f8max, w * (N - n8));

  // Phase 0
  // Plane (3, 1)
  hcpwa::LineSet<8> lines31;
  // Plane (3, 6)
  hcpwa::LineSet<8> lines36;
  // Plane (2, 4)
  hcpwa::LineSet<8> lines24;
  // Plane (2, 7)
  hcpwa::LineSet<8> lines27;
  // Plane (5, 8)
  hcpwa::LineSet<8> lines58;
  // Phase 1
  // Plane (5, 1)
  hcpwa::LineSet<8> lines51;
  // Plane (5, 7)
  hcpwa::LineSet<8> lines57;
  // Plane (8, 4)
  hcpwa::LineSet<8> lines84;
  // Plane (8, 6)
  hcpwa::LineSet<8> lines86;
  // Plane (2, 3)
  hcpwa::LineSet<8> lines23;
  {
    auto resolutions = MinResolutions(f51);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines51.append_range(lines);
  }
  {
    auto resolutions = MinResolutions(f57);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines57.append_range(lines);
  }
  {
    auto resolutions = MinResolutions(f84);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines84.append_range(lines);
  }
  {
    auto resolutions = MinResolutions(f86);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines86.append_range(lines);
  }
  {
      auto resolutions = MinResolutions(f31);
      constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
      auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
      lines31.append_range(lines);
  }
  {
      auto resolutions = MinResolutions(f36);
      constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
      auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
      lines36.append_range(lines);
  }
  {
      auto resolutions = MinResolutions(f24);
      constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
      auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
      lines24.append_range(lines);
  }
  {
      auto resolutions = MinResolutions(f27);
      constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
      auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
      lines27.append_range(lines);
  }
  // In planes lower bounds
  {
    auto resolutions = MinResolutions(f2in_lower);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines24.append_range(lines);
    lines27.append_range(lines);
    lines23.append_range(lines);
  }
  {
    auto resolutions = MinResolutions(f3in_lower);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines31.append_range(lines);
    lines36.append_range(lines);
    lines23.append_range(lines);
  }
  {
    auto resolutions = MinResolutions(f5in_lower);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines51.append_range(lines);
    lines57.append_range(lines);
    lines58.append_range(lines);
  }
  {
    auto resolutions = MinResolutions(f8in_lower);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines84.append_range(lines);
    lines86.append_range(lines);
    lines58.append_range(lines);
  }
  // In planes upper bounds
  {
    auto resolutions = MinResolutions(f2in_upper);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines24.append_range(lines);
    lines27.append_range(lines);
    lines23.append_range(lines);
  }
  {
    auto resolutions = MinResolutions(f3in_upper);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines31.append_range(lines);
    lines36.append_range(lines);
    lines23.append_range(lines);
  }
  {
    auto resolutions = MinResolutions(f5in_upper);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines51.append_range(lines);
    lines57.append_range(lines);
    lines58.append_range(lines);
  }
  {
    auto resolutions = MinResolutions(f8in_upper);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines84.append_range(lines);
    lines86.append_range(lines);
    lines58.append_range(lines);
  }
  // Out planes
  {
    auto resolutions = MinResolutions(f1out_lower);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines51.append_range(lines);
    lines31.append_range(lines);
  }
  {
    auto resolutions = MinResolutions(f4out_lower);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines84.append_range(lines);
    lines24.append_range(lines);
  }
  {
    auto resolutions = MinResolutions(f6out_lower);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines36.append_range(lines);
    lines86.append_range(lines);
  }
  {
    auto resolutions = MinResolutions(f7out_lower);
    constexpr int kDim = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);
    lines27.append_range(lines);
    lines57.append_range(lines);
  }

  // DEBUG: print lines counts
  std::cout << "Lines51 count: " << lines51.size() << "\n";
  std::cout << "Lines57 count: " << lines57.size() << "\n";
  std::cout << "Lines84 count: " << lines84.size() << "\n";
  std::cout << "Lines86 count: " << lines86.size() << "\n";
  std::cout << "Lines58 count: " << lines58.size() << "\n";
  std::cout << "Lines31 count: " << lines31.size() << "\n";
  std::cout << "Lines36 count: " << lines36.size() << "\n";
  std::cout << "Lines24 count: " << lines24.size() << "\n";
  std::cout << "Lines27 count: " << lines27.size() << "\n";
  std::cout << "Lines23 count: " << lines23.size() << "\n";

  // DEBUG: print lines coefficients
  std::cout << "Lines51 coefficients: " << std::format("{}", lines51) << "\n";
  std::cout << "Lines57 coefficients: " << std::format("{}", lines57) << "\n";
  std::cout << "Lines84 coefficients: " << std::format("{}", lines84) << "\n";
  std::cout << "Lines86 coefficients: " << std::format("{}", lines86) << "\n";
  std::cout << "Lines58 coefficients: " << std::format("{}", lines58) << "\n";
  std::cout << "Lines31 coefficients: " << std::format("{}", lines31) << "\n";
  std::cout << "Lines36 coefficients: " << std::format("{}", lines36) << "\n";
  std::cout << "Lines24 coefficients: " << std::format("{}", lines24) << "\n";
  std::cout << "Lines27 coefficients: " << std::format("{}", lines27) << "\n";
  std::cout << "Lines23 coefficients: " << std::format("{}", lines23) << "\n";

  auto polygon51 = hcpwa::SplitAABBWithLines(aabb2d, hcpwa::DimensionCast<2, 8>(lines51, {1 - 1, 5 - 1}));
  auto polygon57 = hcpwa::SplitAABBWithLines(aabb2d, hcpwa::DimensionCast<2, 8>(lines57, {5 - 1, 7 - 1}));
  auto polygon84 = hcpwa::SplitAABBWithLines(aabb2d, hcpwa::DimensionCast<2, 8>(lines84, {4 - 1, 8 - 1}));
  auto polygon86 = hcpwa::SplitAABBWithLines(aabb2d, hcpwa::DimensionCast<2, 8>(lines86, {6 - 1, 8 - 1}));
  auto polygon58 = hcpwa::SplitAABBWithLines(aabb2d, hcpwa::DimensionCast<2, 8>(lines58, {5 - 1, 8 - 1}));
  auto polygon31 = hcpwa::SplitAABBWithLines(aabb2d, hcpwa::DimensionCast<2, 8>(lines31, {1 - 1, 3 - 1}));
  auto polygon36 = hcpwa::SplitAABBWithLines(aabb2d, hcpwa::DimensionCast<2, 8>(lines36, {3 - 1, 6 - 1}));
  auto polygon24 = hcpwa::SplitAABBWithLines(aabb2d, hcpwa::DimensionCast<2, 8>(lines24, {2 - 1, 4 - 1}));
  auto polygon27 = hcpwa::SplitAABBWithLines(aabb2d, hcpwa::DimensionCast<2, 8>(lines27, {2 - 1, 7 - 1}));
  auto polygon23 = hcpwa::SplitAABBWithLines(aabb2d, hcpwa::DimensionCast<2, 8>(lines23, {2 - 1, 3 - 1}));

  // DEBUG: print polygons counts
  std::cout << "Polygon51 count: " << polygon51.size() << "\n";
  std::cout << "Polygon57 count: " << polygon57.size() << "\n";
  std::cout << "Polygon84 count: " << polygon84.size() << "\n";
  std::cout << "Polygon86 count: " << polygon86.size() << "\n";
  std::cout << "Polygon58 count: " << polygon58.size() << "\n";
  std::cout << "Polygon31 count: " << polygon31.size() << "\n";
  std::cout << "Polygon36 count: " << polygon36.size() << "\n";
  std::cout << "Polygon24 count: " << polygon24.size() << "\n";
  std::cout << "Polygon27 count: " << polygon27.size() << "\n";
  std::cout << "Polygon23 count: " << polygon23.size() << "\n";

  // DEBUG: print polygons
  // std::cout << "Polygons51: " << "\n";
  for (int i = 0; i < polygon51.size(); i++) {
    std::cout << "polygon " << i << ": ";
    for (int j = 0; j < polygon51[i].polygon.size(); j++) {
      std::cout << polygon51[i].polygon[j] << " ";
    }
    std::cout << "\n";
  }
  // std::cout << "Polygons57: " << "\n";
  for (int i = 0; i < polygon57.size(); i++) {
    std::cout << "polygon " << i << ": ";
    for (int j = 0; j < polygon57[i].polygon.size(); j++) {
      std::cout << polygon57[i].polygon[j] << " ";
    }
    std::cout << "\n";
  }
  // std::cout << "Polygons84: " << "\n";
  for (int i = 0; i < polygon84.size(); i++) {
    std::cout << "polygon " << i << ": ";
    for (int j = 0; j < polygon84[i].polygon.size(); j++) {
      std::cout << polygon84[i].polygon[j] << " ";
    }
    std::cout << "\n";
  }
  // std::cout << "Polygons86: " << "\n";
  for (int i = 0; i < polygon86.size(); i++) {
    std::cout << "polygon " << i << ": ";
    for (int j = 0; j < polygon86[i].polygon.size(); j++) {
      std::cout << polygon86[i].polygon[j] << " ";
    }
    std::cout << "\n";
  }
  // std::cout << "Polygons58: " << "\n";
  for (int i = 0; i < polygon58.size(); i++) {
    std::cout << "polygon " << i << ": ";
    for (int j = 0; j < polygon58[i].polygon.size(); j++) {
      std::cout << polygon58[i].polygon[j] << " ";
    }
    std::cout << "\n";
  }
  // std::cout << "Polygons31: " << "\n";
  for (int i = 0; i < polygon31.size(); i++) {
    std::cout << "polygon " << i << ": ";
    for (int j = 0; j < polygon31[i].polygon.size(); j++) {
      std::cout << polygon31[i].polygon[j] << " ";
    }
    std::cout << "\n";
  }
  // std::cout << "Polygons36: " << "\n";
  for (int i = 0; i < polygon36.size(); i++) {
    std::cout << "polygon " << i << ": ";
    for (int j = 0; j < polygon36[i].polygon.size(); j++) {
      std::cout << polygon36[i].polygon[j] << " ";
    }
    std::cout << "\n";
  }
  // std::cout << "Polygons24: " << "\n";
  for (int i = 0; i < polygon24.size(); i++) {
    std::cout << "polygon " << i << ": ";
    for (int j = 0; j < polygon24[i].polygon.size(); j++) {
      std::cout << polygon24[i].polygon[j] << " ";
    }
    std::cout << "\n";
  }
  // std::cout << "Polygons27: " << "\n";
  for (int i = 0; i < polygon27.size(); i++) {
    std::cout << "polygon " << i << ": ";
    for (int j = 0; j < polygon27[i].polygon.size(); j++) {
      std::cout << polygon27[i].polygon[j] << " ";
    }
    std::cout << "\n";
  }
  // std::cout << "Polygons23: " << "\n";
  for (int i = 0; i < polygon23.size(); i++) {
    std::cout << "polygon " << i << ": ";
    for (int j = 0; j < polygon23[i].polygon.size(); j++) {
      std::cout << polygon23[i].polygon[j] << " ";
    }
    std::cout << "\n";
  }

  std::cout << "Triangulated polygons with unique vertices: " << "\n";
  // Triangulated polygons with unique vertices
  auto triangles51 = hcpwa::GetTrianglesWithUniqueVertices(aabb2d, polygon51);
  auto triangles57 = hcpwa::GetTrianglesWithUniqueVertices(aabb2d, polygon57);
  auto triangles84 = hcpwa::GetTrianglesWithUniqueVertices(aabb2d, polygon84);
  auto triangles86 = hcpwa::GetTrianglesWithUniqueVertices(aabb2d, polygon86);
  auto triangles58 = hcpwa::GetTrianglesWithUniqueVertices(aabb2d, polygon58);
  auto triangles31 = hcpwa::GetTrianglesWithUniqueVertices(aabb2d, polygon31);
  auto triangles36 = hcpwa::GetTrianglesWithUniqueVertices(aabb2d, polygon36);
  auto triangles24 = hcpwa::GetTrianglesWithUniqueVertices(aabb2d, polygon24);
  auto triangles27 = hcpwa::GetTrianglesWithUniqueVertices(aabb2d, polygon27);
  auto triangles23 = hcpwa::GetTrianglesWithUniqueVertices(aabb2d, polygon23);

  std::cout << "Triangles count: " << triangles51.size() << " " << triangles57.size() << " " << triangles84.size() << " " << triangles86.size() << " " << triangles58.size() << " " << triangles31.size() << " " << triangles36.size() << " " << triangles24.size() << " " << triangles27.size() << " " << triangles23.size() << "\n";
  for (int i = 0; i < triangles51.size(); i++) {
    std::cout << "triangle " << i << ": " << triangles51[i].a << " " << triangles51[i].b << " " << triangles51[i].c << "\n";
  }
  for (int i = 0; i < triangles57.size(); i++) {
    std::cout << "triangle " << i << ": " << triangles57[i].a << " " << triangles57[i].b << " " << triangles57[i].c << "\n";
  }
  for (int i = 0; i < triangles84.size(); i++) {
    std::cout << "triangle " << i << ": " << triangles84[i].a << " " << triangles84[i].b << " " << triangles84[i].c << "\n";
  }
  for (int i = 0; i < triangles86.size(); i++) {
    std::cout << "triangle " << i << ": " << triangles86[i].a << " " << triangles86[i].b << " " << triangles86[i].c << "\n";
  }
  for (int i = 0; i < triangles58.size(); i++) {
    std::cout << "triangle " << i << ": " << triangles58[i].a << " " << triangles58[i].b << " " << triangles58[i].c << "\n";
  }
  for (int i = 0; i < triangles31.size(); i++) {
    std::cout << "triangle " << i << ": " << triangles31[i].a << " " << triangles31[i].b << " " << triangles31[i].c << "\n";
  }
  for (int i = 0; i < triangles36.size(); i++) {
    std::cout << "triangle " << i << ": " << triangles36[i].a << " " << triangles36[i].b << " " << triangles36[i].c << "\n";
  }
  for (int i = 0; i < triangles24.size(); i++) {
    std::cout << "triangle " << i << ": " << triangles24[i].a << " " << triangles24[i].b << " " << triangles24[i].c << "\n";
  }
  for (int i = 0; i < triangles27.size(); i++) {
    std::cout << "triangle " << i << ": " << triangles27[i].a << " " << triangles27[i].b << " " << triangles27[i].c << "\n";
  }
  for (int i = 0; i < triangles23.size(); i++) {
    std::cout << "triangle " << i << ": " << triangles23[i].a << " " << triangles23[i].b << " " << triangles23[i].c << "\n";
  }
  // Calculate prisms
  std::vector<hcpwa::LineSet<8>> prisms51;
  std::vector<hcpwa::LineSet<8>> prisms57;
  std::vector<hcpwa::LineSet<8>> prisms84;
  std::vector<hcpwa::LineSet<8>> prisms86;
  std::vector<hcpwa::LineSet<8>> prisms58;
  std::vector<hcpwa::LineSet<8>> prisms31;
  std::vector<hcpwa::LineSet<8>> prisms36;
  std::vector<hcpwa::LineSet<8>> prisms24;
  std::vector<hcpwa::LineSet<8>> prisms27;
  std::vector<hcpwa::LineSet<8>> prisms23;
  // Phase 0
  for (auto& triangle : triangles31) {
    auto prism = hcpwa::CalcPrism(triangle, {1 - 1, 3 - 1});
    prisms31.push_back(prism);
  }
  for (auto& triangle : triangles36) {
    auto prism = hcpwa::CalcPrism(triangle, {3 - 1, 6 - 1});
    prisms36.push_back(prism);
  }
  for (auto& triangle : triangles24) {
    auto prism = hcpwa::CalcPrism(triangle, {2 - 1, 4 - 1});
    prisms24.push_back(prism);
  }
  for (auto& triangle : triangles27) {
    auto prism = hcpwa::CalcPrism(triangle, {2 - 1, 7 - 1});
    prisms27.push_back(prism);
  }
  for (auto& triangle : triangles58) {
    auto prism = hcpwa::CalcPrism(triangle, {5 - 1, 8 - 1});
    prisms58.push_back(prism);
  }
  // Phase 1
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
  for (auto& triangle : triangles23) {
    auto prism = hcpwa::CalcPrism(triangle, {2 - 1, 3 - 1});
    prisms23.push_back(prism);
  }
  // Helper lambda to compute intersections for a set of prism collections
  auto compute_intersections = [](
      const hcpwa::LineSet<8>& bounds,
      const std::vector<hcpwa::LineSet<8>>& prisms0,
      const std::vector<hcpwa::LineSet<8>>& prisms1,
      const std::vector<hcpwa::LineSet<8>>& prisms2,
      const std::vector<hcpwa::LineSet<8>>& prisms3,
      const std::vector<hcpwa::LineSet<8>>& prisms4,
      std::vector<std::vector<size_t>>& out_indices,
      std::vector<std::vector<hcpwa::Vec<8>>>& out_points) {
    for (size_t idx0 = 0; idx0 < prisms0.size(); idx0++) {
      for (size_t idx1 = 0; idx1 < prisms1.size(); idx1++) {
        for (size_t idx2 = 0; idx2 < prisms2.size(); idx2++) {
          for (size_t idx3 = 0; idx3 < prisms3.size(); idx3++) {
            for (size_t idx4 = 0; idx4 < prisms4.size(); idx4++) {
              hcpwa::LineSet<8> concatenated_prisms;
              // IMPORTANT: bound the 8D polytope so it has vertices
              concatenated_prisms.append_range(bounds);
              concatenated_prisms.append_range(prisms0[idx0]);
              concatenated_prisms.append_range(prisms1[idx1]);
              concatenated_prisms.append_range(prisms2[idx2]);
              concatenated_prisms.append_range(prisms3[idx3]);
              concatenated_prisms.append_range(prisms4[idx4]);
              auto intersection = hcpwa::LinesToPoints(concatenated_prisms);
              if (intersection.size() > 0) {
                out_points.push_back(intersection);
                // Store the indices of the prisms that form the intersection
                std::vector<size_t> prism_indices = {idx0, idx1, idx2, idx3, idx4};
                out_indices.push_back(prism_indices);
              }
            }
          }
        }
      }
    }
  };

  // Intersect all prisms phase 0
  std::vector<std::vector<size_t>> intersection_prism_indices_phase0;
  std::vector<std::vector<hcpwa::Vec<8>>> intersection_points_phase0;
  compute_intersections(aabb_bounds, prisms31, prisms36, prisms24, prisms27, prisms58,
                        intersection_prism_indices_phase0, intersection_points_phase0);

  // Intersect all prisms phase 1
  std::vector<std::vector<size_t>> intersection_prism_indices_phase1;
  std::vector<std::vector<hcpwa::Vec<8>>> intersection_points_phase1;
  compute_intersections(aabb_bounds, prisms51, prisms57, prisms84, prisms86, prisms23,
                        intersection_prism_indices_phase1, intersection_points_phase1);
  
  AreasVerticesResult result;
  // Phase 0
  result.triangles31 = triangles31;
  result.triangles36 = triangles36;
  result.triangles24 = triangles24;
  result.triangles27 = triangles27;
  result.triangles58 = triangles58;
  // Phase 1
  result.triangles51 = triangles51;
  result.triangles57 = triangles57;
  result.triangles84 = triangles84;
  result.triangles86 = triangles86;
  result.triangles23 = triangles23;
  // Intersection points and indices
  result.intersection_points_phase0 = intersection_points_phase0;
  result.intersection_prism_indices_phase0 = intersection_prism_indices_phase0;
  result.intersection_points_phase1 = intersection_points_phase1;
  result.intersection_prism_indices_phase1 = intersection_prism_indices_phase1;

  // DEBUG: print result
  // Print triangle sizes for phase 0
  std::cout << "result.triangles31.size(): " << result.triangles31.size() << '\n';
  std::cout << "result.triangles36.size(): " << result.triangles36.size() << '\n';
  std::cout << "result.triangles24.size(): " << result.triangles24.size() << '\n';
  std::cout << "result.triangles27.size(): " << result.triangles27.size() << '\n';
  std::cout << "result.triangles58.size(): " << result.triangles58.size() << '\n';
  // Print triangle sizes for phase 1
  std::cout << "result.triangles51.size(): " << result.triangles51.size() << '\n';
  std::cout << "result.triangles57.size(): " << result.triangles57.size() << '\n';
  std::cout << "result.triangles84.size(): " << result.triangles84.size() << '\n';
  std::cout << "result.triangles86.size(): " << result.triangles86.size() << '\n';
  std::cout << "result.triangles23.size(): " << result.triangles23.size() << '\n';
  // Print intersection points sizes for phase 0
  std::cout << "result.intersection_points_phase0.size(): " << result.intersection_points_phase0.size() << '\n';
  std::cout << "result.intersection_prism_indices_phase0.size(): " << result.intersection_prism_indices_phase0.size() << '\n';
  // Print intersection points sizes for phase 1
  std::cout << "result.intersection_points_phase1.size(): " << result.intersection_points_phase1.size() << '\n';
  std::cout << "result.intersection_prism_indices_phase1.size(): " << result.intersection_prism_indices_phase1.size() << '\n';
  
  return result;
}

}  // namespace hcpwa

