#include <algo.hpp>
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <format>
#include <mutex>
#include <ostream>
#include <symbolic.hpp>
#include <morph.hpp>
#include <thread>
#include <vector>
#include "cddwrap/cdd.hpp"
#include "types.hpp"
#include "utility.hpp"

namespace hcpwa {

// NOLINTNEXTLINE
using namespace hcpwa::symbols;

PolygonResolutions compute_polygon_resolutions(
    double N, double F, double v, double w, double b51, double b57, double b84,
    double b86, double b31, double b36, double b24, double b27, double f2min,
    double f3min, double f5min, double f8min, double f2max, double f3max,
    double f5max, double f8max, bool verbose) {
  if (verbose) {
    std::cout << "N = " << N << std::endl;
    std::cout << "F = " << F << std::endl;
    std::cout << "v = " << v << std::endl;
    std::cout << "w = " << w << std::endl;
    std::cout << "b51 = " << b51 << std::endl;
    std::cout << "b57 = " << b57 << std::endl;
    std::cout << "b84 = " << b84 << std::endl;
    std::cout << "b86 = " << b86 << std::endl;
    std::cout << "b31 = " << b31 << std::endl;
    std::cout << "b36 = " << b36 << std::endl;
    std::cout << "b24 = " << b24 << std::endl;
    std::cout << "b27 = " << b27 << std::endl;
    std::cout << "f2min = " << f2min << std::endl;
    std::cout << "f3min = " << f3min << std::endl;
    std::cout << "f5min = " << f5min << std::endl;
    std::cout << "f8min = " << f8min << std::endl;
    std::cout << "f2max = " << f2max << std::endl;
    std::cout << "f3max = " << f3max << std::endl;
    std::cout << "f5max = " << f5max << std::endl;
    std::cout << "f8max = " << f8max << std::endl;
  }

  hcpwa::AABB<2> aabb2d
      = {{0, 0}, {static_cast<hcpwa::Float>(N), static_cast<hcpwa::Float>(N)}};

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
  hcpwa::LineSet<8> lines31;
  hcpwa::LineSet<8> lines36;
  hcpwa::LineSet<8> lines24;
  hcpwa::LineSet<8> lines27;
  hcpwa::LineSet<8> lines58;
  // Phase 1
  hcpwa::LineSet<8> lines51;
  hcpwa::LineSet<8> lines57;
  hcpwa::LineSet<8> lines84;
  hcpwa::LineSet<8> lines86;
  hcpwa::LineSet<8> lines23;
  {
    auto resolutions = MinResolutions(f51);
    constexpr int kDim
        = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim, 8>(resolutions);
    lines51.insert(lines51.end(), lines.begin(), lines.end());
  }
  {
    auto resolutions = MinResolutions(f57);
    constexpr int kDim
        = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim, 8>(resolutions);
    lines57.insert(lines57.end(), lines.begin(), lines.end());
  }
  {
    auto resolutions = MinResolutions(f84);
    constexpr int kDim
        = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim, 8>(resolutions);
    lines84.insert(lines84.end(), lines.begin(), lines.end());
  }
  {
    auto resolutions = MinResolutions(f86);
    constexpr int kDim
        = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim, 8>(resolutions);
    lines86.insert(lines86.end(), lines.begin(), lines.end());
  }
  {
    auto resolutions = MinResolutions(f31);
    constexpr int kDim
        = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim, 8>(resolutions);
    lines31.insert(lines31.end(), lines.begin(), lines.end());
  }
  {
    auto resolutions = MinResolutions(f36);
    constexpr int kDim
        = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim, 8>(resolutions);
    lines36.insert(lines36.end(), lines.begin(), lines.end());
  }
  {
    auto resolutions = MinResolutions(f24);
    constexpr int kDim
        = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim, 8>(resolutions);
    lines24.insert(lines24.end(), lines.begin(), lines.end());
  }
  {
    auto resolutions = MinResolutions(f27);
    constexpr int kDim
        = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim, 8>(resolutions);
    lines27.insert(lines27.end(), lines.begin(), lines.end());
  }
  // In planes lower bounds
  {
    auto resolutions = MinResolutions(f2in_lower);
    constexpr int kDim
        = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim, 8>(resolutions);
    lines24.insert(lines24.end(), lines.begin(), lines.end());
    lines27.insert(lines27.end(), lines.begin(), lines.end());
    lines23.insert(lines23.end(), lines.begin(), lines.end());
  }
  {
    auto resolutions = MinResolutions(f3in_lower);
    constexpr int kDim
        = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim, 8>(resolutions);
    lines31.insert(lines31.end(), lines.begin(), lines.end());
    lines36.insert(lines36.end(), lines.begin(), lines.end());
    lines23.insert(lines23.end(), lines.begin(), lines.end());
  }
  {
    auto resolutions = MinResolutions(f5in_lower);
    constexpr int kDim
        = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim, 8>(resolutions);
    lines51.insert(lines51.end(), lines.begin(), lines.end());
    lines57.insert(lines57.end(), lines.begin(), lines.end());
    lines58.insert(lines58.end(), lines.begin(), lines.end());
  }
  {
    auto resolutions = MinResolutions(f8in_lower);
    constexpr int kDim
        = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim, 8>(resolutions);
    lines84.insert(lines84.end(), lines.begin(), lines.end());
    lines86.insert(lines86.end(), lines.begin(), lines.end());
    lines58.insert(lines58.end(), lines.begin(), lines.end());
  }
  // In planes upper bounds
  {
    auto resolutions = MinResolutions(f2in_upper);
    constexpr int kDim
        = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim, 8>(resolutions);
    lines24.insert(lines24.end(), lines.begin(), lines.end());
    lines27.insert(lines27.end(), lines.begin(), lines.end());
    lines23.insert(lines23.end(), lines.begin(), lines.end());
  }
  {
    auto resolutions = MinResolutions(f3in_upper);
    constexpr int kDim
        = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim, 8>(resolutions);
    lines31.insert(lines31.end(), lines.begin(), lines.end());
    lines36.insert(lines36.end(), lines.begin(), lines.end());
    lines23.insert(lines23.end(), lines.begin(), lines.end());
  }
  {
    auto resolutions = MinResolutions(f5in_upper);
    constexpr int kDim
        = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim, 8>(resolutions);
    lines51.insert(lines51.end(), lines.begin(), lines.end());
    lines57.insert(lines57.end(), lines.begin(), lines.end());
    lines58.insert(lines58.end(), lines.begin(), lines.end());
  }
  {
    auto resolutions = MinResolutions(f8in_upper);
    constexpr int kDim
        = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim, 8>(resolutions);
    lines84.insert(lines84.end(), lines.begin(), lines.end());
    lines86.insert(lines86.end(), lines.begin(), lines.end());
    lines58.insert(lines58.end(), lines.begin(), lines.end());
  }
  // Out planes
  {
    auto resolutions = MinResolutions(f1out_lower);
    constexpr int kDim
        = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim, 8>(resolutions);
    lines51.insert(lines51.end(), lines.begin(), lines.end());
    lines31.insert(lines31.end(), lines.begin(), lines.end());
  }
  {
    auto resolutions = MinResolutions(f4out_lower);
    constexpr int kDim
        = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim, 8>(resolutions);
    lines84.insert(lines84.end(), lines.begin(), lines.end());
    lines24.insert(lines24.end(), lines.begin(), lines.end());
  }
  {
    auto resolutions = MinResolutions(f6out_lower);
    constexpr int kDim
        = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim, 8>(resolutions);
    lines36.insert(lines36.end(), lines.begin(), lines.end());
    lines86.insert(lines86.end(), lines.begin(), lines.end());
  }
  {
    auto resolutions = MinResolutions(f7out_lower);
    constexpr int kDim
        = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
    auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim, 8>(resolutions);
    lines27.insert(lines27.end(), lines.begin(), lines.end());
    lines57.insert(lines57.end(), lines.begin(), lines.end());
  }

  if (verbose) {
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
  }

  PolygonResolutions result;
  result.resolution_51 = hcpwa::SplitAABBWithLines(
      aabb2d, hcpwa::DimensionCast<2, 8>(lines51, {1 - 1, 5 - 1}));
  result.resolution_57 = hcpwa::SplitAABBWithLines(
      aabb2d, hcpwa::DimensionCast<2, 8>(lines57, {5 - 1, 7 - 1}));
  result.resolution_84 = hcpwa::SplitAABBWithLines(
      aabb2d, hcpwa::DimensionCast<2, 8>(lines84, {4 - 1, 8 - 1}));
  result.resolution_86 = hcpwa::SplitAABBWithLines(
      aabb2d, hcpwa::DimensionCast<2, 8>(lines86, {6 - 1, 8 - 1}));
  result.resolution_58 = hcpwa::SplitAABBWithLines(
      aabb2d, hcpwa::DimensionCast<2, 8>(lines58, {5 - 1, 8 - 1}));
  result.resolution_31 = hcpwa::SplitAABBWithLines(
      aabb2d, hcpwa::DimensionCast<2, 8>(lines31, {1 - 1, 3 - 1}));
  result.resolution_36 = hcpwa::SplitAABBWithLines(
      aabb2d, hcpwa::DimensionCast<2, 8>(lines36, {3 - 1, 6 - 1}));
  result.resolution_24 = hcpwa::SplitAABBWithLines(
      aabb2d, hcpwa::DimensionCast<2, 8>(lines24, {2 - 1, 4 - 1}));
  result.resolution_27 = hcpwa::SplitAABBWithLines(
      aabb2d, hcpwa::DimensionCast<2, 8>(lines27, {2 - 1, 7 - 1}));
  result.resolution_23 = hcpwa::SplitAABBWithLines(
      aabb2d, hcpwa::DimensionCast<2, 8>(lines23, {2 - 1, 3 - 1}));

  if (verbose) {
    std::cout << "Polygon resolution 51 count: " << result.resolution_51.size()
              << "\n";
    std::cout << "Polygon resolution 57 count: " << result.resolution_57.size()
              << "\n";
    std::cout << "Polygon resolution 84 count: " << result.resolution_84.size()
              << "\n";
    std::cout << "Polygon resolution 86 count: " << result.resolution_86.size()
              << "\n";
    std::cout << "Polygon resolution 58 count: " << result.resolution_58.size()
              << "\n";
    std::cout << "Polygon resolution 31 count: " << result.resolution_31.size()
              << "\n";
    std::cout << "Polygon resolution 36 count: " << result.resolution_36.size()
              << "\n";
    std::cout << "Polygon resolution 24 count: " << result.resolution_24.size()
              << "\n";
    std::cout << "Polygon resolution 27 count: " << result.resolution_27.size()
              << "\n";
    std::cout << "Polygon resolution 23 count: " << result.resolution_23.size()
              << "\n";
    for (int i = 0; i < result.resolution_51.size(); i++) {
      std::cout << "polygon " << i << ": ";
      for (int j = 0; j < result.resolution_51[i].polygon.size(); j++) {
        std::cout << result.resolution_51[i].polygon[j] << " ";
      }
      std::cout << std::endl;
    }
    for (int i = 0; i < result.resolution_57.size(); i++) {
      std::cout << "polygon " << i << ": ";
      for (int j = 0; j < result.resolution_57[i].polygon.size(); j++) {
        std::cout << result.resolution_57[i].polygon[j] << " ";
      }
      std::cout << "\n";
    }
    for (int i = 0; i < result.resolution_84.size(); i++) {
      std::cout << "polygon " << i << ": ";
      for (int j = 0; j < result.resolution_84[i].polygon.size(); j++) {
        std::cout << result.resolution_84[i].polygon[j] << " ";
      }
      std::cout << "\n";
    }
    for (int i = 0; i < result.resolution_86.size(); i++) {
      std::cout << "polygon " << i << ": ";
      for (int j = 0; j < result.resolution_86[i].polygon.size(); j++) {
        std::cout << result.resolution_86[i].polygon[j] << " ";
      }
      std::cout << "\n";
    }
    for (int i = 0; i < result.resolution_58.size(); i++) {
      std::cout << "polygon " << i << ": ";
      for (int j = 0; j < result.resolution_58[i].polygon.size(); j++) {
        std::cout << result.resolution_58[i].polygon[j] << " ";
      }
      std::cout << "\n";
    }
    for (int i = 0; i < result.resolution_31.size(); i++) {
      std::cout << "polygon " << i << ": ";
      for (int j = 0; j < result.resolution_31[i].polygon.size(); j++) {
        std::cout << result.resolution_31[i].polygon[j] << " ";
      }
      std::cout << "\n";
    }
    for (int i = 0; i < result.resolution_36.size(); i++) {
      std::cout << "polygon " << i << ": ";
      for (int j = 0; j < result.resolution_36[i].polygon.size(); j++) {
        std::cout << result.resolution_36[i].polygon[j] << " ";
      }
      std::cout << "\n";
    }
    for (int i = 0; i < result.resolution_24.size(); i++) {
      std::cout << "polygon " << i << ": ";
      for (int j = 0; j < result.resolution_24[i].polygon.size(); j++) {
        std::cout << result.resolution_24[i].polygon[j] << " ";
      }
      std::cout << "\n";
    }
    for (int i = 0; i < result.resolution_27.size(); i++) {
      std::cout << "polygon " << i << ": ";
      for (int j = 0; j < result.resolution_27[i].polygon.size(); j++) {
        std::cout << result.resolution_27[i].polygon[j] << " ";
      }
      std::cout << "\n";
    }
    for (int i = 0; i < result.resolution_23.size(); i++) {
      std::cout << "polygon " << i << ": ";
      for (int j = 0; j < result.resolution_23[i].polygon.size(); j++) {
        std::cout << result.resolution_23[i].polygon[j] << " ";
      }
      std::cout << std::endl;
    }
  }
  return result;
}

TriangulationAndPrismsResult compute_triangulation_and_prisms(
    PolygonResolutions& polygon_resolutions, const hcpwa::AABB<2>& aabb2d,
    bool verbose) {
  auto& polygon_resolution_51 = polygon_resolutions.resolution_51;
  auto& polygon_resolution_57 = polygon_resolutions.resolution_57;
  auto& polygon_resolution_84 = polygon_resolutions.resolution_84;
  auto& polygon_resolution_86 = polygon_resolutions.resolution_86;
  auto& polygon_resolution_58 = polygon_resolutions.resolution_58;
  auto& polygon_resolution_31 = polygon_resolutions.resolution_31;
  auto& polygon_resolution_36 = polygon_resolutions.resolution_36;
  auto& polygon_resolution_24 = polygon_resolutions.resolution_24;
  auto& polygon_resolution_27 = polygon_resolutions.resolution_27;
  auto& polygon_resolution_23 = polygon_resolutions.resolution_23;

  if (verbose) {
    std::cout << "Triangulated polygons with unique vertices: "
              << "\n";
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
  // Triangulated polygons with unique vertices
  auto polygons51
      = hcpwa::GetTrianglesWithUniqueVertices(aabb2d, polygon_resolution_51);
  auto polygons57
      = hcpwa::GetTrianglesWithUniqueVertices(aabb2d, polygon_resolution_57);
  auto polygons84
      = hcpwa::GetTrianglesWithUniqueVertices(aabb2d, polygon_resolution_84);
  auto polygons86
      = hcpwa::GetTrianglesWithUniqueVertices(aabb2d, polygon_resolution_86);
  auto polygons58
      = hcpwa::GetTrianglesWithUniqueVertices(aabb2d, polygon_resolution_58);
  auto polygons31
      = hcpwa::GetTrianglesWithUniqueVertices(aabb2d, polygon_resolution_31);
  auto polygons36
      = hcpwa::GetTrianglesWithUniqueVertices(aabb2d, polygon_resolution_36);
  auto polygons24
      = hcpwa::GetTrianglesWithUniqueVertices(aabb2d, polygon_resolution_24);
  auto polygons27
      = hcpwa::GetTrianglesWithUniqueVertices(aabb2d, polygon_resolution_27);
  auto polygons23
      = hcpwa::GetTrianglesWithUniqueVertices(aabb2d, polygon_resolution_23);

  if (verbose) {
    std::cout << "Triangles count: " << polygons51.size() << " "
              << polygons57.size() << " " << polygons84.size() << " "
              << polygons86.size() << " " << polygons58.size() << " "
              << polygons31.size() << " " << polygons36.size() << " "
              << polygons24.size() << " " << polygons27.size() << " "
              << polygons23.size() << "\n";
    for (int i = 0; i < polygons51.size(); i++) {
      std::cout << "triangle " << i << ": " << polygons51[i].a << " "
                << polygons51[i].b << " " << polygons51[i].c << "\n";
    }
    for (int i = 0; i < polygons57.size(); i++) {
      std::cout << "triangle " << i << ": " << polygons57[i].a << " "
                << polygons57[i].b << " " << polygons57[i].c << "\n";
    }
    for (int i = 0; i < polygons84.size(); i++) {
      std::cout << "triangle " << i << ": " << polygons84[i].a << " "
                << polygons84[i].b << " " << polygons84[i].c << "\n";
    }
    for (int i = 0; i < polygons86.size(); i++) {
      std::cout << "triangle " << i << ": " << polygons86[i].a << " "
                << polygons86[i].b << " " << polygons86[i].c << "\n";
    }
    for (int i = 0; i < polygons58.size(); i++) {
      std::cout << "triangle " << i << ": " << polygons58[i].a << " "
                << polygons58[i].b << " " << polygons58[i].c << "\n";
    }
    for (int i = 0; i < polygons31.size(); i++) {
      std::cout << "triangle " << i << ": " << polygons31[i].a << " "
                << polygons31[i].b << " " << polygons31[i].c << "\n";
    }
    for (int i = 0; i < polygons36.size(); i++) {
      std::cout << "triangle " << i << ": " << polygons36[i].a << " "
                << polygons36[i].b << " " << polygons36[i].c << "\n";
    }
    for (int i = 0; i < polygons24.size(); i++) {
      std::cout << "triangle " << i << ": " << polygons24[i].a << " "
                << polygons24[i].b << " " << polygons24[i].c << "\n";
    }
    for (int i = 0; i < polygons27.size(); i++) {
      std::cout << "triangle " << i << ": " << polygons27[i].a << " "
                << polygons27[i].b << " " << polygons27[i].c << "\n";
    }
    for (int i = 0; i < polygons23.size(); i++) {
      std::cout << "triangle " << i << ": " << polygons23[i].a << " "
                << polygons23[i].b << " " << polygons23[i].c << std::endl;
    }
  }
  // Phase 0
  for (auto& triangle : polygons31) {
    auto prism = hcpwa::CalcPrism(triangle, {1 - 1, 3 - 1});
    prisms31.push_back(prism);
  }
  for (auto& triangle : polygons36) {
    auto prism = hcpwa::CalcPrism(triangle, {3 - 1, 6 - 1});
    prisms36.push_back(prism);
  }
  for (auto& triangle : polygons24) {
    auto prism = hcpwa::CalcPrism(triangle, {2 - 1, 4 - 1});
    prisms24.push_back(prism);
  }
  for (auto& triangle : polygons27) {
    auto prism = hcpwa::CalcPrism(triangle, {2 - 1, 7 - 1});
    prisms27.push_back(prism);
  }
  for (auto& triangle : polygons58) {
    auto prism = hcpwa::CalcPrism(triangle, {5 - 1, 8 - 1});
    prisms58.push_back(prism);
  }
  // Phase 1
  for (auto& triangle : polygons51) {
    auto prism = hcpwa::CalcPrism(triangle, {1 - 1, 5 - 1});
    prisms51.push_back(prism);
  }
  for (auto& triangle : polygons57) {
    auto prism = hcpwa::CalcPrism(triangle, {5 - 1, 7 - 1});
    prisms57.push_back(prism);
  }
  for (auto& triangle : polygons84) {
    auto prism = hcpwa::CalcPrism(triangle, {4 - 1, 8 - 1});
    prisms84.push_back(prism);
  }
  for (auto& triangle : polygons86) {
    auto prism = hcpwa::CalcPrism(triangle, {6 - 1, 8 - 1});
    prisms86.push_back(prism);
  }
  for (auto& triangle : polygons23) {
    auto prism = hcpwa::CalcPrism(triangle, {2 - 1, 3 - 1});
    prisms23.push_back(prism);
  }

  TriangulationAndPrismsResult result;
  result.prisms31 = std::move(prisms31);
  result.prisms36 = std::move(prisms36);
  result.prisms24 = std::move(prisms24);
  result.prisms27 = std::move(prisms27);
  result.prisms58 = std::move(prisms58);
  result.prisms51 = std::move(prisms51);
  result.prisms57 = std::move(prisms57);
  result.prisms84 = std::move(prisms84);
  result.prisms86 = std::move(prisms86);
  result.prisms23 = std::move(prisms23);
  result.triangles31 = std::move(polygons31);
  result.triangles36 = std::move(polygons36);
  result.triangles24 = std::move(polygons24);
  result.triangles27 = std::move(polygons27);
  result.triangles58 = std::move(polygons58);
  result.triangles51 = std::move(polygons51);
  result.triangles57 = std::move(polygons57);
  result.triangles84 = std::move(polygons84);
  result.triangles86 = std::move(polygons86);
  result.triangles23 = std::move(polygons23);
  return result;
}

PolygonPrismsResult compute_prisms_from_polygons(
    PolygonResolutions& polygon_resolutions, const hcpwa::AABB<2>& aabb2d,
    bool verbose) {
  auto& polygons51 = polygon_resolutions.resolution_51;
  auto& polygons57 = polygon_resolutions.resolution_57;
  auto& polygons84 = polygon_resolutions.resolution_84;
  auto& polygons86 = polygon_resolutions.resolution_86;
  auto& polygons58 = polygon_resolutions.resolution_58;
  auto& polygons31 = polygon_resolutions.resolution_31;
  auto& polygons36 = polygon_resolutions.resolution_36;
  auto& polygons24 = polygon_resolutions.resolution_24;
  auto& polygons27 = polygon_resolutions.resolution_27;
  auto& polygons23 = polygon_resolutions.resolution_23;

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

  if (verbose) {
    std::cout << "Triangles count: " << polygons51.size() << " "
              << polygons57.size() << " " << polygons84.size() << " "
              << polygons86.size() << " " << polygons58.size() << " "
              << polygons31.size() << " " << polygons36.size() << " "
              << polygons24.size() << " " << polygons27.size() << " "
              << polygons23.size() << "\n";
  }
  // Phase 0
  for (auto& poly_res : polygons31) {
    auto prism = hcpwa::CalcPrism(poly_res.polygon, {1 - 1, 3 - 1});
    prisms31.push_back(prism);
  }
  for (auto& poly_res : polygons36) {
    auto prism = hcpwa::CalcPrism(poly_res.polygon, {3 - 1, 6 - 1});
    prisms36.push_back(prism);
  }
  for (auto& poly_res : polygons24) {
    auto prism = hcpwa::CalcPrism(poly_res.polygon, {2 - 1, 4 - 1});
    prisms24.push_back(prism);
  }
  for (auto& poly_res : polygons27) {
    auto prism = hcpwa::CalcPrism(poly_res.polygon, {2 - 1, 7 - 1});
    prisms27.push_back(prism);
  }
  for (auto& poly_res : polygons58) {
    auto prism = hcpwa::CalcPrism(poly_res.polygon, {5 - 1, 8 - 1});
    prisms58.push_back(prism);
  }
  // Phase 1
  for (auto& poly_res : polygons51) {
    auto prism = hcpwa::CalcPrism(poly_res.polygon, {1 - 1, 5 - 1});
    prisms51.push_back(prism);
  }
  for (auto& poly_res : polygons57) {
    auto prism = hcpwa::CalcPrism(poly_res.polygon, {5 - 1, 7 - 1});
    prisms57.push_back(prism);
  }
  for (auto& poly_res : polygons84) {
    auto prism = hcpwa::CalcPrism(poly_res.polygon, {4 - 1, 8 - 1});
    prisms84.push_back(prism);
  }
  for (auto& poly_res : polygons86) {
    auto prism = hcpwa::CalcPrism(poly_res.polygon, {6 - 1, 8 - 1});
    prisms86.push_back(prism);
  }
  for (auto& poly_res : polygons23) {
    auto prism = hcpwa::CalcPrism(poly_res.polygon, {2 - 1, 3 - 1});
    prisms23.push_back(prism);
  }

  PolygonPrismsResult result;
  result.prisms31 = std::move(prisms31);
  result.prisms36 = std::move(prisms36);
  result.prisms24 = std::move(prisms24);
  result.prisms27 = std::move(prisms27);
  result.prisms58 = std::move(prisms58);
  result.prisms51 = std::move(prisms51);
  result.prisms57 = std::move(prisms57);
  result.prisms84 = std::move(prisms84);
  result.prisms86 = std::move(prisms86);
  result.prisms23 = std::move(prisms23);
  return result;
}

PhaseIntersectionResult compute_intersection_points(
    const std::vector<hcpwa::LineSet<8>>& prisms31,
    const std::vector<hcpwa::LineSet<8>>& prisms36,
    const std::vector<hcpwa::LineSet<8>>& prisms24,
    const std::vector<hcpwa::LineSet<8>>& prisms27,
    const std::vector<hcpwa::LineSet<8>>& prisms58,
    const std::vector<hcpwa::LineSet<8>>& prisms51,
    const std::vector<hcpwa::LineSet<8>>& prisms57,
    const std::vector<hcpwa::LineSet<8>>& prisms84,
    const std::vector<hcpwa::LineSet<8>>& prisms86,
    const std::vector<hcpwa::LineSet<8>>& prisms23,
    const std::vector<hcpwa::TriangleWithUniqueVertices>& triangles58,
    const std::vector<hcpwa::TriangleWithUniqueVertices>& triangles23,
    hcpwa::Float N, bool verbose) {
  hcpwa::AABB<3> aabb3d = {{0, 0, 0}, {N, N, N}};
  const auto aabb3d_bounds = hcpwa::AABBBounds(aabb3d);

  auto computend = []<int Dim>(
                       std::array<int, Dim> dims,
                       const hcpwa::LineSet<Dim>& bounds,
                       const std::vector<hcpwa::LineSet<8>>& prisms0,
                       const std::vector<hcpwa::LineSet<8>>& prisms1,
                       std::vector<std::vector<size_t>>& out_indices,
                       std::vector<std::vector<hcpwa::Vec<Dim>>>& out_points) {
    for (size_t idx0 = 0; idx0 < prisms0.size(); idx0++) {
      for (size_t idx1 = 0; idx1 < prisms1.size(); idx1++) {
        hcpwa::LineSet<Dim> concatenated_prisms = bounds;

        for (auto& i : prisms0[idx0]) {
          concatenated_prisms.push_back(hcpwa::DimensionCast<Dim, 8>(i, dims));
        }
        for (auto& i : prisms1[idx1]) {
          concatenated_prisms.push_back(hcpwa::DimensionCast<Dim, 8>(i, dims));
        }
        auto intersection = hcpwa::LinesToPoints<Dim>(concatenated_prisms);
        if (intersection.size() > 0) {
          out_points.push_back(intersection);
          // Store the indices of the prisms that form the intersection
          std::vector<size_t> prism_indices = {idx0, idx1};
          out_indices.push_back(prism_indices);
        }
      }
    }
  };

  if (verbose) {
    std::cout << std::format("Prism 31: {}", prisms31.size()) << std::endl;
    std::cout << std::format("Prism 36: {}", prisms36.size()) << std::endl;
    std::cout << std::format("Prism 24: {}", prisms24.size()) << std::endl;
    std::cout << std::format("Prism 27: {}", prisms27.size()) << std::endl;
  }

  std::vector<std::vector<size_t>> intersection_prism_indices_136;
  std::vector<std::vector<hcpwa::Vec<3>>> intersection_points_136;
  computend({0, 2, 5}, aabb3d_bounds, prisms31, prisms36,
            intersection_prism_indices_136, intersection_points_136);

  std::vector<std::vector<size_t>> intersection_prism_indices_247;
  std::vector<std::vector<hcpwa::Vec<3>>> intersection_points_247;
  computend({1, 3, 6}, aabb3d_bounds, prisms24, prisms27,
            intersection_prism_indices_247, intersection_points_247);

  std::vector<std::vector<size_t>> intersection_prism_indices_157;
  std::vector<std::vector<hcpwa::Vec<3>>> intersection_points_157;
  computend({0, 4, 6}, aabb3d_bounds, prisms51, prisms57,
            intersection_prism_indices_157, intersection_points_157);

  std::vector<std::vector<size_t>> intersection_prism_indices_468;
  std::vector<std::vector<hcpwa::Vec<3>>> intersection_points_468;
  computend({3, 5, 7}, aabb3d_bounds, prisms84, prisms86,
            intersection_prism_indices_468, intersection_points_468);

  if (verbose) {
    std::cout << "Intersection counted:" << std::endl;
    std::cout << "\t136 count: " << intersection_points_136.size() << std::endl;
    std::cout << "\t247 count: " << intersection_points_247.size() << std::endl;
    std::cout << "\t58 count: " << triangles58.size() << std::endl;
    std::cout << "\t157 count: " << intersection_points_157.size() << std::endl;
    std::cout << "\t468 count: " << intersection_points_468.size() << std::endl;
    std::cout << "\t23 count: " << triangles23.size() << std::endl;
  }

  std::vector<std::vector<size_t>> intersection_prism_indices_phase0;
  std::vector<std::vector<hcpwa::Vec<8>>> intersection_points_phase0;
  std::vector<std::vector<size_t>> intersection_prism_indices_phase1;
  std::vector<std::vector<hcpwa::Vec<8>>> intersection_points_phase1;

  for (size_t i136 = 0; i136 < intersection_points_136.size(); i136++) {
    for (size_t i247 = 0; i247 < intersection_prism_indices_247.size();
         i247++) {
      // auto t_start = std::chrono::high_resolution_clock::now();
      for (size_t i58 = 0; i58 < triangles58.size(); i58++) {
        intersection_points_phase0.emplace_back();
        intersection_prism_indices_phase0.emplace_back();
        const auto& indices_136 = intersection_prism_indices_136[i136];
        const auto& indices_247 = intersection_prism_indices_247[i247];
        intersection_prism_indices_phase0.back().insert(
            intersection_prism_indices_phase0.back().end(), indices_136.begin(),
            indices_136.end());
        intersection_prism_indices_phase0.back().insert(
            intersection_prism_indices_phase0.back().end(), indices_247.begin(),
            indices_247.end());
        intersection_prism_indices_phase0.back().push_back(i58);

        for (size_t j136 = 0;
             j136 < intersection_prism_indices_136[i136].size(); j136++) {
          for (size_t j247 = 0;
               j247 < intersection_prism_indices_247[i247].size(); j247++) {
            for (size_t j58 = 0; j58 < triangles58[i58].size(); j58++) {
              const auto& v136 = intersection_points_136[i136][j136];
              const auto& v247 = intersection_points_247[i247][j247];
              const auto& v58 = triangles58[i58][j58];
              hcpwa::Vec<8> v = kZeroVec;
              v[0] = v136[0];
              v[1] = v247[0];
              v[2] = v136[1];
              v[3] = v247[1];
              v[4] = v58[0];
              v[5] = v136[2];
              v[6] = v247[2];
              v[7] = v58[1];
              intersection_points_phase0.back().push_back(v);
            }
          }
        }
      }
      // auto t_end = std::chrono::high_resolution_clock::now();
      // std::chrono::duration<double> t_diff = t_end - t_start;
      // std::cout << "[Timing] i136=" << i136 << ", i247=" << i247 << ": " <<
      // t_diff.count() << "s" << std::endl;
    }
  }

  for (size_t i157 = 0; i157 < intersection_points_157.size(); i157++) {
    for (size_t i468 = 0; i468 < intersection_prism_indices_468.size();
         i468++) {
      for (size_t i23 = 0; i23 < triangles23.size(); i23++) {
        intersection_points_phase1.emplace_back();
        intersection_prism_indices_phase1.emplace_back();
        const auto& indices_157 = intersection_prism_indices_157[i157];
        const auto& indices_468 = intersection_prism_indices_468[i468];
        intersection_prism_indices_phase1.back().insert(
            intersection_prism_indices_phase1.back().end(), indices_157.begin(),
            indices_157.end());
        intersection_prism_indices_phase1.back().insert(
            intersection_prism_indices_phase1.back().end(), indices_468.begin(),
            indices_468.end());
        intersection_prism_indices_phase1.back().push_back(i23);

        for (size_t j157 = 0;
             j157 < intersection_prism_indices_157[i157].size(); j157++) {
          for (size_t j468 = 0;
               j468 < intersection_prism_indices_468[i468].size(); j468++) {
            for (size_t j23 = 0; j23 < triangles23[i23].size(); j23++) {
              const auto& v157 = intersection_points_157[i157][j157];
              const auto& v468 = intersection_points_468[i468][j468];
              const auto& v23 = triangles23[i23][j23];
              hcpwa::Vec<8> v = kZeroVec;
              v[0] = v157[0];
              v[1] = v23[0];
              v[2] = v23[1];
              v[3] = v468[0];
              v[4] = v157[1];
              v[5] = v468[1];
              v[6] = v157[2];
              v[7] = v468[2];
              intersection_points_phase1.back().push_back(v);
            }
          }
        }
      }
    }
  }

  PhaseIntersectionResult result;
  result.intersection_prism_indices_phase0
      = std::move(intersection_prism_indices_phase0);
  result.intersection_points_phase0 = std::move(intersection_points_phase0);
  result.intersection_prism_indices_phase1
      = std::move(intersection_prism_indices_phase1);
  result.intersection_points_phase1 = std::move(intersection_points_phase1);
  return result;
}

PhaseIntersectionResult compute_intersection_points(
  const std::vector<hcpwa::LineSet<8>>& prisms31,
  const std::vector<hcpwa::LineSet<8>>& prisms36,
  const std::vector<hcpwa::LineSet<8>>& prisms24,
  const std::vector<hcpwa::LineSet<8>>& prisms27,
  const std::vector<hcpwa::LineSet<8>>& prisms58,
  const std::vector<hcpwa::LineSet<8>>& prisms51,
  const std::vector<hcpwa::LineSet<8>>& prisms57,
  const std::vector<hcpwa::LineSet<8>>& prisms84,
  const std::vector<hcpwa::LineSet<8>>& prisms86,
  const std::vector<hcpwa::LineSet<8>>& prisms23,
  const std::vector<hcpwa::PolygonResolution>& polygons58,
  const std::vector<hcpwa::PolygonResolution>& polygons23,
  hcpwa::Float N, bool verbose) {
hcpwa::AABB<3> aabb3d = {{0, 0, 0}, {N, N, N}};
const auto aabb3d_bounds = hcpwa::AABBBounds(aabb3d);

auto computend = []<int Dim>(
                     std::array<int, Dim> dims,
                     const hcpwa::LineSet<Dim>& bounds,
                     const std::vector<hcpwa::LineSet<8>>& prisms0,
                     const std::vector<hcpwa::LineSet<8>>& prisms1,
                     std::vector<std::vector<size_t>>& out_indices,
                     std::vector<std::vector<hcpwa::Vec<Dim>>>& out_points) {
  for (size_t idx0 = 0; idx0 < prisms0.size(); idx0++) {
    for (size_t idx1 = 0; idx1 < prisms1.size(); idx1++) {
      hcpwa::LineSet<Dim> concatenated_prisms = bounds;

      for (auto& i : prisms0[idx0]) {
        concatenated_prisms.push_back(hcpwa::DimensionCast<Dim, 8>(i, dims));
      }
      for (auto& i : prisms1[idx1]) {
        concatenated_prisms.push_back(hcpwa::DimensionCast<Dim, 8>(i, dims));
      }
      auto intersection = hcpwa::LinesToPoints<Dim>(concatenated_prisms);
      if (intersection.size() > 0) {
        out_points.push_back(intersection);
        // Store the indices of the prisms that form the intersection
        std::vector<size_t> prism_indices = {idx0, idx1};
        out_indices.push_back(prism_indices);
      }
    }
  }
};

if (verbose) {
  std::cout << std::format("Prism 31: {}", prisms31.size()) << std::endl;
  std::cout << std::format("Prism 36: {}", prisms36.size()) << std::endl;
  std::cout << std::format("Prism 24: {}", prisms24.size()) << std::endl;
  std::cout << std::format("Prism 27: {}", prisms27.size()) << std::endl;
}

std::vector<std::vector<size_t>> intersection_prism_indices_136;
std::vector<std::vector<hcpwa::Vec<3>>> intersection_points_136;
computend({0, 2, 5}, aabb3d_bounds, prisms31, prisms36,
          intersection_prism_indices_136, intersection_points_136);

std::vector<std::vector<size_t>> intersection_prism_indices_247;
std::vector<std::vector<hcpwa::Vec<3>>> intersection_points_247;
computend({1, 3, 6}, aabb3d_bounds, prisms24, prisms27,
          intersection_prism_indices_247, intersection_points_247);

std::vector<std::vector<size_t>> intersection_prism_indices_157;
std::vector<std::vector<hcpwa::Vec<3>>> intersection_points_157;
computend({0, 4, 6}, aabb3d_bounds, prisms51, prisms57,
          intersection_prism_indices_157, intersection_points_157);

std::vector<std::vector<size_t>> intersection_prism_indices_468;
std::vector<std::vector<hcpwa::Vec<3>>> intersection_points_468;
computend({3, 5, 7}, aabb3d_bounds, prisms84, prisms86,
          intersection_prism_indices_468, intersection_points_468);

if (verbose) {
  std::cout << "Intersection counted:" << std::endl;
  std::cout << "\t136 count: " << intersection_points_136.size() << std::endl;
  std::cout << "\t247 count: " << intersection_points_247.size() << std::endl;
  std::cout << "\t58 count: " << polygons58.size() << std::endl;
  std::cout << "\t157 count: " << intersection_points_157.size() << std::endl;
  std::cout << "\t468 count: " << intersection_points_468.size() << std::endl;
  std::cout << "\t23 count: " << polygons23.size() << std::endl;
}

std::vector<std::vector<size_t>> intersection_prism_indices_phase0;
std::vector<std::vector<hcpwa::Vec<8>>> intersection_points_phase0;
std::vector<std::vector<size_t>> intersection_prism_indices_phase1;
std::vector<std::vector<hcpwa::Vec<8>>> intersection_points_phase1;

for (size_t i136 = 0; i136 < intersection_points_136.size(); i136++) {
  for (size_t i247 = 0; i247 < intersection_prism_indices_247.size();
       i247++) {
    // auto t_start = std::chrono::high_resolution_clock::now();
    for (size_t i58 = 0; i58 < polygons58.size(); i58++) {
      intersection_points_phase0.emplace_back();
      intersection_prism_indices_phase0.emplace_back();
      const auto& indices_136 = intersection_prism_indices_136[i136];
      const auto& indices_247 = intersection_prism_indices_247[i247];
      intersection_prism_indices_phase0.back().insert(
          intersection_prism_indices_phase0.back().end(), indices_136.begin(),
          indices_136.end());
      intersection_prism_indices_phase0.back().insert(
          intersection_prism_indices_phase0.back().end(), indices_247.begin(),
          indices_247.end());
      intersection_prism_indices_phase0.back().push_back(i58);

      for (size_t j136 = 0;
           j136 < intersection_prism_indices_136[i136].size(); j136++) {
        for (size_t j247 = 0;
             j247 < intersection_prism_indices_247[i247].size(); j247++) {
          for (size_t j58 = 0; j58 < polygons58[i58].polygon.size(); j58++) {
            const auto& v136 = intersection_points_136[i136][j136];
            const auto& v247 = intersection_points_247[i247][j247];
            const auto& v58 = polygons58[i58].polygon[j58];
            hcpwa::Vec<8> v = kZeroVec;
            v[0] = v136[0];
            v[1] = v247[0];
            v[2] = v136[1];
            v[3] = v247[1];
            v[4] = v58[0];
            v[5] = v136[2];
            v[6] = v247[2];
            v[7] = v58[1];
            intersection_points_phase0.back().push_back(v);
          }
        }
      }
    }
    // auto t_end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> t_diff = t_end - t_start;
    // std::cout << "[Timing] i136=" << i136 << ", i247=" << i247 << ": " <<
    // t_diff.count() << "s" << std::endl;
  }
}

for (size_t i157 = 0; i157 < intersection_points_157.size(); i157++) {
  for (size_t i468 = 0; i468 < intersection_prism_indices_468.size();
       i468++) {
    for (size_t i23 = 0; i23 < polygons23.size(); i23++) {
      intersection_points_phase1.emplace_back();
      intersection_prism_indices_phase1.emplace_back();
      const auto& indices_157 = intersection_prism_indices_157[i157];
      const auto& indices_468 = intersection_prism_indices_468[i468];
      intersection_prism_indices_phase1.back().insert(
          intersection_prism_indices_phase1.back().end(), indices_157.begin(),
          indices_157.end());
      intersection_prism_indices_phase1.back().insert(
          intersection_prism_indices_phase1.back().end(), indices_468.begin(),
          indices_468.end());
      intersection_prism_indices_phase1.back().push_back(i23);

      for (size_t j157 = 0;
           j157 < intersection_prism_indices_157[i157].size(); j157++) {
        for (size_t j468 = 0;
             j468 < intersection_prism_indices_468[i468].size(); j468++) {
          for (size_t j23 = 0; j23 < polygons23[i23].polygon.size(); j23++) {
            const auto& v157 = intersection_points_157[i157][j157];
            const auto& v468 = intersection_points_468[i468][j468];
            const auto& v23 = polygons23[i23].polygon[j23];
            hcpwa::Vec<8> v = kZeroVec;
            v[0] = v157[0];
            v[1] = v23[0];
            v[2] = v23[1];
            v[3] = v468[0];
            v[4] = v157[1];
            v[5] = v468[1];
            v[6] = v157[2];
            v[7] = v468[2];
            intersection_points_phase1.back().push_back(v);
          }
        }
      }
    }
  }
}

PhaseIntersectionResult result;
result.intersection_prism_indices_phase0
    = std::move(intersection_prism_indices_phase0);
result.intersection_points_phase0 = std::move(intersection_points_phase0);
result.intersection_prism_indices_phase1
    = std::move(intersection_prism_indices_phase1);
result.intersection_points_phase1 = std::move(intersection_points_phase1);
return result;
}

TriangleAreasVerticesResult compute_triangle_areas_vertices(
    double N, double F, double v, double w, double b51, double b57, double b84,
    double b86, double b31, double b36, double b24, double b27, double f2min,
    double f3min, double f5min, double f8min, double f2max, double f3max,
    double f5max, double f8max, bool verbose) {
  // Compute polygon min resolutions (splits) for each hyperplane
  auto polygon_resolutions = compute_polygon_resolutions(
      N, F, v, w, b51, b57, b84, b86, b31, b36, b24, b27, f2min, f3min, f5min,
      f8min, f2max, f3max, f5max, f8max, verbose);

  hcpwa::AABB<8> aabb
      = {{0, 0, 0, 0, 0, 0, 0, 0},
         {static_cast<hcpwa::Float>(N), static_cast<hcpwa::Float>(N),
          static_cast<hcpwa::Float>(N), static_cast<hcpwa::Float>(N),
          static_cast<hcpwa::Float>(N), static_cast<hcpwa::Float>(N),
          static_cast<hcpwa::Float>(N), static_cast<hcpwa::Float>(N)}};
  hcpwa::AABB<2> aabb2d
      = {{0, 0}, {static_cast<hcpwa::Float>(N), static_cast<hcpwa::Float>(N)}};
  const auto aabb_bounds = hcpwa::AABBBounds(aabb);

  // Triangulate computed polygons and return computed triangles and 8D prisms
  auto triangulation_result
      = compute_triangulation_and_prisms(polygon_resolutions, aabb2d, verbose);

  // Intersects prisms and returns intersection points for poth phases
  auto& prisms31 = triangulation_result.prisms31;
  auto& prisms36 = triangulation_result.prisms36;
  auto& prisms24 = triangulation_result.prisms24;
  auto& prisms27 = triangulation_result.prisms27;
  auto& prisms58 = triangulation_result.prisms58;
  auto& prisms51 = triangulation_result.prisms51;
  auto& prisms57 = triangulation_result.prisms57;
  auto& prisms84 = triangulation_result.prisms84;
  auto& prisms86 = triangulation_result.prisms86;
  auto& prisms23 = triangulation_result.prisms23;
  auto& polygons31 = triangulation_result.triangles31;
  auto& polygons36 = triangulation_result.triangles36;
  auto& polygons24 = triangulation_result.triangles24;
  auto& polygons27 = triangulation_result.triangles27;
  auto& polygons58 = triangulation_result.triangles58;
  auto& polygons51 = triangulation_result.triangles51;
  auto& polygons57 = triangulation_result.triangles57;
  auto& polygons84 = triangulation_result.triangles84;
  auto& polygons86 = triangulation_result.triangles86;
  auto& polygons23 = triangulation_result.triangles23;

  // Intersects prisms and returns intersection points for both phases
  auto intersection_result = compute_intersection_points(
      prisms31, prisms36, prisms24, prisms27, prisms58, prisms51, prisms57,
      prisms84, prisms86, prisms23, polygons58, polygons23,
      static_cast<hcpwa::Float>(N), verbose);
  auto& intersection_points_phase0
      = intersection_result.intersection_points_phase0;
  auto& intersection_prism_indices_phase0
      = intersection_result.intersection_prism_indices_phase0;
  auto& intersection_points_phase1
      = intersection_result.intersection_points_phase1;
  auto& intersection_prism_indices_phase1
      = intersection_result.intersection_prism_indices_phase1;

  TriangleAreasVerticesResult result;
  // Phase 0
  result.triangles31 = polygons31;
  result.triangles36 = polygons36;
  result.triangles24 = polygons24;
  result.triangles27 = polygons27;
  result.triangles58 = polygons58;
  // Phase 1
  result.triangles51 = polygons51;
  result.triangles57 = polygons57;
  result.triangles84 = polygons84;
  result.triangles86 = polygons86;
  result.triangles23 = polygons23;
  // Intersection points and indices
  result.intersection_points_phase0 = intersection_points_phase0;
  result.intersection_prism_indices_phase0 = intersection_prism_indices_phase0;
  result.intersection_points_phase1 = intersection_points_phase1;
  result.intersection_prism_indices_phase1 = intersection_prism_indices_phase1;

  if (verbose) {
    std::cout << "result.triangles31.size(): " << result.triangles31.size()
              << '\n';
    std::cout << "result.triangles36.size(): " << result.triangles36.size()
              << '\n';
    std::cout << "result.triangles24.size(): " << result.triangles24.size()
              << '\n';
    std::cout << "result.triangles27.size(): " << result.triangles27.size()
              << '\n';
    std::cout << "result.triangles58.size(): " << result.triangles58.size()
              << '\n';
    std::cout << "result.triangles51.size(): " << result.triangles51.size()
              << '\n';
    std::cout << "result.triangles57.size(): " << result.triangles57.size()
              << '\n';
    std::cout << "result.triangles84.size(): " << result.triangles84.size()
              << '\n';
    std::cout << "result.triangles86.size(): " << result.triangles86.size()
              << '\n';
    std::cout << "result.triangles23.size(): " << result.triangles23.size()
              << '\n';
    std::cout << "result.intersection_points_phase0.size(): "
              << result.intersection_points_phase0.size() << '\n';
    std::cout << "result.intersection_prism_indices_phase0.size(): "
              << result.intersection_prism_indices_phase0.size() << '\n';
    std::cout << "result.intersection_points_phase1.size(): "
              << result.intersection_points_phase1.size() << '\n';
    std::cout << "result.intersection_prism_indices_phase1.size(): "
              << result.intersection_prism_indices_phase1.size() << '\n';
  }

  return result;
}

PolygonAreasVerticesResult compute_polygon_areas_vertices(
    double N, double F, double v, double w, double b51, double b57, double b84,
    double b86, double b31, double b36, double b24, double b27, double f2min,
    double f3min, double f5min, double f8min, double f2max, double f3max,
    double f5max, double f8max, bool verbose) {
  // Compute polygon min resolutions (splits) for each hyperplane
  auto polygon_resolutions = compute_polygon_resolutions(
      N, F, v, w, b51, b57, b84, b86, b31, b36, b24, b27, f2min, f3min, f5min,
      f8min, f2max, f3max, f5max, f8max, verbose);

  hcpwa::AABB<8> aabb
      = {{0, 0, 0, 0, 0, 0, 0, 0},
         {static_cast<hcpwa::Float>(N), static_cast<hcpwa::Float>(N),
          static_cast<hcpwa::Float>(N), static_cast<hcpwa::Float>(N),
          static_cast<hcpwa::Float>(N), static_cast<hcpwa::Float>(N),
          static_cast<hcpwa::Float>(N), static_cast<hcpwa::Float>(N)}};
  hcpwa::AABB<2> aabb2d
      = {{0, 0}, {static_cast<hcpwa::Float>(N), static_cast<hcpwa::Float>(N)}};
  const auto aabb_bounds = hcpwa::AABBBounds(aabb);

  // Triangulate computed polygons and return computed triangles and 8D prisms
  auto polygon_prisms_result
      = compute_prisms_from_polygons(polygon_resolutions, aabb2d, verbose);

  // Intersects prisms and returns intersection points for poth phases
  auto& prisms31 = polygon_prisms_result.prisms31;
  auto& prisms36 = polygon_prisms_result.prisms36;
  auto& prisms24 = polygon_prisms_result.prisms24;
  auto& prisms27 = polygon_prisms_result.prisms27;
  auto& prisms58 = polygon_prisms_result.prisms58;
  auto& prisms51 = polygon_prisms_result.prisms51;
  auto& prisms57 = polygon_prisms_result.prisms57;
  auto& prisms84 = polygon_prisms_result.prisms84;
  auto& prisms86 = polygon_prisms_result.prisms86;
  auto& prisms23 = polygon_prisms_result.prisms23;
  auto& polygons58 = polygon_resolutions.resolution_58;
  auto& polygons23 = polygon_resolutions.resolution_23;

  // Intersects prisms and returns intersection points for both phases
  auto intersection_result = compute_intersection_points(
      prisms31, prisms36, prisms24, prisms27, prisms58, prisms51, prisms57,
      prisms84, prisms86, prisms23, polygons58, polygons23,
      static_cast<hcpwa::Float>(N), verbose);
  auto& intersection_points_phase0
      = intersection_result.intersection_points_phase0;
  auto& intersection_prism_indices_phase0
      = intersection_result.intersection_prism_indices_phase0;
  auto& intersection_points_phase1
      = intersection_result.intersection_points_phase1;
  auto& intersection_prism_indices_phase1
      = intersection_result.intersection_prism_indices_phase1;

  PolygonAreasVerticesResult result;
  // Intersection points and indices
  result.intersection_points_phase0 = intersection_points_phase0;
  result.intersection_prism_indices_phase0 = intersection_prism_indices_phase0;
  result.intersection_points_phase1 = intersection_points_phase1;
  result.intersection_prism_indices_phase1 = intersection_prism_indices_phase1;

  if (verbose) {
    std::cout << "result.intersection_points_phase0.size(): "
              << result.intersection_points_phase0.size() << '\n';
    std::cout << "result.intersection_prism_indices_phase0.size(): "
              << result.intersection_prism_indices_phase0.size() << '\n';
    std::cout << "result.intersection_points_phase1.size(): "
              << result.intersection_points_phase1.size() << '\n';
    std::cout << "result.intersection_prism_indices_phase1.size(): "
              << result.intersection_prism_indices_phase1.size() << '\n';
  }

  return result;
}

}  // namespace hcpwa
