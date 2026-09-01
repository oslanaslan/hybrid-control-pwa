#include <algo.hpp>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <format>
#include <limits>
#include <mutex>
#include <ostream>
#include <stdexcept>
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

namespace {

// ---- Shared block-decomposition helpers ------------------------------------
//
// Both overloads of compute_intersection_points() need the same thing: repackage
// what computend() already produced into the three coordinate blocks, keeping
// the COMPLETE vertex lists. The 8D assembly loops that follow truncate those
// lists to their first two entries; the block data is the un-truncated record.
//
// Nothing here computes geometry. See algo.hpp for the block structure and
// docs/barycentric_block_reduction_context.md part I for what the truncation
// costs.

BlockBounds3d blockBoundsOf(const std::vector<hcpwa::Vec<3>>& vertices,
                            int coord_count) {
  BlockBounds3d bounds;
  for (int d = 0; d < 3; ++d) {
    bounds.min[d] = 0.0;
    bounds.max[d] = 0.0;
  }
  for (int d = 0; d < coord_count; ++d) {
    bounds.min[d] = std::numeric_limits<double>::infinity();
    bounds.max[d] = -std::numeric_limits<double>::infinity();
  }
  for (const auto& vertex : vertices) {
    for (int d = 0; d < coord_count; ++d) {
      const double value = static_cast<double>(vertex[d]);
      bounds.min[d] = std::min(bounds.min[d], value);
      bounds.max[d] = std::max(bounds.max[d], value);
    }
  }
  return bounds;
}

// Blocks A and B: a fibre product of two planes over their shared coordinate,
// three-dimensional, six defining inequalities.
BlockRegions makePairBlock(
    std::array<int, 3> coords, std::array<int, 2> layer_ids,
    const std::vector<std::vector<size_t>>& indices,
    const std::vector<std::vector<hcpwa::Vec<3>>>& points, const char* label) {
  if (indices.size() != points.size()) {
    throw std::runtime_error(std::format(
        "compute_intersection_points: block {} has {} index tuples but {} "
        "vertex lists",
        label, indices.size(), points.size()));
  }
  BlockRegions block;
  block.coords = coords;
  block.coord_count = 3;
  block.layer_ids = layer_ids;
  block.layer_count = 2;
  block.triangle_ids.reserve(indices.size());
  block.vertices.reserve(points.size());
  block.bounds.reserve(points.size());
  for (size_t k = 0; k < indices.size(); ++k) {
    if (indices[k].size() != 2) {
      throw std::runtime_error(std::format(
          "compute_intersection_points: block {} region {} has {} simplex ids, "
          "expected 2",
          label, k, indices[k].size()));
    }
    // A stored block is the intersection of a 3D box with two prisms;
    // LinesToPoints<3> returns {} unless the result is full dimensional, so
    // every stored block has at least three vertices. Fewer would mean the
    // guarantee in algo.cpp changed.
    if (points[k].size() < 3) {
      throw std::runtime_error(std::format(
          "compute_intersection_points: block {} region {} has {} vertices, "
          "expected at least 3",
          label, k, points[k].size()));
    }
    block.triangle_ids.push_back({indices[k][0], indices[k][1]});
    block.vertices.push_back(points[k]);
    block.bounds.push_back(blockBoundsOf(points[k], 3));
  }
  return block;
}

// Block C: a single plane cell, two-dimensional. A triangle in the triangulated
// path, a general convex polygon in the polygon path, so the vertex count is
// only bounded below.
BlockRegions makeSingleBlock(
    std::array<int, 3> coords, std::array<int, 2> layer_ids,
    const std::vector<std::vector<hcpwa::Vec<2>>>& cells, const char* label) {
  BlockRegions block;
  block.coords = coords;
  block.coord_count = 2;
  block.layer_ids = layer_ids;
  block.layer_count = 1;
  block.triangle_ids.reserve(cells.size());
  block.vertices.reserve(cells.size());
  block.bounds.reserve(cells.size());
  for (size_t k = 0; k < cells.size(); ++k) {
    if (cells[k].size() < 3) {
      throw std::runtime_error(std::format(
          "compute_intersection_points: block {} region {} has {} vertices, "
          "expected at least 3",
          label, k, cells[k].size()));
    }
    std::vector<hcpwa::Vec<3>> vertices;
    vertices.reserve(cells[k].size());
    for (const auto& vertex : cells[k]) {
      hcpwa::Vec<3> point = {0, 0, 0};
      point[0] = vertex[0];
      point[1] = vertex[1];
      vertices.push_back(point);
    }
    block.triangle_ids.push_back({k, 0});
    block.bounds.push_back(blockBoundsOf(vertices, 2));
    block.vertices.push_back(std::move(vertices));
  }
  return block;
}

// Lemma 1 (step7_block_reduction.md section 4): the area index set is the FULL
// Cartesian product of the three block region sets, with no incompatible
// combinations, so M_i = M_A * M_B * M_C exactly. This is the one place that
// identity can be checked against the geometry rather than assumed, and
// everything downstream that maps an area id to a block triple depends on it.
void checkProductStructure(const std::array<BlockRegions, 3>& blocks,
                           size_t area_count, int phase) {
  const size_t m_a = blocks[0].vertices.size();
  const size_t m_b = blocks[1].vertices.size();
  const size_t m_c = blocks[2].vertices.size();
  if (m_a * m_b * m_c != area_count) {
    throw std::runtime_error(std::format(
        "compute_intersection_points: phase {} has {} areas but M_A*M_B*M_C = "
        "{}*{}*{} = {}. Lemma 1 says these must be equal; a mismatch means the "
        "block enumeration and the area assembly disagree.",
        phase, area_count, m_a, m_b, m_c, m_a * m_b * m_c));
  }
}

}  // namespace

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
    hcpwa::Float N, TriangleAreasOptions options) {
  const bool verbose = options.verbose;
  hcpwa::AABB<3> aabb3d = {{0, 0, 0}, {N, N, N}};
  const auto aabb3d_bounds = hcpwa::AABBBounds(aabb3d);

  auto computend = []<int Dim>(
                       std::array<int, Dim> dims,
                       const hcpwa::LineSet<Dim>& bounds,
                       const std::vector<hcpwa::LineSet<8>>& prisms0,
                       const std::vector<hcpwa::LineSet<8>>& prisms1,
                       const char* label,
                       bool verbose,
                       std::vector<std::vector<size_t>>& out_indices,
                       std::vector<std::vector<hcpwa::Vec<Dim>>>& out_points) {
    const std::size_t total_pairs = prisms0.size() * prisms1.size();
    std::size_t processed_pairs = 0;
    std::size_t non_empty_pairs = 0;
    const auto started_at = std::chrono::steady_clock::now();
    auto next_progress_at = started_at;
    auto report_progress = [&](bool force) {
      if (!verbose) {
        return;
      }
      const auto now = std::chrono::steady_clock::now();
      if (!force && now < next_progress_at) {
        return;
      }
      const double percent
          = total_pairs == 0
                ? 100.0
                : 100.0 * static_cast<double>(processed_pairs)
                      / static_cast<double>(total_pairs);
      const double elapsed_seconds
          = std::chrono::duration<double>(now - started_at).count();
      std::cerr << std::format(
          "Intersection {} progress: {:.1f}% ({}/{} pairs), non_empty={}, "
          "elapsed={:.1f}s\n",
          label, percent, processed_pairs, total_pairs, non_empty_pairs,
          elapsed_seconds);
      std::cerr.flush();
      next_progress_at = now + std::chrono::seconds(5);
    };
    report_progress(true);
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
          ++non_empty_pairs;
          out_points.push_back(intersection);
          // Store the indices of the prisms that form the intersection
          std::vector<size_t> prism_indices = {idx0, idx1};
          out_indices.push_back(prism_indices);
        }
        ++processed_pairs;
        if (processed_pairs % 10000 == 0) {
          report_progress(false);
        }
      }
    }
    report_progress(true);
  };

  if (verbose) {
    std::cout << std::format("Prism 31: {}", prisms31.size()) << std::endl;
    std::cout << std::format("Prism 36: {}", prisms36.size()) << std::endl;
    std::cout << std::format("Prism 24: {}", prisms24.size()) << std::endl;
    std::cout << std::format("Prism 27: {}", prisms27.size()) << std::endl;
  }

  std::vector<std::vector<size_t>> intersection_prism_indices_136;
  std::vector<std::vector<hcpwa::Vec<3>>> intersection_points_136;
  computend({0, 2, 5}, aabb3d_bounds, prisms31, prisms36, "136", verbose,
            intersection_prism_indices_136, intersection_points_136);

  std::vector<std::vector<size_t>> intersection_prism_indices_247;
  std::vector<std::vector<hcpwa::Vec<3>>> intersection_points_247;
  computend({1, 3, 6}, aabb3d_bounds, prisms24, prisms27, "247", verbose,
            intersection_prism_indices_247, intersection_points_247);

  std::vector<std::vector<size_t>> intersection_prism_indices_157;
  std::vector<std::vector<hcpwa::Vec<3>>> intersection_points_157;
  computend({0, 4, 6}, aabb3d_bounds, prisms51, prisms57, "157", verbose,
            intersection_prism_indices_157, intersection_points_157);

  std::vector<std::vector<size_t>> intersection_prism_indices_468;
  std::vector<std::vector<hcpwa::Vec<3>>> intersection_points_468;
  computend({3, 5, 7}, aabb3d_bounds, prisms84, prisms86, "468", verbose,
            intersection_prism_indices_468, intersection_points_468);


  // ---- Block decomposition -------------------------------------------------
  //
  // computend already produced exactly the block polytopes the reduction needs:
  // complete vertex lists and, in prism_indices, the pair of simplex ids that
  // generated each one. The prism arrays are built strictly parallel to the
  // triangle arrays (CalcPrism is applied to polygonsXY in order), so idx0 and
  // idx1 are directly usable as triangle ids.

  auto triangle_cells
      = [](const std::vector<hcpwa::TriangleWithUniqueVertices>& triangles) {
          std::vector<std::vector<hcpwa::Vec<2>>> cells;
          cells.reserve(triangles.size());
          for (const auto& triangle : triangles) {
            std::vector<hcpwa::Vec<2>> cell;
            cell.reserve(triangle.size());
            for (size_t v = 0; v < triangle.size(); ++v) {
              cell.push_back(triangle[v]);
            }
            cells.push_back(std::move(cell));
          }
          return cells;
        };

  std::array<BlockRegions, 3> blocks_phase0
      = {makePairBlock({0, 2, 5}, {0, 1}, intersection_prism_indices_136,
                       intersection_points_136, "136"),
         makePairBlock({1, 3, 6}, {2, 3}, intersection_prism_indices_247,
                       intersection_points_247, "247"),
         makeSingleBlock({4, 7, -1}, {4, -1}, triangle_cells(triangles58),
                         "58")};
  std::array<BlockRegions, 3> blocks_phase1
      = {makePairBlock({0, 4, 6}, {0, 1}, intersection_prism_indices_157,
                       intersection_points_157, "157"),
         makePairBlock({3, 5, 7}, {2, 3}, intersection_prism_indices_468,
                       intersection_points_468, "468"),
         makeSingleBlock({1, 2, -1}, {4, -1}, triangle_cells(triangles23),
                         "23")};

  if (verbose) {
    std::cout << "Intersection counted:" << std::endl;
    std::cout << "\t136 count: " << intersection_points_136.size() << std::endl;
    std::cout << "\t247 count: " << intersection_points_247.size() << std::endl;
    std::cout << "\t58 count: " << triangles58.size() << std::endl;
    std::cout << "\t157 count: " << intersection_points_157.size() << std::endl;
    std::cout << "\t468 count: " << intersection_points_468.size() << std::endl;
    std::cout << "\t23 count: " << triangles23.size() << std::endl;

    // Block cardinalities and the pair counts they came from. A block that is
    // non-empty but not full dimensional is dropped by LinesToPoints, so a
    // large gap between pairs tried and blocks kept is worth noticing.
    auto report_block = [](const char* label, const BlockRegions& block,
                           size_t pairs_tried) {
      size_t min_vertices = std::numeric_limits<size_t>::max();
      size_t max_vertices = 0;
      size_t total_vertices = 0;
      for (const auto& vertices : block.vertices) {
        min_vertices = std::min(min_vertices, vertices.size());
        max_vertices = std::max(max_vertices, vertices.size());
        total_vertices += vertices.size();
      }
      if (block.vertices.empty()) {
        min_vertices = 0;
      }
      std::cout << std::format(
          "\tblock {}: regions={} (of {} pairs tried), vertices min/mean/max = "
          "{}/{:.2f}/{}\n",
          label, block.vertices.size(), pairs_tried, min_vertices,
          block.vertices.empty()
              ? 0.0
              : static_cast<double>(total_vertices)
                    / static_cast<double>(block.vertices.size()),
          max_vertices);
    };
    report_block("A  136", blocks_phase0[0], prisms31.size() * prisms36.size());
    report_block("B  247", blocks_phase0[1], prisms24.size() * prisms27.size());
    report_block("C   58", blocks_phase0[2], triangles58.size());
    report_block("A' 157", blocks_phase1[0], prisms51.size() * prisms57.size());
    report_block("B' 468", blocks_phase1[1], prisms84.size() * prisms86.size());
    report_block("C'  23", blocks_phase1[2], triangles23.size());
  }

  std::vector<std::vector<size_t>> intersection_prism_indices_phase0;
  std::vector<std::vector<hcpwa::Vec<8>>> intersection_points_phase0;
  std::vector<std::vector<size_t>> intersection_prism_indices_phase1;
  std::vector<std::vector<hcpwa::Vec<8>>> intersection_points_phase1;

  const std::size_t total_phase0_areas
      = intersection_points_136.size() * intersection_prism_indices_247.size()
        * triangles58.size();
  std::size_t processed_phase0_areas = 0;
  const auto phase0_started_at = std::chrono::steady_clock::now();
  auto next_phase0_progress_at = phase0_started_at;
  auto report_phase0_progress = [&](bool force) {
    if (!verbose) {
      return;
    }
    const auto now = std::chrono::steady_clock::now();
    if (!force && now < next_phase0_progress_at) {
      return;
    }
    const double percent
        = total_phase0_areas == 0
              ? 100.0
              : 100.0 * static_cast<double>(processed_phase0_areas)
                    / static_cast<double>(total_phase0_areas);
    const double elapsed_seconds
        = std::chrono::duration<double>(now - phase0_started_at).count();
    std::cerr << std::format(
        "Phase0 area assembly progress: {:.1f}% ({}/{} areas), vertices={}, "
        "elapsed={:.1f}s\n",
        percent, processed_phase0_areas, total_phase0_areas,
        intersection_points_phase0.size(), elapsed_seconds);
    std::cerr.flush();
    next_phase0_progress_at = now + std::chrono::seconds(5);
  };
  if (options.build_8d_regions) {
    // The loops below take their bounds from the 2-element prism-index lists
    // while indexing the block vertex lists, so they keep only the first two
    // vertices of each 3D block: 12 vertices per area instead of the 108-192
    // an area actually has. The truncation is left in place on purpose. Every
    // consumer of it is either being migrated to the block data or is slated
    // for deletion, and correcting the bounds here would cost roughly 13 GB of
    // vertices for the piecewise path alone. See
    // docs/barycentric_block_reduction_context.md parts I and IV.2.
    std::cerr << "compute_intersection_points: building 8D area vertices with "
                 "the known vertex-set truncation (12 of 108-192 vertices per "
                 "area). Callers that need complete geometry must use the "
                 "block data instead. See "
                 "docs/barycentric_block_reduction_context.md part I.\n";
    std::cerr.flush();
  }

  report_phase0_progress(true);

  for (size_t i136 = 0; i136 < intersection_points_136.size(); i136++) {
    for (size_t i247 = 0; i247 < intersection_prism_indices_247.size();
         i247++) {
      // auto t_start = std::chrono::high_resolution_clock::now();
      for (size_t i58 = 0; i58 < triangles58.size(); i58++) {
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

        // The simplex-id tuple above is always built: the border path indexes
        // areas by it and it costs 5 integers per area. Only the explicit 8D
        // vertex list below is optional, and it is the truncated one.
        if (options.build_8d_regions) {
          intersection_points_phase0.emplace_back();
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
        ++processed_phase0_areas;
        if (processed_phase0_areas % 1000 == 0) {
          report_phase0_progress(false);
        }
      }
      // auto t_end = std::chrono::high_resolution_clock::now();
      // std::chrono::duration<double> t_diff = t_end - t_start;
      // std::cout << "[Timing] i136=" << i136 << ", i247=" << i247 << ": " <<
      // t_diff.count() << "s" << std::endl;
    }
  }
  report_phase0_progress(true);

  const std::size_t total_phase1_areas
      = intersection_points_157.size() * intersection_prism_indices_468.size()
        * triangles23.size();
  std::size_t processed_phase1_areas = 0;
  const auto phase1_started_at = std::chrono::steady_clock::now();
  auto next_phase1_progress_at = phase1_started_at;
  auto report_phase1_progress = [&](bool force) {
    if (!verbose) {
      return;
    }
    const auto now = std::chrono::steady_clock::now();
    if (!force && now < next_phase1_progress_at) {
      return;
    }
    const double percent
        = total_phase1_areas == 0
              ? 100.0
              : 100.0 * static_cast<double>(processed_phase1_areas)
                    / static_cast<double>(total_phase1_areas);
    const double elapsed_seconds
        = std::chrono::duration<double>(now - phase1_started_at).count();
    std::cerr << std::format(
        "Phase1 area assembly progress: {:.1f}% ({}/{} areas), vertices={}, "
        "elapsed={:.1f}s\n",
        percent, processed_phase1_areas, total_phase1_areas,
        intersection_points_phase1.size(), elapsed_seconds);
    std::cerr.flush();
    next_phase1_progress_at = now + std::chrono::seconds(5);
  };
  report_phase1_progress(true);

  for (size_t i157 = 0; i157 < intersection_points_157.size(); i157++) {
    for (size_t i468 = 0; i468 < intersection_prism_indices_468.size();
         i468++) {
      for (size_t i23 = 0; i23 < triangles23.size(); i23++) {
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

        // See the phase-0 loop: the tuple is always built, the truncated 8D
        // vertex list is optional.
        if (options.build_8d_regions) {
          intersection_points_phase1.emplace_back();
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
        ++processed_phase1_areas;
        if (processed_phase1_areas % 1000 == 0) {
          report_phase1_progress(false);
        }
      }
    }
  }
  report_phase1_progress(true);

  // Lemma 1 (step7_block_reduction.md section 4): the area index set is the
  // FULL Cartesian product of the three block region sets, with no
  // incompatible combinations, so M_i = M_A * M_B * M_C exactly. This is the
  // one place that identity can be checked against the geometry rather than
  // assumed, and everything downstream that maps an area id to a block triple
  // depends on it.
  checkProductStructure(blocks_phase0,
                        intersection_prism_indices_phase0.size(), 0);
  checkProductStructure(blocks_phase1,
                        intersection_prism_indices_phase1.size(), 1);

  // The 8D box of an area is the concatenation of its three block boxes.
  // Exact, not an over-approximation: the area is a product and is
  // unconstrained in the coordinates of the other two blocks, so
  //   min_{n in area} n_r = min_{nu in block polytope} nu_r  for r in I_block.
  // Composing them here costs O(M_A*M_B*M_C) writes but no geometry, against
  // O(sum of area vertex counts) for the old per-area scan over 8D vertices.
  auto compose_area_bounds = [](const std::array<BlockRegions, 3>& blocks) {
    const size_t m_a = blocks[0].bounds.size();
    const size_t m_b = blocks[1].bounds.size();
    const size_t m_c = blocks[2].bounds.size();
    std::vector<AreaBounds8d> out;
    out.reserve(m_a * m_b * m_c);
    for (size_t j_a = 0; j_a < m_a; ++j_a) {
      for (size_t j_b = 0; j_b < m_b; ++j_b) {
        for (size_t j_c = 0; j_c < m_c; ++j_c) {
          const std::array<size_t, 3> ids = {j_a, j_b, j_c};
          AreaBounds8d bounds;
          for (int b = 0; b < 3; ++b) {
            const BlockRegions& block = blocks[b];
            const BlockBounds3d& block_bounds = block.bounds[ids[b]];
            for (int d = 0; d < block.coord_count; ++d) {
              bounds.min[block.coords[d]] = block_bounds.min[d];
              bounds.max[block.coords[d]] = block_bounds.max[d];
            }
          }
          out.push_back(bounds);
        }
      }
    }
    return out;
  };

  PhaseIntersectionResult result;
  result.intersection_prism_indices_phase0
      = std::move(intersection_prism_indices_phase0);
  result.intersection_points_phase0 = std::move(intersection_points_phase0);
  result.intersection_prism_indices_phase1
      = std::move(intersection_prism_indices_phase1);
  result.intersection_points_phase1 = std::move(intersection_points_phase1);
  result.area_bounds_phase0 = compose_area_bounds(blocks_phase0);
  result.area_bounds_phase1 = compose_area_bounds(blocks_phase1);
  result.blocks_phase0 = std::move(blocks_phase0);
  result.blocks_phase1 = std::move(blocks_phase1);
  return result;
}

CommonRefinementResult compute_common_refinement_area_vertices(
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
    const std::vector<std::vector<size_t>>& phase0_area_prism_indices,
    const std::vector<std::vector<size_t>>& phase1_area_prism_indices,
    const std::vector<AreaBounds8d>& phase0_area_bounds,
    const std::vector<AreaBounds8d>& phase1_area_bounds,
    hcpwa::Float N, bool verbose) {
  constexpr int kSpaceDim = 8;
  constexpr double kBoundsTol = 1e-9;
  constexpr double kRelativeRowTol = 1e-12;
  hcpwa::AABB<8> aabb
      = {{0, 0, 0, 0, 0, 0, 0, 0}, {N, N, N, N, N, N, N, N}};
  const auto aabb_bounds = hcpwa::AABBBounds(aabb);

  auto require_tuple = [](const std::vector<size_t>& tuple, int phase,
                          size_t area_id) {
    if (tuple.size() != 5) {
      throw std::runtime_error(std::format(
          "compute_common_refinement_area_vertices: phase {} area {} has {} "
          "prism ids, expected 5",
          phase, area_id, tuple.size()));
    }
  };

  auto require_index = [](size_t idx, size_t size, int phase, size_t area_id,
                          int layer_id, const char* prism_name) {
    if (idx >= size) {
      throw std::runtime_error(std::format(
          "compute_common_refinement_area_vertices: phase {} area {} layer {} "
          "uses {} index {} but only {} prisms exist",
          phase, area_id, layer_id, prism_name, idx, size));
    }
  };

  auto append_prism = [](hcpwa::LineSet<8>& dst,
                         const hcpwa::LineSet<8>& src) {
    dst.insert(dst.end(), src.begin(), src.end());
  };

  if (phase0_area_bounds.size() != phase0_area_prism_indices.size()) {
    throw std::runtime_error(std::format(
        "compute_common_refinement_area_vertices: phase 0 has {} area bounds "
        "but {} prism-index tuples",
        phase0_area_bounds.size(), phase0_area_prism_indices.size()));
  }
  if (phase1_area_bounds.size() != phase1_area_prism_indices.size()) {
    throw std::runtime_error(std::format(
        "compute_common_refinement_area_vertices: phase 1 has {} area bounds "
        "but {} prism-index tuples",
        phase1_area_bounds.size(), phase1_area_prism_indices.size()));
  }

  auto boxes_have_full_dimensional_overlap
      = [](const AreaBounds8d& lhs, const AreaBounds8d& rhs, double tol) {
    for (int d = 0; d < kSpaceDim; ++d) {
      const double overlap
          = std::min(lhs.max[d], rhs.max[d])
            - std::max(lhs.min[d], rhs.min[d]);
      if (overlap <= tol) {
        return false;
      }
    }
    return true;
  };

  std::vector<std::array<size_t, 5>> phase0_prism_ids_cache;
  std::vector<std::array<size_t, 5>> phase1_prism_ids_cache;
  phase0_prism_ids_cache.reserve(phase0_area_prism_indices.size());
  phase1_prism_ids_cache.reserve(phase1_area_prism_indices.size());

  // The boxes are supplied by the caller now; no copy, no recomputation.
  const std::vector<AreaBounds8d>& phase0_bounds = phase0_area_bounds;
  const std::vector<AreaBounds8d>& phase1_bounds = phase1_area_bounds;

  for (size_t phase0_area_id = 0; phase0_area_id < phase0_bounds.size();
       ++phase0_area_id) {
    const auto& phase0_tuple = phase0_area_prism_indices[phase0_area_id];
    require_tuple(phase0_tuple, 0, phase0_area_id);
    const std::array<size_t, 5> phase0_prism_ids = {
        phase0_tuple[0], phase0_tuple[1], phase0_tuple[2], phase0_tuple[3],
        phase0_tuple[4]};
    require_index(
        phase0_prism_ids[0], prisms31.size(), 0, phase0_area_id, 0, "31");
    require_index(
        phase0_prism_ids[1], prisms36.size(), 0, phase0_area_id, 1, "36");
    require_index(
        phase0_prism_ids[2], prisms24.size(), 0, phase0_area_id, 2, "24");
    require_index(
        phase0_prism_ids[3], prisms27.size(), 0, phase0_area_id, 3, "27");
    require_index(
        phase0_prism_ids[4], prisms58.size(), 0, phase0_area_id, 4, "58");
    phase0_prism_ids_cache.push_back(phase0_prism_ids);
  }

  for (size_t phase1_area_id = 0; phase1_area_id < phase1_bounds.size();
       ++phase1_area_id) {
    const auto& phase1_tuple = phase1_area_prism_indices[phase1_area_id];
    require_tuple(phase1_tuple, 1, phase1_area_id);
    const std::array<size_t, 5> phase1_prism_ids = {
        phase1_tuple[0], phase1_tuple[1], phase1_tuple[2], phase1_tuple[3],
        phase1_tuple[4]};
    require_index(
        phase1_prism_ids[0], prisms51.size(), 1, phase1_area_id, 0, "51");
    require_index(
        phase1_prism_ids[1], prisms57.size(), 1, phase1_area_id, 1, "57");
    require_index(
        phase1_prism_ids[2], prisms84.size(), 1, phase1_area_id, 2, "84");
    require_index(
        phase1_prism_ids[3], prisms86.size(), 1, phase1_area_id, 3, "86");
    require_index(
        phase1_prism_ids[4], prisms23.size(), 1, phase1_area_id, 4, "23");
    phase1_prism_ids_cache.push_back(phase1_prism_ids);
  }

  auto has_positive_width = [](const AreaBounds8d& bounds, double tol) {
    for (int d = 0; d < kSpaceDim; ++d) {
      if (bounds.max[d] - bounds.min[d] <= tol) {
        return false;
      }
    }
    return true;
  };
  std::vector<size_t> eligible_phase0;
  std::vector<size_t> eligible_phase1;
  eligible_phase0.reserve(phase0_bounds.size());
  eligible_phase1.reserve(phase1_bounds.size());
  for (size_t id = 0; id < phase0_bounds.size(); ++id) {
    if (has_positive_width(phase0_bounds[id], kBoundsTol)) {
      eligible_phase0.push_back(id);
    }
  }
  for (size_t id = 0; id < phase1_bounds.size(); ++id) {
    if (has_positive_width(phase1_bounds[id], kBoundsTol)) {
      eligible_phase1.push_back(id);
    }
  }

  // An area dropped here never reaches the common refinement, so the border LP
  // simply stops constraining it. When the boxes were derived from truncated
  // 8D vertex lists this discarded 99.99% of phase-0 areas without a word.
  // Report it unconditionally: a high drop rate is either a genuinely
  // degenerate arrangement or a regression in how the boxes are built, and
  // both are things we want to hear about rather than discover later.
  auto report_drop_rate = [](int phase, size_t eligible, size_t total) {
    if (total == 0 || eligible == total) {
      return;
    }
    const double dropped_percent = 100.0
                                   * static_cast<double>(total - eligible)
                                   / static_cast<double>(total);
    std::cerr << std::format(
        "compute_common_refinement_area_vertices: phase {} dropped {} of {} "
        "areas ({:.2f}%) as degenerate before pairing.{}\n",
        phase, total - eligible, total, dropped_percent,
        dropped_percent > 50.0
            ? " That is most of the partition; if the per-area boxes are"
              " built from complete block vertex sets this should be rare."
            : "");
    std::cerr.flush();
  };
  report_drop_rate(0, eligible_phase0.size(), phase0_bounds.size());
  report_drop_rate(1, eligible_phase1.size(), phase1_bounds.size());

  auto choose_sweep_axis_by_exact_overlap =
      [&](const std::vector<size_t>& phase0_ids,
          const std::vector<size_t>& phase1_ids) {
        if (phase0_ids.empty() || phase1_ids.empty()) {
          return std::pair<int, std::uint64_t>{0, 0};
        }
        int best_axis = 0;
        std::uint64_t best_overlap = std::numeric_limits<std::uint64_t>::max();
        for (int axis = 0; axis < kSpaceDim; ++axis) {
          std::vector<double> min1;
          std::vector<double> max1;
          min1.reserve(phase1_ids.size());
          max1.reserve(phase1_ids.size());
          for (size_t id : phase1_ids) {
            min1.push_back(phase1_bounds[id].min[axis]);
            max1.push_back(phase1_bounds[id].max[axis]);
          }
          std::sort(min1.begin(), min1.end());
          std::sort(max1.begin(), max1.end());

          std::uint64_t overlap_count = 0;
          for (size_t id : phase0_ids) {
            const double right = phase0_bounds[id].max[axis] - kBoundsTol;
            const double left = phase0_bounds[id].min[axis] + kBoundsTol;
            const std::size_t started
                = static_cast<std::size_t>(std::lower_bound(
                      min1.begin(), min1.end(), right)
                                           - min1.begin());
            const std::size_t ended
                = static_cast<std::size_t>(std::upper_bound(
                      max1.begin(), max1.end(), left)
                                           - max1.begin());
            if (started > ended) {
              overlap_count += static_cast<std::uint64_t>(started - ended);
            }
          }
          if (overlap_count < best_overlap) {
            best_overlap = overlap_count;
            best_axis = axis;
          }
        }
        return std::pair<int, std::uint64_t>{best_axis, best_overlap};
      };

  const auto [sweep_axis, axis_overlap_estimate]
      = choose_sweep_axis_by_exact_overlap(eligible_phase0, eligible_phase1);

  std::vector<size_t> phase1_by_min = eligible_phase1;
  std::vector<size_t> phase1_by_max = eligible_phase1;
  std::sort(phase1_by_min.begin(), phase1_by_min.end(),
            [&](size_t lhs, size_t rhs) {
              if (phase1_bounds[lhs].min[sweep_axis]
                  != phase1_bounds[rhs].min[sweep_axis]) {
                return phase1_bounds[lhs].min[sweep_axis]
                       < phase1_bounds[rhs].min[sweep_axis];
              }
              return lhs < rhs;
            });
  std::sort(phase1_by_max.begin(), phase1_by_max.end(),
            [&](size_t lhs, size_t rhs) {
              if (phase1_bounds[lhs].max[sweep_axis]
                  != phase1_bounds[rhs].max[sweep_axis]) {
                return phase1_bounds[lhs].max[sweep_axis]
                       < phase1_bounds[rhs].max[sweep_axis];
              }
              return lhs < rhs;
            });

  std::vector<double> phase1_min_sorted;
  std::vector<double> phase1_max_sorted;
  phase1_min_sorted.reserve(phase1_by_min.size());
  phase1_max_sorted.reserve(phase1_by_max.size());
  for (size_t id : phase1_by_min) {
    phase1_min_sorted.push_back(phase1_bounds[id].min[sweep_axis]);
  }
  for (size_t id : phase1_by_max) {
    phase1_max_sorted.push_back(phase1_bounds[id].max[sweep_axis]);
  }

  std::atomic<std::uint64_t> sweep_candidates{0};
  std::atomic<std::uint64_t> processed_phase0{0};
  const auto progress_started_at = std::chrono::steady_clock::now();
  // Mutated from every worker below, so it cannot be a plain time_point. Stored
  // as a duration since progress_started_at because steady_clock::time_point is
  // not trivially atomic.
  std::atomic<std::int64_t> next_progress_ns{0};
  std::mutex progress_mutex;
  auto report_box_progress = [&](bool force) {
    if (!verbose) {
      return;
    }
    const auto now = std::chrono::steady_clock::now();
    const std::int64_t now_ns
        = std::chrono::duration_cast<std::chrono::nanoseconds>(
              now - progress_started_at)
              .count();
    if (!force && now_ns < next_progress_ns.load(std::memory_order_relaxed)) {
      return;
    }
    std::lock_guard<std::mutex> progress_lock(progress_mutex);
    if (!force && now_ns < next_progress_ns.load(std::memory_order_relaxed)) {
      return;
    }
    const std::uint64_t processed = processed_phase0.load();
    const std::uint64_t candidates = sweep_candidates.load();
    const double percent = eligible_phase0.empty()
                               ? 100.0
                               : 100.0 * static_cast<double>(processed)
                                     / static_cast<double>(eligible_phase0.size());
    const double elapsed_seconds
        = std::chrono::duration<double>(now - progress_started_at).count();
    std::cerr << std::format(
        "Common refinement box scan: {:.1f}% ({}/{} areas), candidates={}, "
        "elapsed={:.1f}s\n",
        percent, processed, eligible_phase0.size(), candidates, elapsed_seconds);
    std::cerr.flush();
    next_progress_ns.store(
        now_ns + std::chrono::nanoseconds(std::chrono::seconds(5)).count(),
        std::memory_order_relaxed);
  };

  const unsigned int hw = std::thread::hardware_concurrency();
  const std::size_t worker_count = std::min<std::size_t>(
      std::max<unsigned int>(1, hw),
      std::max<std::size_t>(std::size_t{1}, eligible_phase0.size()));
  const std::size_t chunk_size = worker_count == 0
                                     ? 0
                                     : (eligible_phase0.size() + worker_count - 1)
                                           / worker_count;
  std::vector<std::vector<std::pair<size_t, size_t>>> worker_matches(
      worker_count);
  std::vector<std::thread> workers;
  workers.reserve(worker_count);

  report_box_progress(true);
  for (std::size_t worker_id = 0; worker_id < worker_count; ++worker_id) {
    const std::size_t begin = worker_id * chunk_size;
    const std::size_t end = std::min(eligible_phase0.size(), begin + chunk_size);
    workers.emplace_back([&, worker_id, begin, end]() {
      auto& local_matches = worker_matches[worker_id];
      for (std::size_t idx = begin; idx < end; ++idx) {
        const size_t phase0_area_id = eligible_phase0[idx];
        const auto& p0 = phase0_bounds[phase0_area_id];
        const double upper_bound_min = p0.max[sweep_axis] - kBoundsTol;
        const double lower_bound_max = p0.min[sweep_axis] + kBoundsTol;

        const std::size_t prefix_end = static_cast<std::size_t>(
            std::lower_bound(phase1_min_sorted.begin(), phase1_min_sorted.end(),
                             upper_bound_min)
            - phase1_min_sorted.begin());
        const std::size_t suffix_begin = static_cast<std::size_t>(
            std::upper_bound(phase1_max_sorted.begin(), phase1_max_sorted.end(),
                             lower_bound_max)
            - phase1_max_sorted.begin());

        if (prefix_end <= (phase1_max_sorted.size() - suffix_begin)) {
          for (std::size_t i = 0; i < prefix_end; ++i) {
            ++sweep_candidates;
            const size_t phase1_area_id = phase1_by_min[i];
            const auto& p1 = phase1_bounds[phase1_area_id];
            if (p1.max[sweep_axis] <= lower_bound_max) {
              continue;
            }
            if (boxes_have_full_dimensional_overlap(p0, p1, kBoundsTol)) {
              local_matches.emplace_back(phase0_area_id, phase1_area_id);
            }
          }
        } else {
          for (std::size_t i = suffix_begin; i < phase1_max_sorted.size(); ++i) {
            ++sweep_candidates;
            const size_t phase1_area_id = phase1_by_max[i];
            const auto& p1 = phase1_bounds[phase1_area_id];
            if (p1.min[sweep_axis] >= upper_bound_min) {
              continue;
            }
            if (boxes_have_full_dimensional_overlap(p0, p1, kBoundsTol)) {
              local_matches.emplace_back(phase0_area_id, phase1_area_id);
            }
          }
        }

        ++processed_phase0;
        report_box_progress(false);
      }
    });
  }
  for (auto& worker : workers) {
    worker.join();
  }
  report_box_progress(true);

  std::vector<std::pair<size_t, size_t>> candidate_pairs;
  std::size_t total_matches = 0;
  for (const auto& local : worker_matches) {
    total_matches += local.size();
  }
  candidate_pairs.reserve(total_matches);
  for (auto& local : worker_matches) {
    candidate_pairs.insert(candidate_pairs.end(), local.begin(), local.end());
  }
  std::sort(candidate_pairs.begin(), candidate_pairs.end());
  candidate_pairs.erase(
      std::unique(candidate_pairs.begin(), candidate_pairs.end()),
      candidate_pairs.end());

  CommonRefinementResult result;
  std::uint64_t cdd_calls = 0;
  std::size_t processed_candidates = 0;
  auto next_cdd_progress_at = std::chrono::steady_clock::now();
  auto report_cdd_progress = [&](bool force) {
    if (!verbose) {
      return;
    }
    const auto now = std::chrono::steady_clock::now();
    if (!force && now < next_cdd_progress_at) {
      return;
    }
    const double percent = candidate_pairs.empty()
                               ? 100.0
                               : 100.0 * static_cast<double>(processed_candidates)
                                     / static_cast<double>(candidate_pairs.size());
    const double elapsed_seconds
        = std::chrono::duration<double>(now - progress_started_at).count();
    std::cerr << std::format(
        "Common refinement cdd pass: {:.1f}% ({}/{} pairs), cdd_calls={}, "
        "areas={}, elapsed={:.1f}s\n",
        percent, processed_candidates, candidate_pairs.size(), cdd_calls,
        result.areas.size(), elapsed_seconds);
    std::cerr.flush();
    next_cdd_progress_at = now + std::chrono::seconds(5);
  };

  report_cdd_progress(true);
  for (const auto& [phase0_area_id, phase1_area_id] : candidate_pairs) {
    const auto& phase0_prism_ids = phase0_prism_ids_cache[phase0_area_id];
    const auto& phase1_prism_ids = phase1_prism_ids_cache[phase1_area_id];
    hcpwa::LineSet<8> combined_prisms = aabb_bounds;
    append_prism(combined_prisms, prisms31[phase0_prism_ids[0]]);
    append_prism(combined_prisms, prisms36[phase0_prism_ids[1]]);
    append_prism(combined_prisms, prisms24[phase0_prism_ids[2]]);
    append_prism(combined_prisms, prisms27[phase0_prism_ids[3]]);
    append_prism(combined_prisms, prisms58[phase0_prism_ids[4]]);
    append_prism(combined_prisms, prisms51[phase1_prism_ids[0]]);
    append_prism(combined_prisms, prisms57[phase1_prism_ids[1]]);
    append_prism(combined_prisms, prisms84[phase1_prism_ids[2]]);
    append_prism(combined_prisms, prisms86[phase1_prism_ids[3]]);
    append_prism(combined_prisms, prisms23[phase1_prism_ids[4]]);

    for (auto& row : combined_prisms) {
      double normal_scale = 0.0;
      for (int col = 0; col < kSpaceDim; ++col) {
        normal_scale
            = std::max(normal_scale, std::abs(static_cast<double>(row[col])));
      }
      if (normal_scale == 0.0) {
        throw std::runtime_error(
            "compute_common_refinement_area_vertices: inequality has zero "
            "normal");
      }
      const double zero_threshold = kRelativeRowTol * normal_scale;
      for (int col = 0; col <= kSpaceDim; ++col) {
        const double value = static_cast<double>(row[col]);
        row[col] = std::abs(value) <= zero_threshold ? 0.0
                                                     : value / normal_scale;
      }
    }

    ++cdd_calls;
    std::vector<hcpwa::Vec<8>> vertices;
    try {
      vertices = hcpwa::LinesToPoints<8>(combined_prisms);
    } catch (const std::exception& error) {
      std::cerr << "Common refinement cdd failure\n";
      std::cerr << "  reason: " << error.what() << '\n';
      std::cerr << "  sweep_axis: " << sweep_axis << '\n';
      std::cerr << "  sweep_candidates: " << sweep_candidates.load() << '\n';
      std::cerr << "  full_box_matches: " << candidate_pairs.size() << '\n';
      std::cerr << "  cdd_call: " << cdd_calls << '\n';
      std::cerr << "  phase0_area_id: " << phase0_area_id << '\n';
      std::cerr << "  phase1_area_id: " << phase1_area_id << '\n';

      auto print_ids = [](const char* label,
                          const std::array<size_t, 5>& ids) {
        std::cerr << "  " << label << ": [" << ids[0] << ", " << ids[1]
                  << ", " << ids[2] << ", " << ids[3] << ", " << ids[4]
                  << "]\n";
      };
      print_ids("phase0_prism_ids [31,36,24,27,58]", phase0_prism_ids);
      print_ids("phase1_prism_ids [51,57,84,86,23]", phase1_prism_ids);

      auto print_bounds = [&](const char* label, const AreaBounds8d& bounds) {
        std::cerr << "  " << label << "_min: [";
        for (int d = 0; d < kSpaceDim; ++d) {
          std::cerr << (d == 0 ? "" : ", ") << bounds.min[d];
        }
        std::cerr << "]\n";
        std::cerr << "  " << label << "_max: [";
        for (int d = 0; d < kSpaceDim; ++d) {
          std::cerr << (d == 0 ? "" : ", ") << bounds.max[d];
        }
        std::cerr << "]\n";
      };
      print_bounds("phase0_bounds", phase0_bounds[phase0_area_id]);
      print_bounds("phase1_bounds", phase1_bounds[phase1_area_id]);

      size_t non_finite_values = 0;
      double max_abs_coefficient = 0.0;
      double min_nonzero_abs_coefficient
          = std::numeric_limits<double>::infinity();
      for (const auto& row : combined_prisms) {
        for (int col = 0; col <= kSpaceDim; ++col) {
          const double value = static_cast<double>(row[col]);
          if (!std::isfinite(value)) {
            ++non_finite_values;
            continue;
          }
          const double abs_value = std::abs(value);
          max_abs_coefficient = std::max(max_abs_coefficient, abs_value);
          if (abs_value > 0.0) {
            min_nonzero_abs_coefficient
                = std::min(min_nonzero_abs_coefficient, abs_value);
          }
        }
      }
      std::cerr << "  inequality_rows: " << combined_prisms.size() << '\n';
      std::cerr << "  non_finite_values: " << non_finite_values << '\n';
      std::cerr << "  min_nonzero_abs_coefficient: "
                << min_nonzero_abs_coefficient << '\n';
      std::cerr << "  max_abs_coefficient: " << max_abs_coefficient << '\n';
      std::cerr << "  inequalities [a0..a7, constant]:\n";
      for (size_t row_id = 0; row_id < combined_prisms.size(); ++row_id) {
        std::cerr << "    " << row_id << ": [";
        for (int col = 0; col <= kSpaceDim; ++col) {
          std::cerr << (col == 0 ? "" : ", ")
                    << static_cast<double>(combined_prisms[row_id][col]);
        }
        std::cerr << "]\n";
      }
      throw std::runtime_error(std::format(
          "common refinement cdd failed for phase0 area {}, phase1 area {}, "
          "cdd call {}: {}",
          phase0_area_id, phase1_area_id, cdd_calls, error.what()));
    }
    if (!vertices.empty()) {
      CommonRefinementArea area;
      area.phase0_area_id = phase0_area_id;
      area.phase1_area_id = phase1_area_id;
      area.phase0_prism_indices = phase0_prism_ids;
      area.phase1_prism_indices = phase1_prism_ids;
      area.vertices = std::move(vertices);
      result.areas.push_back(std::move(area));
    }
    ++processed_candidates;
    report_cdd_progress(false);
  }
  report_cdd_progress(true);

  if (verbose) {
    std::cout << std::format("Common refinement workers: {}\n", worker_count);
    std::cout << std::format("Common refinement eligible phase0 boxes: {} of {}\n",
                             eligible_phase0.size(), phase0_bounds.size());
    std::cout << std::format("Common refinement eligible phase1 boxes: {} of {}\n",
                             eligible_phase1.size(), phase1_bounds.size());
    std::cout << std::format(
        "Common refinement selected-axis 1D overlap estimate: {}\n",
        axis_overlap_estimate);
    std::cout << std::format("Common refinement sweep axis: {}\n", sweep_axis);
    std::cout << std::format("Common refinement sweep candidates: {}\n",
                             sweep_candidates.load());
    std::cout << std::format("Common refinement 8D-box matches: {}\n",
                             candidate_pairs.size());
    std::cout << std::format("Common refinement cdd calls: {}\n", cdd_calls);
    std::cout << "Common refinement areas count: " << result.areas.size()
              << '\n';
  }
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
  hcpwa::Float N, TriangleAreasOptions options) {
  const bool verbose = options.verbose;
hcpwa::AABB<3> aabb3d = {{0, 0, 0}, {N, N, N}};
const auto aabb3d_bounds = hcpwa::AABBBounds(aabb3d);

auto computend = []<int Dim>(
                     std::array<int, Dim> dims,
                     const hcpwa::LineSet<Dim>& bounds,
                     const std::vector<hcpwa::LineSet<8>>& prisms0,
                     const std::vector<hcpwa::LineSet<8>>& prisms1,
                     const char* label,
                     bool verbose,
                     std::vector<std::vector<size_t>>& out_indices,
                     std::vector<std::vector<hcpwa::Vec<Dim>>>& out_points) {
  const std::size_t total_pairs = prisms0.size() * prisms1.size();
  std::size_t processed_pairs = 0;
  std::size_t non_empty_pairs = 0;
  const auto started_at = std::chrono::steady_clock::now();
  auto next_progress_at = started_at;
  auto report_progress = [&](bool force) {
    if (!verbose) {
      return;
    }
    const auto now = std::chrono::steady_clock::now();
    if (!force && now < next_progress_at) {
      return;
    }
    const double percent
        = total_pairs == 0
              ? 100.0
              : 100.0 * static_cast<double>(processed_pairs)
                    / static_cast<double>(total_pairs);
    const double elapsed_seconds
        = std::chrono::duration<double>(now - started_at).count();
    std::cerr << std::format(
        "Intersection {} progress: {:.1f}% ({}/{} pairs), non_empty={}, "
        "elapsed={:.1f}s\n",
        label, percent, processed_pairs, total_pairs, non_empty_pairs,
        elapsed_seconds);
    std::cerr.flush();
    next_progress_at = now + std::chrono::seconds(5);
  };
  report_progress(true);
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
        ++non_empty_pairs;
        out_points.push_back(intersection);
        // Store the indices of the prisms that form the intersection
        std::vector<size_t> prism_indices = {idx0, idx1};
        out_indices.push_back(prism_indices);
      }
      ++processed_pairs;
      if (processed_pairs % 10000 == 0) {
        report_progress(false);
      }
    }
  }
  report_progress(true);
};

if (verbose) {
  std::cout << std::format("Prism 31: {}", prisms31.size()) << std::endl;
  std::cout << std::format("Prism 36: {}", prisms36.size()) << std::endl;
  std::cout << std::format("Prism 24: {}", prisms24.size()) << std::endl;
  std::cout << std::format("Prism 27: {}", prisms27.size()) << std::endl;
}

std::vector<std::vector<size_t>> intersection_prism_indices_136;
std::vector<std::vector<hcpwa::Vec<3>>> intersection_points_136;
computend({0, 2, 5}, aabb3d_bounds, prisms31, prisms36, "136", verbose,
          intersection_prism_indices_136, intersection_points_136);

std::vector<std::vector<size_t>> intersection_prism_indices_247;
std::vector<std::vector<hcpwa::Vec<3>>> intersection_points_247;
computend({1, 3, 6}, aabb3d_bounds, prisms24, prisms27, "247", verbose,
          intersection_prism_indices_247, intersection_points_247);

std::vector<std::vector<size_t>> intersection_prism_indices_157;
std::vector<std::vector<hcpwa::Vec<3>>> intersection_points_157;
computend({0, 4, 6}, aabb3d_bounds, prisms51, prisms57, "157", verbose,
          intersection_prism_indices_157, intersection_points_157);

std::vector<std::vector<size_t>> intersection_prism_indices_468;
std::vector<std::vector<hcpwa::Vec<3>>> intersection_points_468;
computend({3, 5, 7}, aabb3d_bounds, prisms84, prisms86, "468", verbose,
          intersection_prism_indices_468, intersection_points_468);

  // ---- Block decomposition (polygon path) ----------------------------------
  //
  // Identical in structure to the triangulated path: the block partition is a
  // property of the plane arrangement, not of how each plane is subdivided.
  // The only difference is block C, which here is a general convex polygon cell
  // rather than a triangle, so its vertex count is only bounded below by 3.
  auto polygon_cells
      = [](const std::vector<hcpwa::PolygonResolution>& polygons) {
          std::vector<std::vector<hcpwa::Vec<2>>> cells;
          cells.reserve(polygons.size());
          for (const auto& resolution : polygons) {
            cells.push_back(resolution.polygon);
          }
          return cells;
        };

  std::array<BlockRegions, 3> blocks_phase0
      = {makePairBlock({0, 2, 5}, {0, 1}, intersection_prism_indices_136,
                       intersection_points_136, "136"),
         makePairBlock({1, 3, 6}, {2, 3}, intersection_prism_indices_247,
                       intersection_points_247, "247"),
         makeSingleBlock({4, 7, -1}, {4, -1}, polygon_cells(polygons58), "58")};
  std::array<BlockRegions, 3> blocks_phase1
      = {makePairBlock({0, 4, 6}, {0, 1}, intersection_prism_indices_157,
                       intersection_points_157, "157"),
         makePairBlock({3, 5, 7}, {2, 3}, intersection_prism_indices_468,
                       intersection_points_468, "468"),
         makeSingleBlock({1, 2, -1}, {4, -1}, polygon_cells(polygons23), "23")};

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

const std::size_t total_phase0_areas
    = intersection_points_136.size() * intersection_prism_indices_247.size()
      * polygons58.size();
std::size_t processed_phase0_areas = 0;
const auto phase0_started_at = std::chrono::steady_clock::now();
auto next_phase0_progress_at = phase0_started_at;
auto report_phase0_progress = [&](bool force) {
  if (!verbose) {
    return;
  }
  const auto now = std::chrono::steady_clock::now();
  if (!force && now < next_phase0_progress_at) {
    return;
  }
  const double percent
      = total_phase0_areas == 0
            ? 100.0
            : 100.0 * static_cast<double>(processed_phase0_areas)
                  / static_cast<double>(total_phase0_areas);
  const double elapsed_seconds
      = std::chrono::duration<double>(now - phase0_started_at).count();
  std::cerr << std::format(
      "Phase0 area assembly progress: {:.1f}% ({}/{} areas), vertices={}, "
      "elapsed={:.1f}s\n",
      percent, processed_phase0_areas, total_phase0_areas,
      intersection_points_phase0.size(), elapsed_seconds);
  std::cerr.flush();
  next_phase0_progress_at = now + std::chrono::seconds(5);
};
report_phase0_progress(true);

for (size_t i136 = 0; i136 < intersection_points_136.size(); i136++) {
  for (size_t i247 = 0; i247 < intersection_prism_indices_247.size();
       i247++) {
    // auto t_start = std::chrono::high_resolution_clock::now();
    for (size_t i58 = 0; i58 < polygons58.size(); i58++) {
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

      // The simplex-id tuple above is always built. Only the explicit 8D
      // vertex list below is optional, and it is the truncated one: its two
      // outer loops are bounded by the 2-element prism-index lists while they
      // index the block vertex lists, which hold at least 3 entries each. The
      // innermost loop over the 2D cell uses the right container, which is the
      // tell. See docs/barycentric_block_reduction_context.md part I.
      if (options.build_8d_regions) {
        intersection_points_phase0.emplace_back();
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
      ++processed_phase0_areas;
      if (processed_phase0_areas % 1000 == 0) {
        report_phase0_progress(false);
      }
    }
    // auto t_end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> t_diff = t_end - t_start;
    // std::cout << "[Timing] i136=" << i136 << ", i247=" << i247 << ": " <<
    // t_diff.count() << "s" << std::endl;
  }
}
report_phase0_progress(true);

const std::size_t total_phase1_areas
    = intersection_points_157.size() * intersection_prism_indices_468.size()
      * polygons23.size();
std::size_t processed_phase1_areas = 0;
const auto phase1_started_at = std::chrono::steady_clock::now();
auto next_phase1_progress_at = phase1_started_at;
auto report_phase1_progress = [&](bool force) {
  if (!verbose) {
    return;
  }
  const auto now = std::chrono::steady_clock::now();
  if (!force && now < next_phase1_progress_at) {
    return;
  }
  const double percent
      = total_phase1_areas == 0
            ? 100.0
            : 100.0 * static_cast<double>(processed_phase1_areas)
                  / static_cast<double>(total_phase1_areas);
  const double elapsed_seconds
      = std::chrono::duration<double>(now - phase1_started_at).count();
  std::cerr << std::format(
      "Phase1 area assembly progress: {:.1f}% ({}/{} areas), vertices={}, "
      "elapsed={:.1f}s\n",
      percent, processed_phase1_areas, total_phase1_areas,
      intersection_points_phase1.size(), elapsed_seconds);
  std::cerr.flush();
  next_phase1_progress_at = now + std::chrono::seconds(5);
};
report_phase1_progress(true);

for (size_t i157 = 0; i157 < intersection_points_157.size(); i157++) {
  for (size_t i468 = 0; i468 < intersection_prism_indices_468.size();
       i468++) {
    for (size_t i23 = 0; i23 < polygons23.size(); i23++) {
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

      // See the phase-0 loop above.
      if (options.build_8d_regions) {
        intersection_points_phase1.emplace_back();
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
      ++processed_phase1_areas;
      if (processed_phase1_areas % 1000 == 0) {
        report_phase1_progress(false);
      }
    }
  }
}
report_phase1_progress(true);

checkProductStructure(blocks_phase0, intersection_prism_indices_phase0.size(),
                      0);
checkProductStructure(blocks_phase1, intersection_prism_indices_phase1.size(),
                      1);

PhaseIntersectionResult result;
result.intersection_prism_indices_phase0
    = std::move(intersection_prism_indices_phase0);
result.intersection_points_phase0 = std::move(intersection_points_phase0);
result.intersection_prism_indices_phase1
    = std::move(intersection_prism_indices_phase1);
result.intersection_points_phase1 = std::move(intersection_points_phase1);
result.blocks_phase0 = std::move(blocks_phase0);
result.blocks_phase1 = std::move(blocks_phase1);
// area_bounds_phase{0,1} are left empty on this path: it has no common
// refinement (compute_polygon_areas_vertices does not build one, and the global
// approximator's border conditions use the box corners instead), so composing
// 1.6M eight-dimensional boxes here would be pure waste.
return result;
}

TriangleAreasVerticesResult compute_triangle_areas_vertices(
    double N, double F, double v, double w, double b51, double b57, double b84,
    double b86, double b31, double b36, double b24, double b27, double f2min,
    double f3min, double f5min, double f8min, double f2max, double f3max,
    double f5max, double f8max, TriangleAreasOptions options) {
  const bool verbose = options.verbose;
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
      static_cast<hcpwa::Float>(N), options);
  auto& intersection_points_phase0
      = intersection_result.intersection_points_phase0;
  auto& intersection_prism_indices_phase0
      = intersection_result.intersection_prism_indices_phase0;
  auto& intersection_points_phase1
      = intersection_result.intersection_points_phase1;
  auto& intersection_prism_indices_phase1
      = intersection_result.intersection_prism_indices_phase1;

  // The border-condition LP needs vertices of the common refinement between the
  // two phase partitions, not just vertices of each phase partition separately.
  // Reuse the exact phase-area prism tuples above so the common-refinement cells
  // inherit the same [31,36,24,27,58] and [51,57,84,86,23] indexing contracts.
  CommonRefinementResult common_refinement
      = compute_common_refinement_area_vertices(
          prisms31, prisms36, prisms24, prisms27, prisms58, prisms51, prisms57,
          prisms84, prisms86, prisms23, intersection_prism_indices_phase0,
          intersection_prism_indices_phase1,
          intersection_result.area_bounds_phase0,
          intersection_result.area_bounds_phase1,
          static_cast<hcpwa::Float>(N), verbose);

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
  result.common_refinement = std::move(common_refinement);
  result.blocks_phase0 = std::move(intersection_result.blocks_phase0);
  result.blocks_phase1 = std::move(intersection_result.blocks_phase1);

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
    std::cout << "result.common_refinement.areas.size(): "
              << result.common_refinement.areas.size() << '\n';
  }

  return result;
}

PolygonAreasVerticesResult compute_polygon_areas_vertices(
    double N, double F, double v, double w, double b51, double b57, double b84,
    double b86, double b31, double b36, double b24, double b27, double f2min,
    double f3min, double f5min, double f8min, double f2max, double f3max,
    double f5max, double f8max, TriangleAreasOptions options) {
  const bool verbose = options.verbose;
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
      static_cast<hcpwa::Float>(N), options);
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
  result.blocks_phase0 = std::move(intersection_result.blocks_phase0);
  result.blocks_phase1 = std::move(intersection_result.blocks_phase1);

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
