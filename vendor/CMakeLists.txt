FetchContent_Declare(
  linalg
  GIT_REPOSITORY https://github.com/sgorsten/linalg.git
  GIT_TAG main
)

FetchContent_MakeAvailable(linalg)

add_library(linalg INTERFACE)
target_include_directories(linalg INTERFACE ${linalg_SOURCE_DIR})

FetchContent_Declare(
  cddlib
  URL https://github.com/cddlib/cddlib/archive/refs/heads/master.zip
)
FetchContent_Populate(cddlib)
file(GLOB_RECURSE CDD_SRC_FILES CONFIGURE_DEPENDS
"${cddlib_SOURCE_DIR}/lib-src/*.c"
)
add_library(cddlib STATIC ${CDD_SRC_FILES})
target_include_directories(cddlib PUBLIC ${cddlib_SOURCE_DIR}/lib-src)
set_target_properties(cddlib PROPERTIES LINKER_LANGUAGE C)

FetchContent_Declare(
  function2
  GIT_REPOSITORY https://github.com/Naios/function2.git
  GIT_TAG master  # or a specific commit/tag
)

FetchContent_MakeAvailable(function2)

# HIGHS
FetchContent_Declare(
  highs
  GIT_REPOSITORY https://github.com/ERGO-Code/HiGHS.git
  GIT_TAG master
)
FetchContent_MakeAvailable(highs)

# EIGEN
# Disable all testing globally (also helps other deps)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)

# Eigen-specific knobs (names vary by Eigen version, but these are commonly honored)
set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)

FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 5.0.1
)
FetchContent_MakeAvailable(eigen)
# FetchContent_Populate(eigen)

# add_library(project_eigen INTERFACE)
# target_include_directories(project_eigen INTERFACE ${eigen_SOURCE_DIR})

# spdlog
set(SPDLOG_USE_STD_FORMAT on)
FetchContent_Declare(
  spdlog
  GIT_REPOSITORY https://github.com/gabime/spdlog.git
  GIT_TAG v1.14.0
)
FetchContent_MakeAvailable(spdlog)

# google or-tools
# include(FetchContent)
# FetchContent_Declare(
#   ortools
#   GIT_REPOSITORY https://github.com/google/or-tools.git
#   GIT_TAG        stable # or a specific version like v9.8
# )
# FetchContent_MakeAvailable(ortools)

# add_library(spdlog INTERFACE)
# target_include_directories(spdlog INTERFACE ${spdlog_SOURCE_DIR}/include)