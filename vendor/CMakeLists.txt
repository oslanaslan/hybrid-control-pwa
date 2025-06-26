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