cmake_minimum_required(VERSION 3.27)

project(
    hcpwa
    VERSION 1.0
    LANGUAGES C CXX
)

include(FetchContent)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_USE_RESPONSE_FILE_FOR_INCLUDES OFF)

add_subdirectory(vendor)
add_subdirectory(${PROJECT_NAME})
add_subdirectory(solver)

set(TESTS_OPTION "${PROJECT_NAME}_BUILD_TESTS")
option(${TESTS_OPTION} "build tests" ON)

if (${TESTS_OPTION})
    add_subdirectory(test)
endif()
