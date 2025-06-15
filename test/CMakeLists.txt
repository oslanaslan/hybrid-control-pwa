FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

include(GoogleTest)

file(GLOB SUBDIRS RELATIVE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/*)

foreach(TEST_DIR ${SUBDIRS})
  if (IS_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/${TEST_DIR})
    set(TEST_NAME ${TEST_DIR})
    set(TEST_DIR ${CMAKE_CURRENT_LIST_DIR}/${TEST_NAME})
    set(TEST_TARGET "test_${TEST_NAME}")
    if (EXISTS "${TEST_DIR}/CmakeLists.txt")
      message("Test subdirectory: ${TEST_NAME}")
      add_subdirectory(${TEST_NAME})
    else()
      file(GLOB_RECURSE TEST_SRC_FILES CONFIGURE_DEPENDS
      "./**/*.cpp"
      )
      add_executable(${TEST_TARGET} ${TEST_SRC_FILES})
      message("Test target: ${TEST_NAME}")
      target_compile_features(${TEST_TARGET} PUBLIC cxx_std_20)
      set_target_properties(${TEST_TARGET} PROPERTIES CXX_EXTENSIONS OFF)
      
      target_include_directories(${TEST_TARGET} PUBLIC ${CMAKE_CURRENT_LIST_DIR})
      
      target_link_libraries(${TEST_TARGET} PUBLIC ${PROJECT_NAME} GTest::gtest_main)
      
      gtest_discover_tests(${TEST_TARGET})
    endif()
  endif()
endforeach()