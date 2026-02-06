file(GLOB_RECURSE SRC_FILES CONFIGURE_DEPENDS
"src/**.cpp"
)

add_library(${PROJECT_NAME} STATIC ${SRC_FILES})

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_23)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_LIST_DIR}/include)

target_link_libraries(${PROJECT_NAME} PUBLIC linalg cddlib function2 highs eigen spdlog)