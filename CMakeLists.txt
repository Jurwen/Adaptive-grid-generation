cmake_minimum_required(VERSION 3.14)

project(adaptive_mesh_refinement LANGUAGES CXX C)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(cli11)
include(Eigen3)
include(nlohmann-json)
include(mtet)
include(implicit_functions)
include(convex_hull_membership)
include(smallvector)
include(unordered_dense)

# project library
file(GLOB SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/src/*.cpp ${CMAKE_CURRENT_LIST_DIR}/src/*.h)
add_library(adaptive_mesh_refinement STATIC ${SRC_FILES})
target_include_directories(adaptive_mesh_refinement PUBLIC
     ${CMAKE_CURRENT_LIST_DIR}/src)
target_link_libraries(adaptive_mesh_refinement
    PUBLIC 
    Eigen3::Eigen
    mtet::mtet 
    convex_hull_membership
    nlohmann_json::nlohmann_json
    smallvector::smallvector)

add_executable(isosurfacing "app/isosurfacing.cpp")

#if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
#  src(FeatureSummary)
#  feature_summary(WHAT ALL)
#endif()

target_compile_features(isosurfacing PRIVATE cxx_std_20)
target_include_directories(isosurfacing PUBLIC ${CMAKE_CURRENT_LIST_DIR}/src)
target_link_libraries(isosurfacing PRIVATE
    adaptive_mesh_refinement
    CLI11::CLI11
    mtet::mtet 
    implicit_functions::implicit_functions
    convex_hull_membership
    unordered_dense::unordered_dense
    smallvector::smallvector
    )