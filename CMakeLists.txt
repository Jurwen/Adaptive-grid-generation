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
    implicit_functions::implicit_functions
    nlohmann_json::nlohmann_json
    smallvector::smallvector
    unordered_dense::unordered_dense)

add_executable(gridgen "app/gridgen")

target_compile_features(gridgen PRIVATE cxx_std_20)
target_include_directories(gridgen PUBLIC ${CMAKE_CURRENT_LIST_DIR}/src)
target_link_libraries(gridgen PRIVATE
    adaptive_mesh_refinement
    CLI11::CLI11
    mtet::mtet 
    convex_hull_membership
    unordered_dense::unordered_dense
    smallvector::smallvector
    )


option(GRID_GEN_TESTS "Build adaptive grid gen unit tests" ON)
if (GRID_GEN_TESTS)
 include(CTest)
 enable_testing()
 include(Catch2)

 file(GLOB TEST_FILES "${CMAKE_CURRENT_LIST_DIR}/tests/*.cpp")
 add_executable(grid_gen_tests ${TEST_FILES})
 target_link_libraries(grid_gen_tests
    PRIVATE 
    adaptive_mesh_refinement
    Catch2::Catch2
    mtet::mtet
    unordered_dense::unordered_dense
    smallvector::smallvector)
 target_compile_features(grid_gen_tests PRIVATE cxx_std_20)
 target_compile_definitions(grid_gen_tests PRIVATE CATCH_CONFIG_ENABLE_BENCHMARKING)
target_compile_definitions(
    grid_gen_tests
    PRIVATE                   
    -DTEST_FILE="${CMAKE_CURRENT_LIST_DIR}/data"
)
 catch_discover_tests(grid_gen_tests)
endif()
