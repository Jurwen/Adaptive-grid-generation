cmake_minimum_required(VERSION 3.14)

project(adaptive_mesh_refinement LANGUAGES CXX C)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(cli11)
include(nlohmann-json)
include(mtet)
include(implicit_functions)
include(unordered_dense)

add_executable(isosurfacing "app/isosurfacing.cpp")
target_compile_features(isosurfacing PRIVATE cxx_std_20)
target_include_directories(isosurfacing PUBLIC ${CMAKE_CURRENT_LIST_DIR}/include)
target_link_libraries(isosurfacing PRIVATE
    CLI11::CLI11 mtet::mtet implicit_functions::implicit_functions
    unordered_dense::unordered_dense)