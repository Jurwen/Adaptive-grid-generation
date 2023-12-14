cmake_minimum_required(VERSION 3.14)

project(adaptive_mesh_refinement LANGUAGES CXX C)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(cli11)
include(nlohmann-json)
include(mtet)
include(implicit_functions)

# project library
file(GLOB SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/src/*.cpp ${CMAKE_CURRENT_LIST_DIR}/src/*.h)

add_library(amr STATIC ${SRC_FILES})
target_include_directories(amr PUBLIC ${CMAKE_CURRENT_LIST_DIR}/src)
target_compile_features(amr PRIVATE cxx_std_17)
set_target_properties(amr PROPERTIES LINKER_LANGUAGE CXX)

add_executable(isosurfacing "app/isosurfacing.cpp")
target_compile_features(isosurfacing PRIVATE cxx_std_17)
target_link_libraries(isosurfacing PRIVATE amr
    CLI11::CLI11 mtet::mtet implicit_functions::implicit_functions)