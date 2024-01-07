cmake_minimum_required(VERSION 3.14)

project(adaptive_mesh_refinement LANGUAGES CXX C)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(cli11)
include(nlohmann-json)
include(mtet)
include(implicit_functions)
include(convex_hull_membership)
include(smallvector)
include(unordered_dense)

# Visual Studio compiler with static runtime libraries
if(MSVC AND MT)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /MTd")
endif()

find_package(GUROBI REQUIRED)

include_directories(${GUROBI_INCLUDE_DIRS})

add_executable(isosurfacing "app/isosurfacing.cpp")

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  include(FeatureSummary)
  feature_summary(WHAT ALL)
endif()

target_compile_features(isosurfacing PRIVATE cxx_std_20)
target_include_directories(isosurfacing PUBLIC ${CMAKE_CURRENT_LIST_DIR}/include)
target_link_libraries(isosurfacing PRIVATE
    CLI11::CLI11 
    mtet::mtet 
    implicit_functions::implicit_functions
    convex_hull_membership
    unordered_dense::unordered_dense
    smallvector::smallvector
    ${GUROBI_LIBRARY})

if(CXX)
  set(CMAKE_CXX_STANDARD 11)
  target_link_libraries(isosurfacing optimized ${GUROBI_CXX_LIBRARY}
                        debug ${GUROBI_CXX_DEBUG_LIBRARY})
endif()
