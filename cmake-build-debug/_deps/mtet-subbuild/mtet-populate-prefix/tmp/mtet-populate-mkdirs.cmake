# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/Users/yiwenju/Dev/adaptive-mesh-refinement/cmake-build-debug/_deps/mtet-src"
  "/Users/yiwenju/Dev/adaptive-mesh-refinement/cmake-build-debug/_deps/mtet-build"
  "/Users/yiwenju/Dev/adaptive-mesh-refinement/cmake-build-debug/_deps/mtet-subbuild/mtet-populate-prefix"
  "/Users/yiwenju/Dev/adaptive-mesh-refinement/cmake-build-debug/_deps/mtet-subbuild/mtet-populate-prefix/tmp"
  "/Users/yiwenju/Dev/adaptive-mesh-refinement/cmake-build-debug/_deps/mtet-subbuild/mtet-populate-prefix/src/mtet-populate-stamp"
  "/Users/yiwenju/Dev/adaptive-mesh-refinement/cmake-build-debug/_deps/mtet-subbuild/mtet-populate-prefix/src"
  "/Users/yiwenju/Dev/adaptive-mesh-refinement/cmake-build-debug/_deps/mtet-subbuild/mtet-populate-prefix/src/mtet-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/yiwenju/Dev/adaptive-mesh-refinement/cmake-build-debug/_deps/mtet-subbuild/mtet-populate-prefix/src/mtet-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/yiwenju/Dev/adaptive-mesh-refinement/cmake-build-debug/_deps/mtet-subbuild/mtet-populate-prefix/src/mtet-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
