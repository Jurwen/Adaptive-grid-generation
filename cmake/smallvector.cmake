if (TARGET smallvector::smallvector)
    return()
endif()

message(STATUS "Third-party (external): creating target 'smallvector::smallvector'")

include(CPM)
CPMAddPackage(
    NAME smallvector
    GITHUB_REPOSITORY thelink2012/SmallVector
    GIT_TAG febc8cb7b1a83d902b86dd1612feb7c86c690186
    DOWNLOAD_ONLY YES
)

add_library(smallvector INTERFACE)
target_include_directories(smallvector INTERFACE ${smallvector_SOURCE_DIR})
set_target_properties(smallvector PROPERTIES SYSTEM ON)
set_target_properties(smallvector PROPERTIES FOLDER third_party)
add_library(smallvector::smallvector ALIAS smallvector)
