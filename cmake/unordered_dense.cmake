if(TARGET unordered_dense::unordered_dense)
    return()
endif()

message(STATUS "Third-party (external): creating target 'unordered_dense::unordered_dense'")

include(CPM)
CPMAddPackage(
  NAME unordered_dense
  GITHUB_REPOSITORY martinus/unordered_dense
  GIT_TAG v4.1.2
)

set_target_properties(unordered_dense PROPERTIES SYSTEM ON)
set_target_properties(unordered_dense PROPERTIES FOLDER third_party)
