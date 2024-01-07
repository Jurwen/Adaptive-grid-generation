if (TARGET convex_hull_membership)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    convex_hull_membership
    GIT_REPOSITORY https://github.com/qnzhou/convex_hull_membership.git
    GIT_TAG main
    )

FetchContent_MakeAvailable(convex_hull_membership)
