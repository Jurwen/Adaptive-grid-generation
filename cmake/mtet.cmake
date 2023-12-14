if (TARGET mtet::mtet)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    mtet
    GIT_REPOSITORY https://github.com/qnzhou/mtet.git
    GIT_TAG main
    )

FetchContent_MakeAvailable(mtet)
