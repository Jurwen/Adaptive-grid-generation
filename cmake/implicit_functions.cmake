if (TARGET implicit_functions::implicit_functions)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    implicit_functions
    GIT_REPOSITORY https://github.com/duxingyi-charles/implicit_functions.git
    #GIT_REPOSITORY https://github.com/Jurwen/implicit_functions.git
    GIT_TAG main
    )

FetchContent_MakeAvailable(implicit_functions)
