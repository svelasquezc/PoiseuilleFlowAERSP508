cmake_minimum_required(VERSION 3.14.0)
project(PoiseuilleFlow
 VERSION 0.1.0
 LANGUAGES C CXX
 DESCRIPTION "AERSP 508 Project 1: Poiseuille Flow 1D"
 )

# Only do these if this is the main project, and not if it is included through add_subdirectory
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)

    # Let's ensure -std=c++xx instead of -std=g++xx
    set(CMAKE_CXX_EXTENSIONS OFF)

    # Let's nicely support folders in IDEs
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)


    find_package(Eigen3 NO_MODULE)
    if (NOT TARGET Eigen3::Eigen)

        include(FetchContent)
        FetchContent_Declare(
        Eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0
        GIT_SHALLOW TRUE
        GIT_PROGRESS TRUE)
        set(EIGEN_BUILD_DOC OFF)
        # note: To disable eigen tests,
        # you should put this code in a add_subdirectory to avoid to change
        # BUILD_TESTING for your own project too since variables are directory
        # scoped
        set(BUILD_TESTING OFF)
        set(EIGEN_BUILD_PKGCONFIG OFF)
        set( OFF)
        FetchContent_MakeAvailable(Eigen)

    endif (NOT TARGET Eigen3::Eigen)


    # Testing only available if this is the main app
    # Note this needs to be done in the main CMakeLists
    # since it calls enable_testing, which must be in the
    # main CMakeLists.
    include(CTest)

endif(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)

add_subdirectory(src)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
