cmake_minimum_required(VERSION 3.14.0)
project(
    PhaseBehavior
    VERSION 0.1.0
    LANGUAGES C CXX
    DESCRIPTION "Multicomponent Liquid/Gas Phase Behavior Model"
    )


# FetchContent added in CMake 3.11, downloads during the configure step
# FetchContent_MakeAvailable was added in CMake 3.14; simpler usage
include(FetchContent)

find_package(Eigen3 NO_MODULE)
if (NOT TARGET Eigen3::Eigen)

    
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
    option(BUILD_TESTING OFF)
    set(EIGEN_BUILD_PKGCONFIG OFF)
    set( OFF)
    FetchContent_MakeAvailable(Eigen)

endif (NOT TARGET Eigen3::Eigen)

# Only do these if this is the main project, and not if it is included through add_subdirectory
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)

    # Let's ensure -std=c++xx instead of -std=g++xx
    set(CMAKE_CXX_EXTENSIONS OFF)

    # Let's nicely support folders in IDEs
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)

    # Docs only available if this is the main app
    find_package(Doxygen)
    if(Doxygen_FOUND)
        add_subdirectory(docs)
    else()
        message(STATUS "Doxygen not found, not building docs")
    endif()

    # Testing only available if this is the main app
    # Note this needs to be done in the main CMakeLists
    # since it calls enable_testing, which must be in the
    # main CMakeLists.
    include(CTest)

endif()

add_subdirectory(src)

if (BUILD_TESTING)
    add_subdirectory(tests)
    add_subdirectory(Data)
endif()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
