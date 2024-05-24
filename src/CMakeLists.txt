file(GLOB HEADERS CMAKE_CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/include/*.hpp")

add_library(PhaseBehavior INTERFACE)

# We need this directory, and users of our library will need it too
target_include_directories(PhaseBehavior INTERFACE ../include)

# All users of this library will need at least C++17
target_compile_features(PhaseBehavior INTERFACE cxx_std_17)

#
target_link_libraries(PhaseBehavior INTERFACE Eigen3::Eigen)

# Tests need to be added as executables first
add_executable(PhaseEnvelope PhaseEnvelope/main.cpp)

# Tests need to be added as executables first
add_executable(BlackOilEstimation BlackOilEstimation/main.cpp)

# C++17 syntax
target_compile_features(PhaseEnvelope PRIVATE cxx_std_17)
target_compile_features(BlackOilEstimation PRIVATE cxx_std_17)

# Should be linked to the main library, as well as the Catch2 testing library
target_link_libraries(PhaseEnvelope PRIVATE PhaseBehavior)
target_link_libraries(BlackOilEstimation PRIVATE PhaseBehavior)

# IDEs should put the headers in a nice place
source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  PREFIX "Header Files"
  FILES ${HEADERS}
  )