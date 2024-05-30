# Tests need to be added as executables first
add_executable(PhaseEnvelope PhaseEnvelope.cpp)

# Tests need to be added as executables first
add_executable(BlackOilEstimation BlackOilEstimation.cpp)

# Tests need to be added as executables first
add_executable(OffendingComposition OffendingComposition.cpp)

# C++17 syntax
target_compile_features(PhaseEnvelope PRIVATE cxx_std_17)
target_compile_features(BlackOilEstimation PRIVATE cxx_std_17)
target_compile_features(OffendingComposition PRIVATE cxx_std_17)

# Should be linked to the main library, as well as the Catch2 testing library
target_link_libraries(PhaseEnvelope PRIVATE PhaseBehavior)
target_link_libraries(BlackOilEstimation PRIVATE PhaseBehavior)
target_link_libraries(OffendingComposition PRIVATE PhaseBehavior)

# IDEs should put the headers in a nice place
source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  PREFIX "Header Files"
  FILES ${HEADERS}
  )