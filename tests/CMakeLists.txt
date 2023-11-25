# Testing library
FetchContent_Declare(
  catch
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.4.0)
FetchContent_MakeAvailable(catch)
# Adds Catch2::Catch2

# Tests need to be added as executables first
add_executable(testcubicroots TestCubicRootsFind.cpp)

# C++17 syntax
target_compile_features(testcubicroots PRIVATE cxx_std_17)

# Should be linked to the main library, as well as the Catch2 testing library
target_link_libraries(testcubicroots PRIVATE PhaseBehavior Catch2::Catch2WithMain)

# If you register a test, then ctest and make test will run it.
# You can also run examples and check the output, as well.
add_test(NAME testcubicroots COMMAND testcubicroots) # Command can be a target

# Tests need to be added as executables first
add_executable(TestMixture TestCreateMixture.cpp)

# C++17 syntax
target_compile_features(TestMixture PRIVATE cxx_std_17)

# Should be linked to the main library, as well as the Catch2 testing library
target_link_libraries(TestMixture PRIVATE PhaseBehavior Catch2::Catch2WithMain)

# If you register a test, then ctest and make test will run it.
# You can also run examples and check the output, as well.
add_test(NAME TestMixtureTest COMMAND TestMixture) # Command can be a target

# Tests need to be added as executables first
add_executable(TestNonRandomMixingRule TestNonRandomMixingRule.cpp)

# C++17 syntax
target_compile_features(TestNonRandomMixingRule PRIVATE cxx_std_17)

# Should be linked to the main library, as well as the Catch2 testing library
target_link_libraries(TestNonRandomMixingRule PRIVATE PhaseBehavior Catch2::Catch2WithMain)

# If you register a test, then ctest and make test will run it.
# You can also run examples and check the output, as well.
add_test(NAME TestNonRandomMixingRule COMMAND TestNonRandomMixingRule) # Command can be a target

# Tests need to be added as executables first
add_executable(TestPR TestCreatePengRobinson.cpp)

# C++17 syntax
target_compile_features(TestPR PRIVATE cxx_std_17)

# Should be linked to the main library, as well as the Catch2 testing library
target_link_libraries(TestPR PRIVATE PhaseBehavior Catch2::Catch2WithMain)

# If you register a test, then ctest and make test will run it.
# You can also run examples and check the output, as well.
add_test(NAME TestPR COMMAND TestPR) # Command can be a target