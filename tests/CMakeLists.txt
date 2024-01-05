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
add_executable(TestRandomMixingRule TestRandomMixingRule.cpp)

# C++17 syntax
target_compile_features(TestRandomMixingRule PRIVATE cxx_std_17)

# Should be linked to the main library, as well as the Catch2 testing library
target_link_libraries(TestRandomMixingRule PRIVATE PhaseBehavior Catch2::Catch2WithMain)

# If you register a test, then ctest and make test will run it.
# You can also run examples and check the output, as well.
add_test(NAME TestRandomMixingRule COMMAND TestRandomMixingRule) # Command can be a target

# Tests need to be added as executables first
add_executable(TestPR TestCreatePengRobinson.cpp)

# C++17 syntax
target_compile_features(TestPR PRIVATE cxx_std_17)

# Should be linked to the main library, as well as the Catch2 testing library
target_link_libraries(TestPR PRIVATE PhaseBehavior Catch2::Catch2WithMain)

# If you register a test, then ctest and make test will run it.
# You can also run examples and check the output, as well.
add_test(NAME TestPR COMMAND TestPR) # Command can be a target

# Tests need to be added as executables first
add_executable(TestSRK TestCreateSoaveRedlichKwong.cpp)

# C++17 syntax
target_compile_features(TestSRK PRIVATE cxx_std_17)

# Should be linked to the main library, as well as the Catch2 testing library
target_link_libraries(TestSRK PRIVATE PhaseBehavior Catch2::Catch2WithMain)

# If you register a test, then ctest and make test will run it.
# You can also run examples and check the output, as well.
add_test(NAME TestSRK COMMAND TestSRK) # Command can be a target

# Tests need to be added as executables first
add_executable(TestNR TestNewtonRaphson.cpp)

# C++17 syntax
target_compile_features(TestNR PRIVATE cxx_std_17)

# Should be linked to the main library, as well as the Catch2 testing library
target_link_libraries(TestNR PRIVATE PhaseBehavior Catch2::Catch2WithMain)

# If you register a test, then ctest and make test will run it.
# You can also run examples and check the output, as well.
add_test(NAME TestNR COMMAND TestNR) # Command can be a target

# Tests need to be added as executables first
add_executable(TestBisection TestBisection.cpp)

# C++17 syntax
target_compile_features(TestBisection PRIVATE cxx_std_17)

# Should be linked to the main library, as well as the Catch2 testing library
target_link_libraries(TestBisection PRIVATE PhaseBehavior Catch2::Catch2WithMain)

# If you register a test, then ctest and make test will run it.
# You can also run examples and check the output, as well.
add_test(NAME TestBisection COMMAND TestBisection) # Command can be a target

# Tests need to be added as executables first
add_executable(TestVLE TestRachfordRice.cpp)

# C++17 syntax
target_compile_features(TestVLE PRIVATE cxx_std_17)

# Should be linked to the main library, as well as the Catch2 testing library
target_link_libraries(TestVLE PRIVATE PhaseBehavior Catch2::Catch2WithMain)

# If you register a test, then ctest and make test will run it.
# You can also run examples and check the output, as well.
add_test(NAME TestVLE COMMAND TestVLE) # Command can be a target

# Tests need to be added as executables first
add_executable(TestSSM TestSSM.cpp)

# C++17 syntax
target_compile_features(TestSSM PRIVATE cxx_std_17)

# Should be linked to the main library, as well as the Catch2 testing library
target_link_libraries(TestSSM PRIVATE PhaseBehavior Catch2::Catch2WithMain)

# If you register a test, then ctest and make test will run it.
# You can also run examples and check the output, as well.
add_test(NAME TestSSM COMMAND TestSSM) # Command can be a target