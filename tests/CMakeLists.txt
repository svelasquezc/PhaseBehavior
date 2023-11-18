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
add_test(NAME testcubicrootstest COMMAND testcubicroots) # Command can be a target

add_subdirectory(TestMixture)