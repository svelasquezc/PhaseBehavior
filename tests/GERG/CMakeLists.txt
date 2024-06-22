# Tests need to be added as executables first
add_executable(TestBinaryReducingParameters TestBinaryReducingParameters.cpp)

# C++17 syntax
target_compile_features(TestBinaryReducingParameters PRIVATE cxx_std_17)

# Should be linked to the main library, as well as the Catch2 testing library
target_link_libraries(TestBinaryReducingParameters PRIVATE PhaseBehavior Catch2::Catch2WithMain)

# If you register a test, then ctest and make test will run it.
# You can also run examples and check the output, as well.
add_test(NAME TestBinaryReducingParameters COMMAND TestBinaryReducingParameters) # Command can be a target
