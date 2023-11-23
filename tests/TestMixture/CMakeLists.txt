# Tests need to be added as executables first
add_executable(TestMixture TestCreateMixture.cpp)

# C++17 syntax
target_compile_features(TestMixture PRIVATE cxx_std_17)

# Should be linked to the main library, as well as the Catch2 testing library
target_link_libraries(TestMixture PRIVATE PhaseBehavior Catch2::Catch2WithMain)

# If you register a test, then ctest and make test will run it.
# You can also run examples and check the output, as well.
add_test(NAME TestMixtureTest COMMAND TestMixture) # Command can be a target

# Copy the PVT File to the same folder as the executable file.
add_custom_command(TARGET TestMixture POST_BUILD    # Adds a post-build event to MyTest
    COMMAND ${CMAKE_COMMAND} -E copy_if_different   # which executes "cmake - E copy_if_different..."
    "${CMAKE_CURRENT_SOURCE_DIR}/PVT.csv"           # <--this is also in-file
    "${CMAKE_CURRENT_SOURCE_DIR}/InteractionCoefficients.csv"           # <--this is also in-file
    $<TARGET_FILE_DIR:TestMixture>)                 # <--this is out-file path