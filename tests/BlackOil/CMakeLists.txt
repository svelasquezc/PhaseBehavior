# Tests need to be added as executables first
add_executable(TestBlackOil TestBlackOilEstimation.cpp)

# C++17 syntax
target_compile_features(TestBlackOil PRIVATE cxx_std_17)

# Should be linked to the main library, as well as the Catch2 testing library
target_link_libraries(TestBlackOil PRIVATE PhaseBehavior Catch2::Catch2WithMain)

# If you register a test, then ctest and make test will run it.
# You can also run examples and check the output, as well.
add_test(NAME TestBlackOil COMMAND TestBlackOil) # Command can be a target

# Copy the PVT File to the same folder as the executable file.
add_custom_command(TARGET TestBlackOil POST_BUILD    # Adds a post-build event to MyTest
COMMAND ${CMAKE_COMMAND} -E copy_if_different   # which executes "cmake - E copy_if_different..."
"${CMAKE_CURRENT_SOURCE_DIR}/BOMHandoutData.csv" # <--this is also in-file
$<TARGET_FILE_DIR:TestBlackOil>)                 # <--this is out-file path