add_custom_target(CopyData ALL DEPENDS TestMixture)

# Copy the PVT File to the same folder as the executable file.
add_custom_command(TARGET CopyData POST_BUILD    # Adds a post-build event to MyTest
COMMAND ${CMAKE_COMMAND} -E copy_if_different   # which executes "cmake - E copy_if_different..."
"${CMAKE_CURRENT_SOURCE_DIR}/PVT.csv"           # <--this is also in-file
"${CMAKE_CURRENT_SOURCE_DIR}/InteractionCoefficients.csv"           # <--this is also in-file
$<TARGET_FILE_DIR:TestMixture>)                 # <--this is out-file path