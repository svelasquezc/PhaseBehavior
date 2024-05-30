if (BUILD_EXAMPLES)
    add_custom_target(CopyDataMain ALL DEPENDS PhaseEnvelope)
    # Copy the PVT File to the same folder as the executable file.
    add_custom_command(TARGET CopyDataMain POST_BUILD    # Adds a post-build event to MyTest
    COMMAND ${CMAKE_COMMAND} -E copy_if_different   # which executes "cmake - E copy_if_different..."
    "${CMAKE_CURRENT_SOURCE_DIR}/PVT.csv"           # <--this is also in-file
    "${CMAKE_CURRENT_SOURCE_DIR}/PVTPR.csv"           # <--this is also in-file
    "${CMAKE_CURRENT_SOURCE_DIR}/PVTSRK.csv"           # <--this is also in-file
    "${CMAKE_CURRENT_SOURCE_DIR}/InteractionCoefficients.csv"           # <--this is also in-file
    $<TARGET_FILE_DIR:PhaseEnvelope>)                 # <--this is out-file path
endif()


if(ENABLE_TESTING)

    add_custom_target(CopyDataTest ALL DEPENDS TestMixture)
    # Copy the PVT File to the same folder as the executable file.
    add_custom_command(TARGET CopyDataTest POST_BUILD    # Adds a post-build event to MyTest
    COMMAND ${CMAKE_COMMAND} -E copy_if_different   # which executes "cmake - E copy_if_different..."
    "${CMAKE_CURRENT_SOURCE_DIR}/PVT.csv"           # <--this is also in-file
    "${CMAKE_CURRENT_SOURCE_DIR}/PVTPR.csv"           # <--this is also in-file
    "${CMAKE_CURRENT_SOURCE_DIR}/PVTSRK.csv"           # <--this is also in-file
    "${CMAKE_CURRENT_SOURCE_DIR}/InteractionCoefficients.csv"           # <--this is also in-file
    $<TARGET_FILE_DIR:TestMixture>)                 # <--this is out-file path

    add_custom_target(CopyDataBOM ALL DEPENDS TestBlackOil)
    # Copy the PVT File to the same folder as the executable file.
    add_custom_command(TARGET CopyDataBOM POST_BUILD    # Adds a post-build event to MyTest
    COMMAND ${CMAKE_COMMAND} -E copy_if_different   # which executes "cmake - E copy_if_different..."
    "${CMAKE_CURRENT_SOURCE_DIR}/PVTPR.csv"           # <--this is also in-file
    "${CMAKE_CURRENT_SOURCE_DIR}/InteractionCoefficients.csv"           # <--this is also in-file
    $<TARGET_FILE_DIR:TestBlackOil>)                 # <--this is out-file path

endif()