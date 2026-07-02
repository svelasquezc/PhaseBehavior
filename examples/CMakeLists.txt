function(create_phase_behavior_exe executable_name sources)
  add_executable(${executable_name} ${sources})
  # C++17 syntax
  target_compile_features(${executable_name} PRIVATE cxx_std_17)
  # Should be linked to the main library
  target_link_libraries(${executable_name} PRIVATE PhaseBehavior)
endfunction(create_phase_behavior_exe)

create_phase_behavior_exe(PhaseEnvelope PhaseEnvelope.cpp)
create_phase_behavior_exe(PhaseEnvelopeSPE3SI PhaseEnvelopeSPE3SI.cpp)
create_phase_behavior_exe(BlackOilEstimation BlackOilEstimation.cpp)
create_phase_behavior_exe(OffendingComposition OffendingComposition.cpp)
create_phase_behavior_exe(DensityValidation DensityValidation.cpp)
create_phase_behavior_exe(C1C3C8twoPhaseFlash C1C3C8twoPhaseFlash.cpp)
create_phase_behavior_exe(C1C3C8ternaryDiagram C1C3C8ternaryDiagram.cpp)
create_phase_behavior_exe(C1toC6PhaseEnvelope C1toC6PhaseEnvelope.cpp)
create_phase_behavior_exe(C1toC6twoPhaseFlash C1toC6twoPhaseFlash.cpp)
create_phase_behavior_exe(C1toC6LiqFracPlot C1toC6LiqFracPlot.cpp)
create_phase_behavior_exe(C1toC6SPE3LiqFracPlot C1toC6-SPE3-LiqFracPlot.cpp)
create_phase_behavior_exe(C1C7pTwoPhaseFlash C1C7pTwoPhaseFlash.cpp)
create_phase_behavior_exe(SPE3twoPhaseFlash SPE3twoPhaseFlash.cpp)
create_phase_behavior_exe(C1C3C8twoPhaseFlashGraph C1C3C8twoPhaseFlashGraph.cpp)
create_phase_behavior_exe(BinaryMixturesPXPlot BinaryMixturesPXPlot.cpp)
create_phase_behavior_exe(GibbsFreeEnergyBinary GibbsFreeEnergyBinary.cpp)
create_phase_behavior_exe(FlashSPE3 FlashSPE3.cpp)
create_phase_behavior_exe(PhaseEnvelopeCH4C7 PhaseEnvelopeCH4-C7+.cpp)
create_phase_behavior_exe(C2C3PXPlot C2C3PXPlot.cpp)
create_phase_behavior_exe(C2C3PhaseEnvelope C2C3PhaseEnvelope.cpp)
create_phase_behavior_exe(C3nC4PXPlot C3nC4PXPlot.cpp)
create_phase_behavior_exe(C3nC4PhaseEnvelope C3nC4PhaseEnvelope.cpp)
create_phase_behavior_exe(FlashSPE3SI FlashSPE3SI.cpp)
create_phase_behavior_exe(CO2H2OtwoPhaseFlash CO2H2OtwoPhaseFlash.cpp)
create_phase_behavior_exe(C1C3C8PhaseEnvelope C1C3C8PhaseEnvelope.cpp)
create_phase_behavior_exe(C2nC5PhaseEnvelope C2nC5PhaseEnvelope.cpp)
create_phase_behavior_exe(SPE3TimeEnvelope SPE3TimeEnvelope.cpp)
create_phase_behavior_exe(SPE3UpscaledLiqFraqPlot SPE3UpscaledLiqFraqPlot.cpp)
create_phase_behavior_exe(SPE3CVDLiquidDropout SPE3CVDLiquidDropout.cpp)
find_package (Threads)

target_link_libraries(DensityValidation PRIVATE PhaseBehavior Threads::Threads)

target_include_directories(SPE3TimeEnvelope PRIVATE ${CMAKE_SOURCE_DIR}/external)
target_include_directories(SPE3UpscaledLiqFraqPlot PRIVATE ${CMAKE_SOURCE_DIR}/external)
target_link_libraries(SPE3TimeEnvelope PRIVATE PhaseBehavior Threads::Threads)

# IDEs should put the headers in a nice place
source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  PREFIX "Header Files"
  FILES ${HEADERS}
  )

source_group(
  TREE "${PROJECT_SOURCE_DIR}/external"
  PREFIX "Header Files"
  FILES "${PROJECT_SOURCE_DIR}/external/rapidcsv.h"
)