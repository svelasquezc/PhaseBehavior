function(create_phase_behavior_exe executable_name sources)
  add_executable(${executable_name} ${sources})
  # C++17 syntax
  target_compile_features(${executable_name} PRIVATE cxx_std_17)
  # Should be linked to the main library
  target_link_libraries(${executable_name} PRIVATE PhaseBehavior)
endfunction(create_phase_behavior_exe)

create_phase_behavior_exe(PhaseEnvelope PhaseEnvelope.cpp)
create_phase_behavior_exe(BlackOilEstimation BlackOilEstimation.cpp)
create_phase_behavior_exe(OffendingComposition OffendingComposition.cpp)
create_phase_behavior_exe(DensityValidation DensityValidation.cpp)
create_phase_behavior_exe(C1C3C8twoPhaseFlash C1C3C8twoPhaseFlash.cpp)
create_phase_behavior_exe(C1toC6twoPhaseFlash C1toC6twoPhaseFlash.cpp)
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
find_package (Threads)

target_link_libraries(DensityValidation PRIVATE PhaseBehavior Threads::Threads)

# IDEs should put the headers in a nice place
source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  PREFIX "Header Files"
  FILES ${HEADERS}
  )