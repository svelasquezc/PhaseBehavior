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
create_phase_behavior_exe(BinaryMixturesPXPlot BinaryMixturesPXPlot.cpp)

find_package (Threads)

target_link_libraries(DensityValidation PRIVATE PhaseBehavior Threads::Threads)

# IDEs should put the headers in a nice place
source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  PREFIX "Header Files"
  FILES ${HEADERS}
  )