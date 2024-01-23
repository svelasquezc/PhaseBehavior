#include <iostream>

#include <PhaseBehavior/Utilities/Input.hpp>
#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>
#include <PhaseBehavior/PhaseEnvelope.hpp>



int main(){
    auto mixture = PhaseBehavior::Input::createMixtureFromFile<true>("PVTPR.csv", "InteractionCoefficients.csv");


    auto temperatureResult = PhaseBehavior::PhaseEnvelope<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 0);


    return 0;
}