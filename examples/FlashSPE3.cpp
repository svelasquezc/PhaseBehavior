#include <iostream>
#include <tuple>

#include <PhaseBehavior/Utilities/Types.hpp>
#include <PhaseBehavior/Component.hpp>
#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/Utilities/Input.hpp>
#include <PhaseBehavior/Phase.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>
#include <PhaseBehavior/PhaseEquilibrium.hpp>

using PR = PhaseBehavior::EoS::PR::PengRobinson;
using SRK = PhaseBehavior::EoS::SRK::SoaveRedlichKwong;

int main(){


    auto mixture = PhaseBehavior::Input::createMixtureFromFile<false>("PVT.csv", "InteractionCoefficients.csv");

    auto casePressure = static_cast<double>(2175.57); // psia (equal to 150 bar)
    auto caseTemperature = static_cast<double>(630); // R (equal to 400 K)

    auto result = PhaseBehavior::VaporLiquidEquilibrium::isothermalTwoPhaseFlash<PR>(mixture, casePressure, caseTemperature);   
    
    if (result == PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable){
        auto liquidPhase = PhaseBehavior::Phase::LiquidLikePhase(mixture);
        liquidPhase.molarVolume(mixture.compressibility("liquid"), casePressure, caseTemperature);
        auto densityL = liquidPhase.density();
        std::cout << "Liquid Density: "<< densityL << std::endl;
        auto gasPhase = PhaseBehavior::Phase::VaporLikePhase(mixture);
        gasPhase.molarVolume(mixture.compressibility("vapor"), casePressure, caseTemperature);
        auto densityG = gasPhase.density();
        std::cout << "Gas Density: "<< densityG << std::endl;
        std::cout << "liquid Composition: (";
        double pseudoMolarWL = 0.0;
        for (auto c : mixture){
            std::cout << c.composition("liquid")<<",";
            pseudoMolarWL += c.composition("liquid")*c->molarWeight()/1000;
        }
        std::cout << ")"<<std::endl << "gas Composition: (";
        double pseudoMolarWG = 0.0;
        for (auto c : mixture){
            std::cout << c.composition("vapor")<<",";
            pseudoMolarWG += c.composition("vapor")*c->molarWeight()/1000;
        }
        std::cout <<")"<< std::endl;

        auto gasMolarDensity = densityG/pseudoMolarWG;
        auto liquidMolarDensity = densityL/pseudoMolarWL;

        auto liquidVolumeFraction = 1/(1 + (mixture.molarFraction("vapor")/mixture.molarFraction("liquid"))*(liquidMolarDensity/gasMolarDensity));
    
        double COMSOLInitialDensityL = 1200.0;
        double COMSOLInitialDensityV = 200.0;


        double COMSOLInitialFractionL = (liquidMolarDensity*liquidVolumeFraction + gasMolarDensity*(1-liquidVolumeFraction) - COMSOLInitialDensityV)/(COMSOLInitialDensityL - COMSOLInitialDensityV);


        std::cout << "Gas Molar Density: " << gasMolarDensity << std::endl;
        std::cout << "Liquid Molar Density: " << liquidMolarDensity << std::endl;
        std::cout << "Liquid Volume Fraction: " << liquidVolumeFraction << std::endl;
        std::cout << "COMSOL Initial Liquid Fraction: " << COMSOLInitialFractionL << std::endl;
    };
}