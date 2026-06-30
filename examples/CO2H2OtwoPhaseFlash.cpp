#include <iostream>
#include <tuple>

#include <PhaseBehavior/Component.hpp>
#include <PhaseBehavior/GERG/Components.hpp>
#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/Phase.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>
#include <PhaseBehavior/PhaseEquilibrium.hpp>

using PR = PhaseBehavior::EoS::PR::PengRobinson;
using SRK = PhaseBehavior::EoS::SRK::SoaveRedlichKwong;


int main(){
    
    using PhaseBehavior::Component, PhaseBehavior::VaporLiquidEquilibrium::isothermalTwoPhaseFlash;

    auto CO2 = PhaseBehavior::EoS::GERG::Components::CO2;
    auto H2O = PhaseBehavior::EoS::GERG::Components::H2O;

    PhaseBehavior::Mixture CO2H2O{{CO2, 0.3}, {H2O, 0.7}};

    auto casePressure = 3e2 /*[kPa]*/;
    auto caseTemperature = 380.0 /*[K]*/;


    auto result = isothermalTwoPhaseFlash<PR>(CO2H2O, casePressure, caseTemperature);


    if (result == PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable){
        auto liquidPhase = PhaseBehavior::Phase::LiquidLikePhase(CO2H2O);
        liquidPhase.molarVolume(CO2H2O.compressibility(PhaseBehavior::PhaseName::liquid), casePressure, caseTemperature);
        auto densityL = liquidPhase.density();
        std::cout << "Liquid Density: "<< densityL << std::endl;
        auto gasPhase = PhaseBehavior::Phase::VaporLikePhase(CO2H2O);
        gasPhase.molarVolume(CO2H2O.compressibility(PhaseBehavior::PhaseName::vapor), casePressure, caseTemperature);
        auto densityG = gasPhase.density();
        std::cout << "Gas Density: "<< densityG << std::endl;
        std::cout << "liquid Composition: (";
        double pseudoMolarWL = 0.0;
        for (auto c : CO2H2O){
            std::cout << c.composition(PhaseBehavior::PhaseName::liquid)<<",";
            pseudoMolarWL += c.composition(PhaseBehavior::PhaseName::liquid)*c->molarWeight()/1000;
        }
        std::cout << ")"<<std::endl << "gas Composition: (";
        double pseudoMolarWG = 0.0;
        for (auto c : CO2H2O){
            std::cout << c.composition(PhaseBehavior::PhaseName::vapor)<<",";
            pseudoMolarWG += c.composition(PhaseBehavior::PhaseName::vapor)*c->molarWeight()/1000;
        }
        std::cout <<")"<< std::endl;

        auto gasMolarDensity = densityG/pseudoMolarWG;
        auto liquidMolarDensity = densityL/pseudoMolarWL;

        auto liquidVolumeFraction = 1/(1 + (CO2H2O.molarFraction(PhaseBehavior::PhaseName::vapor)/CO2H2O.molarFraction(PhaseBehavior::PhaseName::liquid))*(liquidMolarDensity/gasMolarDensity));
    
        double COMSOLInitialDensityL = 45000.0;
        double COMSOLInitialDensityV = 10.0;


        double COMSOLInitialFractionL = (liquidMolarDensity*liquidVolumeFraction + gasMolarDensity*(1-liquidVolumeFraction) - COMSOLInitialDensityV)/(COMSOLInitialDensityL - COMSOLInitialDensityV);


        std::cout << "Gas Molar Density: " << gasMolarDensity << std::endl;
        std::cout << "Liquid Molar Density: " << liquidMolarDensity << std::endl;
        std::cout << "Liquid Volume Fraction: " << liquidVolumeFraction << std::endl;
        std::cout << "COMSOL Initial Liquid Fraction: " << COMSOLInitialFractionL << std::endl;
    };
    return 0;
}