#include <iostream>
#include <tuple>

#include <PhaseBehavior/Component.hpp>
#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/Phase.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>
#include <PhaseBehavior/PhaseEquilibrium.hpp>

using PR = PhaseBehavior::EoS::PR::PengRobinson;
using SRK = PhaseBehavior::EoS::SRK::SoaveRedlichKwong;


int main(){

    using PhaseBehavior::Component, PhaseBehavior::VaporLiquidEquilibrium::isothermalTwoPhaseFlash;
    
    auto C1 = PhaseBehavior::Component {"CH4",      666.4,     343.33,   0.0988, 16.043,  0.0104};
    auto C7plus = PhaseBehavior::Component {"C7+",  305.2,     1112,     0.0679, 142.285,  0.4898};

    PhaseBehavior::Mixture mixture{{C1, 0.5}, {C7plus, 0.5}};

    auto casePressure = 2300 /*[psia]*/;
    //auto casePressure = 5.3e3 /*[kPa]*/;
    //auto casePressure = 1.5e3 /*[kPa]*/;
    auto caseTemperature = 720 /*[R]*/;


    auto result = isothermalTwoPhaseFlash<PR>(mixture, casePressure, caseTemperature);


    if (result == PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable){
        auto liquidPhase = PhaseBehavior::Phase::LiquidLikePhase(mixture);
        liquidPhase.molarVolume(mixture.compressibility(PhaseBehavior::PhaseName::liquid), casePressure, caseTemperature);
        auto densityL = liquidPhase.density();
        std::cout << "Liquid Density: "<< densityL << std::endl;
        auto gasPhase = PhaseBehavior::Phase::VaporLikePhase(mixture);
        gasPhase.molarVolume(mixture.compressibility(PhaseBehavior::PhaseName::vapor), casePressure, caseTemperature);
        auto densityG = gasPhase.density();
        std::cout << "Gas Density: "<< densityG << std::endl;
        std::cout << "liquid Composition: (";
        double pseudoMolarWL = 0.0;
        for (auto c : mixture){
            std::cout << c.composition(PhaseBehavior::PhaseName::liquid)<<",";
            pseudoMolarWL += c.composition(PhaseBehavior::PhaseName::liquid)*c->molarWeight()/1000;
        }
        std::cout << ")"<<std::endl << "gas Composition: (";
        double pseudoMolarWG = 0.0;
        for (auto c : mixture){
            std::cout << c.composition(PhaseBehavior::PhaseName::vapor)<<",";
            pseudoMolarWG += c.composition(PhaseBehavior::PhaseName::vapor)*c->molarWeight()/1000;
        }
        std::cout <<")"<< std::endl;

        auto gasMolarDensity = densityG/pseudoMolarWG;
        auto liquidMolarDensity = densityL/pseudoMolarWL;

        auto liquidVolumeFraction = 1/(1 + (mixture.molarFraction(PhaseBehavior::PhaseName::vapor)/mixture.molarFraction(PhaseBehavior::PhaseName::liquid))*(liquidMolarDensity/gasMolarDensity));
    
        double COMSOLInitialDensityL = 13000.0;
        double COMSOLInitialDensityV = 100.0;


        double COMSOLInitialFractionL = (liquidMolarDensity*liquidVolumeFraction + gasMolarDensity*(1-liquidVolumeFraction) - COMSOLInitialDensityV)/(COMSOLInitialDensityL - COMSOLInitialDensityV);


        std::cout << "Gas Molar Density: " << gasMolarDensity << std::endl;
        std::cout << "Liquid Molar Density: " << liquidMolarDensity << std::endl;
        std::cout << "Liquid Volume Fraction: " << liquidVolumeFraction << std::endl;
        std::cout << "COMSOL Initial Liquid Fraction: " << COMSOLInitialFractionL << std::endl;
    };
    return 0;
}