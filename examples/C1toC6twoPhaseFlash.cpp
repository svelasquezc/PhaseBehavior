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

    auto C1 = Component {"C1",      4599.2,    190.56,   6.1477929425e-3, 16.04,  0.0114};
    auto C2 = Component {"C2",      4872.2,    305.32,   4.8402710552e-3, 30.07,  0.0990};
    auto C3 = Component {"C3",      4251.2,    369.89,   4.5355587810e-3, 44.10, 0.1521};
    auto C4 = Component {"C4",      3796.3,    425.34,   4.3890449440e-3, 58.12, 0.1995};
    auto C5 = Component {"C5",      3368.8,    469.89,   4.3103448281e-3, 72.15,  0.2514};
    auto C6 = Component {"C6",      3012.3,    507.56,   4.2885324643e-3, 86.18, 0.2994};

    PhaseBehavior::Mixture C1toC6{{C1, 1.0/6.0}, {C2, 1.0/6.0}, {C3, 1.0/6.0}, {C4, 1.0/6.0}, {C5, 1.0/6.0}, {C6, 1.0/6.0}};

    auto casePressure = 5e3 /*[kPa]*/;
    //auto casePressure = 5.3e3 /*[kPa]*/;
    //auto casePressure = 1.5e3 /*[kPa]*/;
    auto caseTemperature = 400.0 /*[K]*/;


    auto result = isothermalTwoPhaseFlash<PR>(C1toC6, casePressure, caseTemperature);


    if (result == PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable){
        auto liquidPhase = PhaseBehavior::Phase::LiquidLikePhase(C1toC6);
        liquidPhase.molarVolume(C1toC6.compressibility("liquid"), casePressure, caseTemperature);
        auto densityL = liquidPhase.density();
        std::cout << "Liquid Density: "<< densityL << std::endl;
        auto gasPhase = PhaseBehavior::Phase::VaporLikePhase(C1toC6);
        gasPhase.molarVolume(C1toC6.compressibility("vapor"), casePressure, caseTemperature);
        auto densityG = gasPhase.density();
        std::cout << "Gas Density: "<< densityG << std::endl;
        std::cout << "liquid Composition: (";
        double pseudoMolarWL = 0.0;
        for (auto c : C1toC6){
            std::cout << c.composition("liquid")<<",";
            pseudoMolarWL += c.composition("liquid")*c->molarWeight()/1000;
        }
        std::cout << ")"<<std::endl << "gas Composition: (";
        double pseudoMolarWG = 0.0;
        for (auto c : C1toC6){
            std::cout << c.composition("vapor")<<",";
            pseudoMolarWG += c.composition("vapor")*c->molarWeight()/1000;
        }
        std::cout <<")"<< std::endl;

        auto gasMolarDensity = densityG/pseudoMolarWG;
        auto liquidMolarDensity = densityL/pseudoMolarWL;

        auto liquidVolumeFraction = 1/(1 + (C1toC6.molarFraction("vapor")/C1toC6.molarFraction("liquid"))*(liquidMolarDensity/gasMolarDensity));
    
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