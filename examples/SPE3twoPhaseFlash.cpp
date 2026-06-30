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
    using PhaseBehavior::Component;

    auto CO2 = Component{"CO2",       7384.3,    304.39,   2.1385799828e-3, 44.01,   0.2667};
    auto C1 = Component {"C1",        4594.7,    190.74,   6.1477929425e-3, 16.043,  0.0104};
    auto C2 = Component {"C2",        4871.1,    305.51,   4.8402710552e-3, 30.07,  0.0979};
    auto C3 = Component {"C3",        4247.2,    370.03,   4.5355587810e-3, 44.097, 0.1522};
    auto iC4 = Component{"iC4",       3639.7,    407.817,  4.4571224820e-3, 58.1222,   0.183};
    auto nC4 = Component {"nC4",       3796.3,    425.34,  4.3890449440e-3, 58.12, 0.1995};
    auto iC5 = Component{"iC5",       3381.2,    460.35,   4.2373122060e-3, 72.14878,  0.227};
    auto nC5 = Component {"nC5",      3368.8,    469.89,   4.3103448281e-3, 72.15,  0.2514};
    auto C6 = Component {"C6",        3012.3,    507.56,   4.2885324643e-3, 86.18, 0.2994};
    auto C7p = Component{"C7+",       2104.3,    617.78,   4.238848831e-3, 142.285, 0.4898};

    PhaseBehavior::Mixture SPE3{{CO2, 0.0031}, {C1, 0.6192}, {C2, 0.1408}, {C3, 0.0835}, {iC4, 0.0097}, {nC4, 0.0341}, {iC5, 0.0084}, {nC5, 0.0148}, {C6, 0.0179}, {C7p, 0.0685}};

    SPE3.interactionCoefficient("CO2", "C1", 0.105); SPE3.interactionCoefficient("C1", "CO2", 0.105);
    SPE3.interactionCoefficient("CO2", "C2", 0.13); SPE3.interactionCoefficient("C2", "CO2", 0.13);
    SPE3.interactionCoefficient("CO2", "C3", 0.125); SPE3.interactionCoefficient("C3", "CO2", 0.125);
    SPE3.interactionCoefficient("CO2", "iC4", 0.12); SPE3.interactionCoefficient("iC4", "CO2", 0.12);
    SPE3.interactionCoefficient("CO2", "nC4", 0.115); SPE3.interactionCoefficient("nC4", "CO2", 0.115);
    SPE3.interactionCoefficient("CO2", "iC5", 0.115); SPE3.interactionCoefficient("iC5", "CO2", 0.115);
    SPE3.interactionCoefficient("CO2", "nC5", 0.115); SPE3.interactionCoefficient("nC5", "CO2", 0.115);
    SPE3.interactionCoefficient("CO2", "C6", 0.115); SPE3.interactionCoefficient("C6", "CO2", 0.115);
    SPE3.interactionCoefficient("CO2", "C7+", 0.115); SPE3.interactionCoefficient("C7+", "CO2", 0.115);


    auto casePressure = 15e3 /*[kPa]*/;
    //auto casePressure = 5.3e3 /*[kPa]*/;
    //auto casePressure = 1.5e3 /*[kPa]*/;
    auto caseTemperature = 325.0 /*[K]*/;


    auto result = isothermalTwoPhaseFlash<PR>(SPE3, casePressure, caseTemperature);


    if (result == PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable){
        auto liquidPhase = PhaseBehavior::Phase::LiquidLikePhase(SPE3);
        liquidPhase.molarVolume(SPE3.compressibility(PhaseBehavior::PhaseName::liquid), casePressure, caseTemperature);
        auto densityL = liquidPhase.density();
        std::cout << "Liquid Density: "<< densityL << std::endl;
        auto gasPhase = PhaseBehavior::Phase::VaporLikePhase(SPE3);
        gasPhase.molarVolume(SPE3.compressibility(PhaseBehavior::PhaseName::vapor), casePressure, caseTemperature);
        auto densityG = gasPhase.density();
        std::cout << "Gas Density: "<< densityG << std::endl;
        std::cout << "liquid Composition: (";
        double pseudoMolarWL = 0.0;
        for (auto c : SPE3){
            std::cout << c.composition(PhaseBehavior::PhaseName::liquid)<<",";
            pseudoMolarWL += c.composition(PhaseBehavior::PhaseName::liquid)*c->molarWeight()/1000;
        }
        std::cout << ")"<<std::endl << "gas Composition: (";
        double pseudoMolarWG = 0.0;
        for (auto c : SPE3){
            std::cout << c.composition(PhaseBehavior::PhaseName::vapor)<<",";
            pseudoMolarWG += c.composition(PhaseBehavior::PhaseName::vapor)*c->molarWeight()/1000;
        }
        std::cout <<")"<< std::endl;

        auto gasMolarDensity = densityG/pseudoMolarWG;
        auto liquidMolarDensity = densityL/pseudoMolarWL;

        auto liquidVolumeFraction = 1/(1 + (SPE3.molarFraction(PhaseBehavior::PhaseName::vapor)/SPE3.molarFraction(PhaseBehavior::PhaseName::liquid))*(liquidMolarDensity/gasMolarDensity));
    
        double COMSOLInitialDensityL = 20000.0;
        double COMSOLInitialDensityV = 200.0;


        double COMSOLInitialFractionL = (liquidMolarDensity*liquidVolumeFraction + gasMolarDensity*(1-liquidVolumeFraction) - COMSOLInitialDensityV)/(COMSOLInitialDensityL - COMSOLInitialDensityV);


        std::cout << "Gas Molar Density: " << gasMolarDensity << std::endl;
        std::cout << "Liquid Molar Density: " << liquidMolarDensity << std::endl;
        std::cout << "Liquid Volume Fraction: " << liquidVolumeFraction << std::endl;
        std::cout << "COMSOL Initial Liquid Fraction: " << COMSOLInitialFractionL << std::endl;
    };
    return 0;
}