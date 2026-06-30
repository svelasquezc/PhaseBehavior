#include <iostream>
#include <tuple>
#include <fstream>
#include <sstream>

#include <PhaseBehavior/Component.hpp>
#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/Phase.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>
#include <PhaseBehavior/PhaseEquilibrium.hpp>

using PR = PhaseBehavior::EoS::PR::PengRobinson;
using SRK = PhaseBehavior::EoS::SRK::SoaveRedlichKwong;


int main(){

    std::ofstream file("FlashData.csv");
    
    using PhaseBehavior::Component, PhaseBehavior::VaporLiquidEquilibrium::isothermalTwoPhaseFlash;

    auto C1 = Component {"CH4",      4599.2,    190.56,   6.1477929425e-3, 16.04 ,  0.0114};
    auto C3 = Component {"C3H8",     4251.2,    369.89,   4.5355587810e-3, 44.097,  0.1521};
    auto C8 = Component {"n-C8H18",  2497.0,    569.32,   4.2571306948e-3, 114.23,  0.3950};

    PhaseBehavior::Mixture C1C3C8{{C1, 1.0/3.0}, {C3, 1.0/3.0}, {C8, 1.0/3.0}};

    //auto casePressure = 9.5e3 /*[kPa]*/;
    //auto casePressure = 5.3e3 /*[kPa]*/;
    auto casePressure = 9.5e2 /*[kPa]*/;
    auto caseTemperature = 400.0 /*[K]*/;

    double initialPressure = 101.325;

    casePressure = initialPressure;

    file << "Pressure;DensityL;DensityV;xC1;xC3;xC8;yC1;yC3;yC8;" <<std::endl;

    while(casePressure <= 10000){

        auto result = isothermalTwoPhaseFlash<PR>(C1C3C8, casePressure, caseTemperature);

        if (result == PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable){

            file << casePressure <<";";

            auto liquidPhase = PhaseBehavior::Phase::LiquidLikePhase(C1C3C8);
            liquidPhase.molarVolume(C1C3C8.compressibility(PhaseBehavior::PhaseName::liquid), casePressure, caseTemperature);
            auto densityL = liquidPhase.density();
            file << densityL << ";";
            auto gasPhase = PhaseBehavior::Phase::VaporLikePhase(C1C3C8);
            gasPhase.molarVolume(C1C3C8.compressibility(PhaseBehavior::PhaseName::vapor), casePressure, caseTemperature);
            auto densityG = gasPhase.density();
            file << densityG << ";";
            double pseudoMolarWL = 0.0;
            for (auto c : C1C3C8){
                file << c.composition(PhaseBehavior::PhaseName::liquid)<<";";
            }
            double pseudoMolarWG = 0.0;
            for (auto c : C1C3C8){
                file << c.composition(PhaseBehavior::PhaseName::vapor)<<";";
            }
            file << std::endl;
        };

        casePressure += 100;
    }

    file.close();
    return 0;
}