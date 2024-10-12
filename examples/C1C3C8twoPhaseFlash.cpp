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

    auto C1 = Component {"CH4",      4599.2,    190.56,   6.1477929425e-3, 16.04 ,  0.0114};
    auto C3 = Component {"C3H8",     4251.2,    369.89,   4.5355587810e-3, 44.097,  0.1521};
    auto C8 = Component {"n-C8H18",  2497.0,    569.32,   4.2571306948e-3, 114.23,  0.3950};

    PhaseBehavior::Mixture C1C3C8{{C1, 1.0/3.0}, {C3, 1.0/3.0}, {C8, 1.0/3.0}};

    auto result = isothermalTwoPhaseFlash<PR>(C1C3C8, 7.35e2 /*[kPa]*/, 350 /*[K]*/);


    if (result == PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable){
        auto liquidPhase = PhaseBehavior::Phase::LiquidLikePhase(C1C3C8);
        liquidPhase.molarVolume(C1C3C8.compressibility("liquid"), 7.35e2, 350);
        auto densityL = liquidPhase.density();
        std::cout << "Liquid Density: "<< densityL << std::endl;
        auto gasPhase = PhaseBehavior::Phase::VaporLikePhase(C1C3C8);
        gasPhase.molarVolume(C1C3C8.compressibility("vapor"), 7.35e2, 350);
        auto densityG = gasPhase.density();
        std::cout << "Gas Density: "<< densityG << std::endl;
        std::cout << "liquid Composition: (";
        double pseudoMolarWL = 0.0;
        for (auto c : C1C3C8){
            std::cout << c.composition("liquid")<<",";
            pseudoMolarWL += c.composition("liquid")*c->molarWeight()/1000;
        }
        std::cout << ")"<<std::endl << "gas Composition: (";
        double pseudoMolarWG = 0.0;
        for (auto c : C1C3C8){
            std::cout << c.composition("vapor")<<",";
            pseudoMolarWG += c.composition("vapor")*c->molarWeight()/1000;
        }
        std::cout <<")"<< std::endl;

        std::cout << "Gas Molar Density: " << densityG/pseudoMolarWG << std::endl;
        std::cout << "Liquid Molar Density: " << densityL/pseudoMolarWL << std::endl;
    };
    return 0;
}