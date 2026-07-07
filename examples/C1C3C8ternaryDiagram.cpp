#include <iostream>
#include <tuple>
#include <fstream>
#include <sstream>

#include <PhaseBehavior/Component.hpp>
#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/Phase.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>
#include <PhaseBehavior/PhaseEquilibrium.hpp>
#include <array>

using PR = PhaseBehavior::EoS::PR::PengRobinson;
using SRK = PhaseBehavior::EoS::SRK::SoaveRedlichKwong;


int main(){

    std::ofstream file("C1C3C8TernaryDiagram.csv");
    
    using PhaseBehavior::Component, PhaseBehavior::VaporLiquidEquilibrium::isothermalTwoPhaseFlash;

    auto C1 = Component {"CH4",      4599.2,    190.56,   6.1477929425e-3, 16.04 ,  0.0114};
    auto C3 = Component {"C3H8",     4251.2,    369.89,   4.5355587810e-3, 44.097,  0.1521};
    auto C8 = Component {"n-C8H18",  2497.0,    569.32,   4.2571306948e-3, 114.23,  0.3950};

    //auto casePressure = 9.5e3 /*[kPa]*/;
    //auto casePressure = 5.3e3 /*[kPa]*/;
    auto casePressure = 120e2 /*[kPa]*/;
    auto caseTemperature = 400.0 /*[K]*/;

    double initialComposition = 1e-5;
    auto endComposition = 1 - 2e-5;
    auto deltaComposition = (endComposition - initialComposition)/10000;

    file << "DensityL;DensityV;xC1;xC3;xC8;yC1;yC3;yC8;" <<std::endl;

    auto zC1 = (1 - 1e-5)*0.7;
    auto zC8 = (1 - 1e-5)*0.3;

    for(auto zC3 = initialComposition; zC3 < endComposition; zC3 += deltaComposition){

        auto totalComp = zC1 + zC3 + zC8;

        zC1 /= totalComp;
        zC3 /= totalComp;
        zC8 /= totalComp;

        PhaseBehavior::Mixture C1C3C8{{C1, zC1}, {C3, zC3}, {C8, zC8}};

        auto result = isothermalTwoPhaseFlash<PR>(C1C3C8, casePressure, caseTemperature);

        if (result == PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable){

            auto liquidPhase = PhaseBehavior::Phase::LiquidLikePhase(C1C3C8);
            liquidPhase.molarVolume(C1C3C8.compressibility(PhaseBehavior::PhaseName::liquid), casePressure, caseTemperature);
            auto densityL = liquidPhase.density();
            file << densityL << ";";
            auto gasPhase = PhaseBehavior::Phase::VaporLikePhase(C1C3C8);
            gasPhase.molarVolume(C1C3C8.compressibility(PhaseBehavior::PhaseName::vapor), casePressure, caseTemperature);
            auto densityG = gasPhase.density();
            file << densityG << ";";
            double pseudoMolarWL = 0.0;
            std::array<double, 3> x, y;
            std::size_t i = 0;
            for (auto c : C1C3C8){
                x[i] = c.composition(PhaseBehavior::PhaseName::liquid);
                file << x[i]<<";";
                ++i;
            }
            double pseudoMolarWG = 0.0;
            i=0;
            for (auto c : C1C3C8){
                y[i] = c.composition(PhaseBehavior::PhaseName::vapor);
                file << y[i]<<";";
                ++i;
            }

            zC1 = (x[0] + y[0])/2;
            zC3 = (x[1] + y[1])/2;
            zC8 = (x[2] + y[2])/2;

            file << std::endl;
        }else{
            break;
        }
        
    }

    file.close();
    return 0;
}