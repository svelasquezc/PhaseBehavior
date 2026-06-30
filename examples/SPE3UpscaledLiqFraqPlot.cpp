#include <iostream>
#include <fstream>
#include <tuple>

#include <PhaseBehavior/Component.hpp>
#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/Phase.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>
#include <PhaseBehavior/PhaseEquilibrium.hpp>

#include <rapidcsv.h>

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

    PhaseBehavior::Mixture mixture{{C1, 1.0/6.0}, {C2, 1.0/6.0}, {C3, 1.0/6.0}, {C4, 1.0/6.0}, {C5, 1.0/6.0}, {C6, 1.0/6.0}};

    auto caseTemperature = 350.0 /*[K]*/;

    rapidcsv::Document doc("SPE3_UpscaledProps_th90.csv");

    std::vector<int> time = doc.GetColumn<int>("Time (s)");
    std::vector<double> avgP = doc.GetColumn<double>("Average Pressure (Pa)");
    std::vector<double> xC1 = doc.GetColumn<double>("C1 Composition (1)");
    std::vector<double> xC2 = doc.GetColumn<double>("C2 Composition (1)");
    std::vector<double> xC3 = doc.GetColumn<double>("C3 Composition (1)");
    std::vector<double> xC4 = doc.GetColumn<double>("C4 Composition (1)");
    std::vector<double> xC5 = doc.GetColumn<double>("C5 Composition (1)");
    std::vector<double> xC6 = doc.GetColumn<double>("C6 Composition (1)");

    std::ofstream outputFile;
    outputFile.open("SPE3-Upscaled-liquidFraction.csv");
    outputFile << "Pressure,LiquidFraction"<<std::endl;

    std::size_t pos=0;
    for(const auto& pressure : avgP){
        PhaseBehavior::Mixture mixture {{C1, xC1[pos]}, {C2, xC2[pos]}, {C3, xC3[pos]}, {C4, xC4[pos]}, {C5, xC5[pos]}, {C6, xC6[pos]}};

        auto casePressure = pressure/1000.0;

        auto result = isothermalTwoPhaseFlash<PR>(mixture, casePressure, caseTemperature);

        if (result == PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable){
            auto liquidPhase = PhaseBehavior::Phase::LiquidLikePhase(mixture);
            liquidPhase.molarVolume(mixture.compressibility(PhaseBehavior::PhaseName::liquid), casePressure, caseTemperature);
            auto densityL = liquidPhase.density();
            auto gasPhase = PhaseBehavior::Phase::VaporLikePhase(mixture);
            gasPhase.molarVolume(mixture.compressibility(PhaseBehavior::PhaseName::vapor), casePressure, caseTemperature);
            auto densityG = gasPhase.density();
            auto pseudoMolarWL = 0.0;
            for (auto c : mixture){
                pseudoMolarWL += c.composition(PhaseBehavior::PhaseName::liquid)*c->molarWeight()/1000;
            }
            auto pseudoMolarWG = 0.0;
            for (auto c : mixture){
                pseudoMolarWG += c.composition(PhaseBehavior::PhaseName::vapor)*c->molarWeight()/1000;
            }

            auto gasMolarDensity = densityG/pseudoMolarWG;
            auto liquidMolarDensity = densityL/pseudoMolarWL;

            auto liquidVolumeFraction = 1/(1 + (mixture.molarFraction(PhaseBehavior::PhaseName::vapor)/mixture.molarFraction(PhaseBehavior::PhaseName::liquid))*(liquidMolarDensity/gasMolarDensity));
        
            outputFile << casePressure << "," << liquidVolumeFraction << std::endl;
        }else{
            outputFile << casePressure << "," << 0.0 << std::endl;
        }
        result = isothermalTwoPhaseFlash<PR>(mixture, casePressure, caseTemperature);
        ++pos;
    };
    return 0;
}