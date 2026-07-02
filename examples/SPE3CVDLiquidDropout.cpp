#include <iostream>
#include <fstream>
#include <tuple>

#include <PhaseBehavior/Component.hpp>
#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/Phase.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>
#include <PhaseBehavior/PhaseEquilibrium.hpp>

using PR = PhaseBehavior::EoS::PR::PengRobinson;
using SRK = PhaseBehavior::EoS::SRK::SoaveRedlichKwong;


int main(){
    
    using PhaseBehavior::Component, PhaseBehavior::PhaseName, PhaseBehavior::VaporLiquidEquilibrium::isothermalTwoPhaseFlash;

    std::ofstream outputFile;
    outputFile.open("SPE3-CVD-liquidFraction.csv");
    outputFile << "Pressure,LiquidFraction,zC1,zC2,zC3,zC4,zC5,zC6"<<std::endl;

    auto C1 = Component {"C1",      4599.2,    190.56,   6.1477929425e-3, 16.04,  0.0114};
    auto C2 = Component {"C2",      4872.2,    305.32,   4.8402710552e-3, 30.07,  0.0990};
    auto C3 = Component {"C3",      4251.2,    369.89,   4.5355587810e-3, 44.10, 0.1521};
    auto C4 = Component {"C4",      3796.3,    425.34,   4.3890449440e-3, 58.12, 0.1995};
    auto C5 = Component {"C5",      3368.8,    469.89,   4.3103448281e-3, 72.15,  0.2514};
    auto C6 = Component {"C6",      3012.3,    507.56,   4.2885324643e-3, 86.18, 0.2994};

    PR eos;

    auto initialPressure = 14e3 /*[kPa]*/;
    auto endPressure = 1e2 /*[kPa]*/;
    auto caseTemperature = 350.0 /*[K]*/;

    double totalMoles = 1;

    std::array<double, 6> composition = {0.6223, 0.1408, 0.0835, 0.0438, 0.0232, 0.0864};

    PhaseBehavior::Mixture initialMixture {{C1, composition[0]}, {C2, composition[1]}, {C3, composition[2]}, {C4, composition[3]}, {C5, composition[4]}, {C6, composition[5]}};
    
    eos(initialMixture, initialPressure, caseTemperature);
    initialMixture.compressibility(PhaseName::global, eos.selectedCompressibility());
    auto fluid = PhaseBehavior::Phase::singlePhaseIdentification(initialMixture, initialMixture.compressibility(PhaseName::global), initialPressure, caseTemperature, eos);

    auto pvtMolecularWeight = fluid->molecularWeight();
    auto pvtDensity = fluid->density();

    auto PVTCellVolume = totalMoles * pvtMolecularWeight/pvtDensity;

    auto currentMoles = totalMoles;

    for(auto pressure = initialPressure; pressure >= endPressure; pressure -= 5.0){

        PhaseBehavior::Mixture mixture {{C1, composition[0]}, {C2, composition[1]}, {C3, composition[2]}, {C4, composition[3]}, {C5, composition[4]}, {C6, composition[5]}};

        std::array<double, 6> molesPerComponent = {currentMoles*composition[0], currentMoles*composition[1], currentMoles*composition[2], currentMoles*composition[3], currentMoles*composition[4], currentMoles*composition[5]};

        auto result = isothermalTwoPhaseFlash<PR>(mixture, pressure, caseTemperature);

        if (result == PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable){
            auto liquidPhase = PhaseBehavior::Phase::LiquidLikePhase(mixture);
            liquidPhase.molarVolume(mixture.compressibility(PhaseBehavior::PhaseName::liquid), pressure, caseTemperature);
            auto liquidDensity = liquidPhase.density();
            auto gasPhase = PhaseBehavior::Phase::VaporLikePhase(mixture);
            gasPhase.molarVolume(mixture.compressibility(PhaseBehavior::PhaseName::vapor), pressure, caseTemperature);
            auto gasDensity = gasPhase.density();
            
            std::array<double, 6> x, y;
            std::size_t i = 0;
            auto pseudoMolarWL = 0.0;
            for (auto c : mixture){
                x[i] = c.composition(PhaseBehavior::PhaseName::liquid);
                pseudoMolarWL += x[i]*c->molarWeight()/1000; // [kg/mol]
                ++i;
            }

            auto pseudoMolarWG = 0.0;
            i=0;
            for (auto c : mixture){
                y[i] = c.composition(PhaseBehavior::PhaseName::vapor);
                pseudoMolarWG += y[i]*c->molarWeight()/1000; // [kg/mol]
                ++i;
            }

            auto gasMolarDensity = gasDensity/pseudoMolarWG;
            auto liquidMolarDensity = liquidDensity/pseudoMolarWL;

            auto liquidVolumeFraction = 1/(1 + (mixture.molarFraction(PhaseBehavior::PhaseName::vapor)/mixture.molarFraction(PhaseBehavior::PhaseName::liquid))*(liquidMolarDensity/gasMolarDensity));

            auto oilMoles = currentMoles*mixture.molarFraction(PhaseBehavior::PhaseName::liquid);

            auto gasPlusExcessMoles = currentMoles*mixture.molarFraction(PhaseBehavior::PhaseName::vapor);
            auto gasPlusExcessVolume = gasPlusExcessMoles*pseudoMolarWG/gasDensity;
            auto oilVolume = oilMoles*pseudoMolarWL/liquidDensity;

            auto gasVolume = PVTCellVolume - oilVolume;
            auto excessVolume = gasPlusExcessVolume - gasVolume;

            auto gasMoles = gasVolume*gasDensity/pseudoMolarWG;
            auto excessMoles = excessVolume*gasDensity/pseudoMolarWG;

            std::array<double, 6> excessMolesPerComponent = {excessMoles*y[0], excessMoles*y[1], excessMoles*y[2], excessMoles*y[3], excessMoles*y[4], excessMoles*y[5]};

            currentMoles -= excessMoles;

            molesPerComponent = {molesPerComponent[0] - excessMolesPerComponent[0], molesPerComponent[1] - excessMolesPerComponent[1], molesPerComponent[2] - excessMolesPerComponent[2], molesPerComponent[3] - excessMolesPerComponent[3], molesPerComponent[4] - excessMolesPerComponent[4], molesPerComponent[5] - excessMolesPerComponent[5]};

            composition = {molesPerComponent[0]/currentMoles, molesPerComponent[1]/currentMoles, molesPerComponent[2]/currentMoles, molesPerComponent[3]/currentMoles, molesPerComponent[4]/currentMoles, molesPerComponent[5]/currentMoles};
        
            outputFile << pressure << "," << liquidVolumeFraction << "," << composition[0] << "," << composition[1] << "," << composition[2] << "," << composition[3] << "," << composition[4] << "," << composition[5] << std::endl;
        }else{
            eos(mixture, pressure, caseTemperature);
            mixture.compressibility(PhaseName::global, eos.selectedCompressibility());
            auto fluid = PhaseBehavior::Phase::singlePhaseIdentification(mixture, mixture.compressibility(PhaseName::global), pressure, caseTemperature, eos);

            auto molecularWeight = fluid->molecularWeight();
            auto density = fluid->density();

            auto newVolume = currentMoles * molecularWeight/density;
            auto excessVolume = newVolume - PVTCellVolume;
            auto excessMoles = excessVolume * density/molecularWeight;
            std::array<double, 6> excessMolesPerComponent = {excessMoles*composition[0], excessMoles*composition[1], excessMoles*composition[2], excessMoles*composition[3], excessMoles*composition[4], excessMoles*composition[5]};
            currentMoles -= excessMoles;

            molesPerComponent = {molesPerComponent[0] - excessMolesPerComponent[0], molesPerComponent[1] - excessMolesPerComponent[1], molesPerComponent[2] - excessMolesPerComponent[2], molesPerComponent[3] - excessMolesPerComponent[3], molesPerComponent[4] - excessMolesPerComponent[4], molesPerComponent[5] - excessMolesPerComponent[5]};

            composition = {molesPerComponent[0]/currentMoles, molesPerComponent[1]/currentMoles, molesPerComponent[2]/currentMoles, molesPerComponent[3]/currentMoles, molesPerComponent[4]/currentMoles, molesPerComponent[5]/currentMoles};
            outputFile << pressure << "," << 0.0 << "," << composition[0] << "," << composition[1] << "," << composition[2] << "," << composition[3] << "," << composition[4] << "," << composition[5] << std::endl;
        }
    };
    outputFile.close();
    return 0;
}