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

    PhaseBehavior::Component CO2   = {"CO2",      7384.3,    304.39,   2.1385799828e-3,  44.01,   0.2667};
    PhaseBehavior::Component C1    = {"CH4",      4594.7,    190.74,   6.1477929425e-3,  16.043,  0.0104};
    PhaseBehavior::Component C2    = {"C2H6",     4871.1,    305.51,   4.8402710552e-3,  30.070,  0.0979};
    PhaseBehavior::Component C3    = {"C3H8",     4247.2,    370.03,   4.5355587810e-3,  44.097,  0.1522};
    PhaseBehavior::Component i_C4  = {"i-C4H10",  3639.7,    408.03,   4.4571224820e-3,  58.123,  0.1852};
    PhaseBehavior::Component n_C4  = {"n-C4H10",  3796.3,    425.34,   4.3890449440e-3,  58.123,  0.1995};
    PhaseBehavior::Component i_C5  = {"i-C5H12",  3381.2,    460.61,   4.2373122060e-3,  72.150,  0.2280};
    PhaseBehavior::Component n_C5  = {"n-C5H12",  3368.8,    469.89,   4.3103448281e-3,  72.150,  0.2514};
    PhaseBehavior::Component n_C6  = {"n-C6H14",  3012.3,    507.56,   4.2885324643e-3,  86.177,  0.2994};
    PhaseBehavior::Component C7    = {"C7+",      2104.3,    617.78,   4.3103448273e-3,  142.285, 0.4898};

    PhaseBehavior::Mixture mixture {{CO2, 0.0031}, {C1, 0.6192}, {C2, 0.1408}, {C3, 0.0835}, {i_C4, 0.0097}, {n_C4, 0.0341}, {i_C5, 0.0084}, {n_C5, 0.0148}, {n_C6, 0.0179}, {C7, 0.0685}};

    mixture.interactionCoefficient("CO2", "CH4", 0.105);
    mixture.interactionCoefficient("CH4", "CO2", 0.105);
    mixture.interactionCoefficient("CO2", "C2H6", 0.13);
    mixture.interactionCoefficient("C2H6", "CO2", 0.13);
    mixture.interactionCoefficient("CO2", "C3H8", 0.125);
    mixture.interactionCoefficient("C3H8", "CO2", 0.125);
    mixture.interactionCoefficient("CO2", "i-C4H10", 0.12);
    mixture.interactionCoefficient("i-C4H10", "CO2", 0.12);
    mixture.interactionCoefficient("CO2", "n-C4H10", 0.115);
    mixture.interactionCoefficient("n-C4H10", "CO2", 0.115);
    mixture.interactionCoefficient("CO2", "i-C5H12", 0.115);
    mixture.interactionCoefficient("i-C5H12", "CO2", 0.115);
    mixture.interactionCoefficient("CO2", "n-C5H12", 0.115);
    mixture.interactionCoefficient("n-C5H12", "CO2", 0.115);
    mixture.interactionCoefficient("CO2", "n-C6H14", 0.115);
    mixture.interactionCoefficient("n-C6H14", "CO2", 0.115);
    mixture.interactionCoefficient("CO2", "C7+", 0.115);
    mixture.interactionCoefficient("C7+", "CO2", 0.115);
    

    auto casePressure = static_cast<double>(15000); // kPa (equal to 150 bar)
    auto caseTemperature = static_cast<double>(325); // K (equal to 400 K)

    auto result = PhaseBehavior::VaporLiquidEquilibrium::isothermalTwoPhaseFlash<PR>(mixture, casePressure, caseTemperature);   
    
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
        double COMSOLInitialDensityV = 200.0;


        double COMSOLInitialFractionL = (liquidMolarDensity*liquidVolumeFraction + gasMolarDensity*(1-liquidVolumeFraction) - COMSOLInitialDensityV)/(COMSOLInitialDensityL - COMSOLInitialDensityV);


        std::cout << "Gas Molar Density: " << gasMolarDensity << std::endl;
        std::cout << "Liquid Molar Density: " << liquidMolarDensity << std::endl;
        std::cout << "Liquid Volume Fraction: " << liquidVolumeFraction << std::endl;
        std::cout << "COMSOL Initial Liquid Fraction: " << COMSOLInitialFractionL << std::endl;
    };
}