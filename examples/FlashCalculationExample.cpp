#include <tuple>
#include <iostream>
#include <PhaseBehavior/Component.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>
#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/PhaseEquilibrium.hpp>

int main(){

    using PhaseBehavior::VaporLiquidEquilibrium::isothermalTwoPhaseFlash;
    using PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult;
    using PhaseBehavior::EoS::PR::PengRobinson;
    PengRobinson eos;
                                  /*name*/  /*Pc (psia)*/        /*Tc (R)*/      /*Vc (ft3/lbm)*/  /*Mw (lbm/lbmol)*/ /*accentricF*/ /*shift*/
    PhaseBehavior::Component C1   {"C1",    667.2,           -116.59 + 459.67,       0.0991,          16.043,              0.008,     0};
    PhaseBehavior::Component nC5  {"nC5" ,  489.38,           385.61 + 459.67,       0.0675,          72.151,              0.2510,    0};
    PhaseBehavior::Component N2   {"N2",    492.32,          -232.51 + 459.67,       0.0514,          28.014,              0.0400,    0};
    PhaseBehavior::Component nC10 {"nC10",  305.68,          652.01 + 459.67,        0.0679,          142.285,             0.4900,    0};

    auto mixture = PhaseBehavior::Mixture(std::pair(C1, 0.2), std::pair(nC5, 0.2),  std::pair(N2, 0.2), std::pair(nC10, 0.4));

    mixture.interactionCoefficient("N2", "nC10", 0.2);
    mixture.interactionCoefficient("nC10", "N2", 0.2);

    auto [result, phases] = isothermalTwoPhaseFlash<PengRobinson>(mixture, 500, 160 + 459.67);
    std::cout << "Phase 1 Density: " << phases[0]->density() << std::endl;
    if (result == PhaseStabilityResult::Unstable){
        std::cout << "Phase 2 Density: " << phases[1]->density() << std::endl;

    }

    for (const auto& component: mixture){
        std::cout << component.pure().name() << " xi is "<< component.composition("liquid") << " and yi is " << component.composition("vapor") <<std::endl;
        std::cout << component.pure().name() << " ki is " << component.equilibriumCoefficient() << std::endl;
    }

}