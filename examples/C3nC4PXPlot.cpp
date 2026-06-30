#include <iostream>
#include <fstream>
#include <tuple>

#include <PhaseBehavior/Utilities/Types.hpp>
#include <PhaseBehavior/Utilities/Input.hpp>
#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>
#include <PhaseBehavior/PhaseEquilibrium.hpp>

int main(){

    using PhaseBehavior::Component, PhaseBehavior::EoS::PR::PengRobinson, PhaseBehavior::VaporLiquidEquilibrium::isothermalTwoPhaseFlash;

    double temperature = 258; // [K]
    double initialPressure = 101.325; // [kPa] 
    double finalPressure = 6000.0; // [kPa]
    auto deltaPressure = (finalPressure - initialPressure)/1000.0;
    auto pressure = initialPressure;

                                    /*[kPa]*/  /*[K]*/    /*  [m3/kg]  */ /*[g/mol]*/ /*[-]*/
    auto C3 = Component {"C3",        4251.3,    369.89,   4.5355587810e-3, 44.097, 0.1521};
    auto C4 = Component {"C4",        3796.3,    425.34,   4.3890449440e-3, 58.12,  0.1995};


    auto binaryC3C4 = PhaseBehavior::Mixture(std::pair(C3, 0.8), std::pair(C4, 0.2));
    

    std::ofstream binaryC3C4File;
    binaryC3C4File.open("binaryC3C4.csv");

    binaryC3C4File << "pressure;xi;yi"<<std::endl;

    while(pressure < finalPressure){

        {
            auto result = isothermalTwoPhaseFlash<PengRobinson>(binaryC3C4, pressure, temperature);
            if (result == PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable){
                binaryC3C4File << pressure << ";" << binaryC3C4[0].composition(PhaseBehavior::PhaseName::liquid) << ";" << binaryC3C4[0].composition(PhaseBehavior::PhaseName::vapor) << std::endl;
                auto newZ = (binaryC3C4[0].composition(PhaseBehavior::PhaseName::liquid) + binaryC3C4[0].composition(PhaseBehavior::PhaseName::vapor))/2;
                binaryC3C4[0].composition(newZ);
                binaryC3C4[1].composition(1 - newZ);

            }
        }
        
        pressure += deltaPressure;
    }

    binaryC3C4File.close();

    return 0;
}