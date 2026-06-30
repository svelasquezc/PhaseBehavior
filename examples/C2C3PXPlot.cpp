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

    double temperature = 274.79; // [K]
    double initialPressure = 101.325; // [kPa] 
    double finalPressure = 150000.0; // [kPa]
    auto deltaPressure = (finalPressure - initialPressure)/1000.0;
    auto pressure = initialPressure;

                                    /*[kPa]*/  /*[K]*/    /*  [m3/kg]  */ /*[g/mol]*/ /*[-]*/
    auto C2 = Component {"C2",        4872.2,    305.32,   4.8402710552e-3, 30.07,  0.0990};
    auto C3 = Component {"C3",        4251.3,    369.89,   4.5355587810e-3, 44.097, 0.1521};


    auto binaryC2 = PhaseBehavior::Mixture(std::pair(C2, 0.5), std::pair(C3, 0.5));
    

    std::ofstream binaryC2File;
    binaryC2File.open("binaryC2C3.csv");

    binaryC2File << "pressure;xi;yi"<<std::endl;

    while(pressure < finalPressure){

        {
            auto result = isothermalTwoPhaseFlash<PengRobinson>(binaryC2, pressure, temperature);
            if (result == PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable){
                binaryC2File << pressure << ";" << binaryC2[0].composition(PhaseBehavior::PhaseName::liquid) << ";" << binaryC2[0].composition(PhaseBehavior::PhaseName::vapor) << std::endl;
                auto newZ = (binaryC2[0].composition(PhaseBehavior::PhaseName::liquid) + binaryC2[0].composition(PhaseBehavior::PhaseName::vapor))/2;
                binaryC2[0].composition(newZ);
                binaryC2[1].composition(1 - newZ);

            }
        }
        
        pressure += deltaPressure;
    }

    binaryC2File.close();

    return 0;
}