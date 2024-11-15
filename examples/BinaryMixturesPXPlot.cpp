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

    double temperature = 344.261; // [K]
    double initialPressure = 101.325; // [kPa] 
    double finalPressure = 150000.0; // [kPa]
    auto deltaPressure = (finalPressure - initialPressure)/1000.0;
    auto pressure = initialPressure;

                                    /*[kPa]*/  /*[K]*/    /*  [m3/kg]  */ /*[g/mol]*/ /*[-]*/
    Component CO2      = {"CO2",      7380,    304.1282,  2.1385799828e-3, 44.0095,   0.239};                                
    Component CH4      = {"CH4",      4599,    190.564,   6.1477929425e-3, 16.04246,  0.011};
    Component N2       = { "N2",      3390,    126.192,   3.1918384511e-3, 28.0134,   0.039};
    Component n_C10H22 = {"n-C10H22", 2110,    617.7,     4.2855559307e-3, 142.28168, 0.4884};


    auto binaryCO2 = PhaseBehavior::Mixture(std::pair(CO2, 0.5), std::pair(n_C10H22, 0.5));
    auto binaryCH4 = PhaseBehavior::Mixture(std::pair(CH4, 0.5), std::pair(n_C10H22, 0.5));
    auto binaryN2  = PhaseBehavior::Mixture(std::pair(N2, 0.5),  std::pair(n_C10H22, 0.5));
    

    std::ofstream binaryCO2File, binaryN2File, binaryCH4File;
    binaryCO2File.open("binaryCO2.csv");
    binaryCH4File.open("binaryCH4.csv");
    binaryN2File.open("binaryN2.csv");

    binaryCO2File << "pressure;xi;yi"<<std::endl;
    binaryCH4File << "pressure;xi;yi"<<std::endl;
    binaryN2File << "pressure;xi;yi"<<std::endl;

    while(pressure < finalPressure){

        {
            auto result = isothermalTwoPhaseFlash<PengRobinson>(binaryCO2, pressure, temperature);
            if (result == PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable){
                binaryCO2File << pressure*0.145038 << ";" << binaryCO2[0].composition("liquid") << ";" << binaryCO2[0].composition("vapor") << std::endl;
                auto newZ = (binaryCO2[0].composition("liquid") + binaryCO2[0].composition("vapor"))/2;
                binaryCO2[0].composition(newZ);
                binaryCO2[1].composition(1 - newZ);

            }
        }
        
        {
            auto result = isothermalTwoPhaseFlash<PengRobinson>(binaryCH4, pressure, temperature);
            if (result == PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable){
                binaryCH4File << pressure*0.145038 << ";" << binaryCH4[0].composition("liquid") << ";" << binaryCH4[0].composition("vapor") << std::endl;
                auto newZ = (binaryCH4[0].composition("liquid") + binaryCH4[0].composition("vapor"))/2;
                binaryCH4[0].composition(newZ);
                binaryCH4[1].composition(1 - newZ);
            }
            
        }
        
        {
            auto result = isothermalTwoPhaseFlash<PengRobinson>(binaryN2,  pressure, temperature);
            if (result == PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable){
                binaryN2File << pressure*0.145038 << ";" << binaryN2[0].composition("liquid") << ";" << binaryN2[0].composition("vapor") << std::endl;
                auto newZ = (binaryCH4[0].composition("liquid") + binaryCH4[0].composition("vapor"))/2;
                binaryN2[0].composition(newZ);
                binaryN2[1].composition(1 - newZ);
            }
        }
        
        pressure += deltaPressure;
    }

    binaryCO2File.close();
    binaryCH4File.close();
    binaryN2File.close();
    

    return 0;
}