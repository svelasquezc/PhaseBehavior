#include <iostream>
#include <fstream>
#include <tuple>

#include <PhaseBehavior/Component.hpp>
#include <PhaseBehavior/Utilities/Constants.hpp>
#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>
#include <PhaseBehavior/PhaseEquilibrium.hpp>

using PhaseBehavior::Component, PhaseBehavior::Mixture, PhaseBehavior::EoS::PR::PengRobinson;
using PhaseBehavior::Constants::SI::universalGasesConstant;
using PhaseBehavior::VaporLiquidEquilibrium::isothermalTwoPhaseFlash;

int main(){
    PengRobinson prEoS;

    double pressure = 10000 ; //[kPa]
    double temperature = 350; // [K]

                                    /*[kPa]*/  /*[K]*/    /*  [m3/kg]  */ /*[g/mol]*/ /*[-]*/
    Component CH4      = {"CH4",      4599,    190.564,   6.1477929425e-3, 16.04246,  0.011};
    Component nC5      = {"n-C5H12",  3369,    469.7,     4.3103448281e-3, 72.14878,  0.251};

    Mixture mix = {{CH4, 0.99}, {nC5, 0.01}};
    mix.interactionCoefficient("CH4", "n-C5H12", 0.0236);
    mix.interactionCoefficient("n-C5H12", "CH4", 0.0236);

    auto compnC5 = 0.01;
    auto compCH4 = 1 - compnC5;
    auto deltaComp = (1-0.01)/100;

    auto RT = universalGasesConstant*temperature;

    auto ZCH4 = prEoS(CH4, pressure, temperature);
    auto ZnC5 = prEoS(nC5, pressure, temperature);

    auto fCH4 = prEoS.fugacity(CH4, ZCH4, pressure, temperature);
    std::cout << "fugacity CH4: "<< fCH4 <<std::endl;
    auto fnC5 = prEoS.fugacity(nC5, ZnC5, pressure, temperature);

    std::cout << "fugacity nC5: "<< fnC5 <<std::endl;

    auto pureGibbsCH4 = RT*std::log(fCH4);
    auto pureGibbsnC5 = RT*std::log(fnC5);

    std::ofstream molarGFile, deltaMolarGMixFile;

    molarGFile.open("molarG.csv");
    deltaMolarGMixFile.open("DeltaG.csv");

    molarGFile <<"xnC5;molarG"<<std::endl;
    deltaMolarGMixFile <<"xnC5;deltaGMix"<<std::endl;

    std::cout << mix[0].pure().name() <<std::endl;

    while (compnC5<1){
        prEoS(mix, pressure, temperature);
        mix.compressibility("global", prEoS.selectedCompressibility());
        prEoS.fugacities(mix, "global");

        auto DeltaMolarGMix = RT*(mix[0].composition()*std::log(mix[0].fugacity("global",pressure)/fCH4) + mix[1].composition()*std::log(mix[1].fugacity("global",pressure)/fnC5));
        deltaMolarGMixFile <<compnC5<<";"<< DeltaMolarGMix<<std::endl;

        auto MolarG = mix[0].composition()*pureGibbsCH4 + mix[1].composition()*pureGibbsnC5 + DeltaMolarGMix;
        molarGFile <<compnC5<<";"<< MolarG<<std::endl;

        compnC5 += deltaComp;
        compCH4 = 1 - compnC5;
        mix[0].composition(compCH4);
        mix[1].composition(compnC5);
    }
    molarGFile.close();
    deltaMolarGMixFile.close();
    
    mix = {{CH4, 0.5}, {nC5, 0.5}};
    mix.interactionCoefficient("CH4", "n-C5H12", 0.0236);
    mix.interactionCoefficient("n-C5H12", "CH4", 0.0236);
    auto result = isothermalTwoPhaseFlash<PengRobinson>(mix, pressure, temperature);
    if (result == PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable){
        std::cout <<"x_i: " << mix[1].composition("liquid") <<", y_i: "<< mix[1].composition("vapor") <<std::endl;
    }

}