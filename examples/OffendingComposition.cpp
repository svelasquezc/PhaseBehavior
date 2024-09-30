#include <iostream>
#include <fstream>
#include <chrono>

#include <PhaseBehavior/Utilities/Types.hpp>
#include <PhaseBehavior/Utilities/Input.hpp>
#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/BlackOil/CVD.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>

using NP_t = PhaseBehavior::Types::NumericalPrecision;

using PR = PhaseBehavior::EoS::PR::PengRobinson;

struct PVTTableTest{
    NP_t P;
    NP_t Bo;
    NP_t Bg;
    NP_t Rs;
    NP_t Rv;
};


int main(){
    auto mixture = PhaseBehavior::Input::createMixtureFromFile<true>("PVTPR.csv", "InteractionCoefficients.csv");

    std::vector<NP_t> compVec = {5.7408778068051594e-05 , 0.0021097837089117978 , 0.012749626350470398 , 0.059403031413822083 , 0.020486151550549885 , 0.098657511602738612 , 0.04245825280734121 , 0.083270306130554667 , 0.12841091612898009 , 0.55239701152856324};

    for (std::size_t i = 0 ; i< mixture.size(); ++i){
        mixture[i].composition(compVec[i]);
    }

    auto result = PhaseBehavior::VaporLiquidEquilibrium::isothermalTwoPhaseFlash<PR>(mixture, 14.7, 60.0+460);
}