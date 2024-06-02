#include <iostream>
#include <tuple>
#include <fstream>
#include <sstream>

#include <PhaseBehavior/Component.hpp>
#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/Phase.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>

using PR = PhaseBehavior::EoS::PR::PengRobinson;
using SRK = PhaseBehavior::EoS::SRK::SoaveRedlichKwong;


template<typename EoS>
void densityCalculation(PhaseBehavior::Component component1, PhaseBehavior::Component comp2, double composition1, double composition2, double temperature){
    EoS eos;

    auto mixture = PhaseBehavior::Mixture(std::pair(component1, composition1), std::pair(component2, composition2));
    
    for (int pressure=0; pressure<= 2900; pressure+=5){
        eos(mixture, static_cast<double>(pressure), temperature);
        mixture.compressibility["global"]
    }
}

int main(){
    PhaseBehavior::Component H2  {"H2",  188.549059,59.724,194,2.01588,-0.216};
    PhaseBehavior::Component CH4 {"CH4", 668.624, 343.08, 0.0985, 16.0425, 0.011};
    PhaseBehavior::Component N2  {"N2",  492.8092253, 227.142, 0.051, 28.0134, 0.039};
    PhaseBehavior::Component CO2 {"CO2", 1070.3785, 547.524, 0.03426, 44.0095, 0.239};

    return 0;
}