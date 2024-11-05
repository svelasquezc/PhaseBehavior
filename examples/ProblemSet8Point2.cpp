#include <tuple>
#include <iostream>
#include <PhaseBehavior/Component.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>
#include <PhaseBehavior/Mixture.hpp>

int main(){

    PhaseBehavior::EoS::PR::PengRobinson eos;
                                  /*name*/  /*Pc (psia)*/  /*Tc (R)*/  /*Vc (ft3/lbm)*/  /*Mw (lbm/lbmol)*/ /*accentricF*/
    PhaseBehavior::Component CO2  {"CO2" ,  1070.3785,     547.524,      0.03426,          44.0095,          0.239};
    PhaseBehavior::Component nC4  {"nC4" ,  550.6,         765.62,       0.1995,           58.123,           0.0703};
    PhaseBehavior::Component nC7  {"nC7" ,  400,           972.27,       0.06843,          100.20,           0.349};
    PhaseBehavior::Component nC10 {"nC10",  304.21,        1111.77,      0.06865,          142.2817,         0.4890};

    auto mixture = PhaseBehavior::Mixture(std::pair(CO2, 0.9),std::pair(nC4, 0.05), std::pair(nC7, 0.03), std::pair(nC10, 0.02));

    mixture.interactionCoefficient(CO2.name(), nC4.name(), 0.1471);
    mixture.interactionCoefficient(nC4.name(), CO2.name(), 0.1471);
    mixture.interactionCoefficient(CO2.name(), nC7.name(), 0.1136);
    mixture.interactionCoefficient(nC7.name(), CO2.name(), 0.1136);
    mixture.interactionCoefficient(CO2.name(), nC10.name(), 0.1377);
    mixture.interactionCoefficient(nC10.name(),CO2.name(), 0.1377);

    eos(mixture, 1500 /*psia*/, 100 + 459.67 /* R */);
    mixture.compressibility("global", eos.selectedCompressibility());
    eos.fugacities(mixture, "global");

    for (const auto& component: mixture){
        std::cout << component.pure().name() << " fugacity is "<< component.fugacity("global", 1500.0) << " psia" <<std::endl;
    }

}