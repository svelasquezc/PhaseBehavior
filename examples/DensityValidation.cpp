#include <iostream>
#include <tuple>
#include <fstream>
#include <sstream>
#include <thread>

#include <PhaseBehavior/Component.hpp>
#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/Phase.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>

using PR = PhaseBehavior::EoS::PR::PengRobinson;
using SRK = PhaseBehavior::EoS::SRK::SoaveRedlichKwong;


template<typename EoS>
void densityCalculation(PhaseBehavior::Component component1, PhaseBehavior::Component component2, double composition1, double composition2, double temperature){
    EoS eos;

    auto mixture = PhaseBehavior::Mixture(std::pair(component1, composition1), std::pair(component2, composition2));
    
    std::stringstream iss;

    iss << "Density"<<"-"<< eos.name() <<"-"<< composition1*100 << "%" << component1.name() << "+" << composition2*100 << "%" << component2.name() << "-" << static_cast<int>(temperature/1.8) <<"K"<< ".csv";

    std::ofstream file(iss.str());

    file << "Pressure;Density"<<std::endl;

    for (int pressure=0; pressure<= 2901; pressure+=5){
        eos(mixture, static_cast<double>(pressure), temperature);
        mixture.compressibility("global", eos.selectedCompressibility());

        auto gas = PhaseBehavior::Phase::VaporLikePhase(mixture, "global");
        gas.molarVolume(mixture.compressibility("global"), pressure, temperature);
        auto density = gas.density();
        file << pressure*0.00689476 << ";" << density*16.01846 << std::endl;
    }
    file.close();
}

int main(){
    PhaseBehavior::Component H2  {"H2",  188.549059,59.724,194,2.01588,-0.216};
    PhaseBehavior::Component CH4 {"CH4", 668.624, 343.08, 0.0985, 16.0425, 0.011};
    PhaseBehavior::Component N2  {"N2",  492.8092253, 227.142, 0.051, 28.0134, 0.039};
    PhaseBehavior::Component CO2 {"CO2", 1070.3785, 547.524, 0.03426, 44.0095, 0.239};

    std::thread t1(densityCalculation<PR>, H2, CH4, 0.5, 0.5, 585);
    std::thread t2(densityCalculation<SRK>, H2, CH4, 0.5, 0.5, 585);
    std::thread t3(densityCalculation<PR>, H2, CH4, 0.5, 0.5, 630);
    std::thread t4(densityCalculation<SRK>, H2, CH4, 0.5, 0.5, 630);

    std::thread t5(densityCalculation<PR>, H2, CH4, 0.05, 0.95, 540);
    std::thread t6(densityCalculation<SRK>, H2, CH4, 0.05, 0.95, 540);
    std::thread t7(densityCalculation<PR>, H2, CH4, 0.05, 0.95, 585);
    std::thread t8(densityCalculation<SRK>, H2, CH4, 0.05, 0.95, 585);

    std::thread t11(densityCalculation<PR>, H2, N2, 0.5, 0.5, 585);
    std::thread t12(densityCalculation<SRK>, H2, N2, 0.5, 0.5, 585);
    std::thread t13(densityCalculation<PR>, H2, N2, 0.5, 0.5, 630);
    std::thread t14(densityCalculation<SRK>, H2, N2, 0.5, 0.5, 630);

    std::thread t15(densityCalculation<PR>, H2, CO2, 0.1, 0.9, 545.4);
    std::thread t16(densityCalculation<SRK>, H2, CO2, 0.1, 0.9, 545.4);
    std::thread t17(densityCalculation<PR>, H2, CO2, 0.1, 0.9, 581.4);
    std::thread t18(densityCalculation<SRK>, H2, CO2, 0.1, 0.9, 581.4);

    t1.join();
    t2.join();
    t3.join();
    t4.join();

    t5.join();
    t6.join();
    t7.join();
    t8.join();

    t11.join();
    t12.join();
    t13.join();
    t14.join();

    t15.join();
    t16.join();
    t17.join();
    t18.join();

    return 0;
}