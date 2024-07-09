#include <iostream>
#include <tuple>
#include <fstream>
#include <sstream>
#include <thread>

#include <PhaseBehavior/Component.hpp>
#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/Phase.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>
#include <PhaseBehavior/GERG/GERG2008.hpp>
#include <PhaseBehavior/GERG/Components.hpp>

using PR = PhaseBehavior::EoS::PR::PengRobinson;
using SRK = PhaseBehavior::EoS::SRK::SoaveRedlichKwong;
using GERG2008 = PhaseBehavior::EoS::GERG::GERG2008;


template<typename EoS>
void densityCalculation(PhaseBehavior::Component component1, PhaseBehavior::Component component2, double composition1, double composition2, double temperature){
    EoS eos;

    auto mixture = PhaseBehavior::Mixture(std::pair(component1, composition1), std::pair(component2, composition2));
    
    std::stringstream iss;

    iss << "Density"<<"-"<< eos.name() <<"-"<< composition1*100 << "%" << component1.name() << "+" << composition2*100 << "%" << component2.name() << "-" << static_cast<int>(temperature) <<"K"<< ".csv";

    std::ofstream file(iss.str());

    file << "Pressure;Density"<<std::endl;

    for (int pressure=1; pressure <= 20000; pressure+=100){
        eos(mixture, static_cast<double>(pressure), temperature);
        mixture.compressibility("global", eos.selectedCompressibility());

        auto gas = PhaseBehavior::Phase::VaporLikePhase(mixture, "global");
        gas.molarVolume(mixture.compressibility("global"), pressure, temperature);
        auto density = gas.density();
        file << static_cast<double>(pressure)/1000.0 << ";" << density << std::endl;
    }
    file.close();
}

int main(){
    

    auto H2 = PhaseBehavior::EoS::GERG::Components::H2;
    auto CH4 = PhaseBehavior::EoS::GERG::Components::CH4;
    auto N2 = PhaseBehavior::EoS::GERG::Components::N2;
    auto CO2 = PhaseBehavior::EoS::GERG::Components::CO2;

    std::thread t1(densityCalculation<PR>, H2, CH4, 0.5, 0.5, 325);
    std::thread t2(densityCalculation<SRK>, H2, CH4, 0.5, 0.5, 325);
    std::thread t3(densityCalculation<PR>, H2, CH4, 0.5, 0.5, 350);
    std::thread t4(densityCalculation<SRK>, H2, CH4, 0.5, 0.5, 350);

    std::thread t5(densityCalculation<PR>, H2, CH4, 0.05, 0.95, 300);
    std::thread t6(densityCalculation<SRK>, H2, CH4, 0.05, 0.95, 300);
    std::thread t7(densityCalculation<PR>, H2, CH4, 0.05, 0.95, 325);
    std::thread t8(densityCalculation<SRK>, H2, CH4, 0.05, 0.95, 325);

    std::thread t11(densityCalculation<PR>, H2, N2, 0.5, 0.5, 325);
    std::thread t12(densityCalculation<SRK>, H2, N2, 0.5, 0.5, 325);
    std::thread t13(densityCalculation<PR>, H2, N2, 0.5, 0.5, 350);
    std::thread t14(densityCalculation<SRK>, H2, N2, 0.5, 0.5, 350);

    std::thread t15(densityCalculation<PR>, H2, CO2, 0.1, 0.9, 303);
    std::thread t16(densityCalculation<SRK>, H2, CO2, 0.1, 0.9, 303);
    std::thread t17(densityCalculation<PR>, H2, CO2, 0.1, 0.9, 323);
    std::thread t18(densityCalculation<SRK>, H2, CO2, 0.1, 0.9, 323);

    std::thread t9(densityCalculation<GERG2008>, H2, CH4, 0.5, 0.5, 325);
    std::thread t10(densityCalculation<GERG2008>, H2, CH4, 0.5, 0.5, 350);
    std::thread t19(densityCalculation<GERG2008>, H2, CH4, 0.05, 0.95, 300);
    std::thread t20(densityCalculation<GERG2008>, H2, CH4, 0.05, 0.95, 325);
    std::thread t21(densityCalculation<GERG2008>, H2, N2, 0.5, 0.5, 325);
    std::thread t22(densityCalculation<GERG2008>, H2, N2, 0.5, 0.5, 350);
    std::thread t23(densityCalculation<GERG2008>, H2, CO2, 0.1, 0.9, 303);
    std::thread t24(densityCalculation<GERG2008>, H2, CO2, 0.1, 0.9, 323);


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

    t9.join();
    t10.join();
    t19.join();
    t20.join();
    t21.join();
    t22.join();
    t23.join();
    t24.join();

    return 0;
}