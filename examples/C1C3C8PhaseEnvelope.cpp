#include <iostream>
#include <fstream>
#include <chrono>

#include <PhaseBehavior/Utilities/Input.hpp>
#include <PhaseBehavior/Component.hpp>
#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>
#include <PhaseBehavior/PhaseEnvelope.hpp>



int main(){
    auto C1 = PhaseBehavior::Component {"CH4",      4599.2,    190.56,   6.1477929425e-3, 16.04 ,  0.0114};
    auto C3 = PhaseBehavior::Component {"C3H8",     4251.2,    369.89,   4.5355587810e-3, 44.097,  0.1521};
    auto C8 = PhaseBehavior::Component {"n-C8H18",  2497.0,    569.32,   4.2571306948e-3, 114.23,  0.3950};

    PhaseBehavior::Mixture mixture{{C1, 1.0/3.0}, {C3, 1.0/3.0}, {C8, 1.0/3.0}};

    auto start = std::chrono::system_clock::now();
    auto envelope = PhaseBehavior::PhaseEnvelope<PhaseBehavior::EoS::PR::PengRobinson>(mixture);
    auto points = envelope.bruteForce(100, 1000.0, 101.3529, 250);
    auto end = std::chrono::system_clock::now();
    std::ofstream envelopeFile;
    envelopeFile.open("C1C3C8-envelope.csv");
    envelopeFile << "Pressure,BubbleTemperature,DewTemperature"<<std::endl;
    for (auto [pressure, bubbleTemp, dewTemp] : points){
        envelopeFile << pressure << "," << bubbleTemp << "," << dewTemp << std::endl;
    }
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    envelopeFile.close();
    std::cout << "File generated in " << elapsed.count() << " seconds" << std::endl;
    return 0;
}