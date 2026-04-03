#include <iostream>
#include <fstream>
#include <chrono>

#include <PhaseBehavior/Utilities/Input.hpp>
#include <PhaseBehavior/Component.hpp>
#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>
#include <PhaseBehavior/PhaseEnvelope.hpp>



int main(){
    auto C3 = PhaseBehavior::Component {"C3",        4251.3,    369.89,   4.5355587810e-3, 44.097, 0.1521};
    auto C4 = PhaseBehavior::Component {"C4",        3796.3,    425.34,   4.3890449440e-3, 58.12,  0.1995};

    PhaseBehavior::Mixture mixture{{C3, 0.9}, {C4, 0.1}};

    auto start = std::chrono::system_clock::now();
    auto envelope = PhaseBehavior::PhaseEnvelope<PhaseBehavior::EoS::PR::PengRobinson>(mixture);
    auto points = envelope.bruteForce(100, 1000.0, 101.3529, 250);
    auto end = std::chrono::system_clock::now();
    std::ofstream envelopeFile;
    envelopeFile.open("C3C4-envelope.csv");
    envelopeFile << "Pressure,BubbleTemperature,DewTemperature"<<std::endl;
    for (auto [pressure, bubbleTemp, dewTemp] : points){
        envelopeFile << pressure << "," << bubbleTemp << "," << dewTemp << std::endl;
    }
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    envelopeFile.close();
    std::cout << "File generated in " << elapsed.count() << " seconds" << std::endl;
    return 0;
}