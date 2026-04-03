#include <iostream>
#include <fstream>
#include <chrono>

#include <PhaseBehavior/Utilities/Input.hpp>
#include <PhaseBehavior/Component.hpp>
#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>
#include <PhaseBehavior/PhaseEnvelope.hpp>



int main(){
    auto C2 = PhaseBehavior::Component {"C2",        4872.2,    305.32,   4.8402710552e-3, 30.07,  0.0990};
    auto C3 = PhaseBehavior::Component {"C3",        4251.3,    369.89,   4.5355587810e-3, 44.097, 0.1521};

    PhaseBehavior::Mixture mixture{{C2, 0.1}, {C3, 0.9}};

    auto start = std::chrono::system_clock::now();
    auto envelope = PhaseBehavior::PhaseEnvelope<PhaseBehavior::EoS::PR::PengRobinson>(mixture);
    auto points = envelope.bruteForce(100, 1000.0, 101.3529, 150);
    auto end = std::chrono::system_clock::now();
    std::ofstream envelopeFile;
    envelopeFile.open("C2C3-envelope.csv");
    envelopeFile << "Pressure,BubbleTemperature,DewTemperature"<<std::endl;
    for (auto [pressure, bubbleTemp, dewTemp] : points){
        envelopeFile << pressure << "," << bubbleTemp << "," << dewTemp << std::endl;
    }
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    envelopeFile.close();
    std::cout << "File generated in " << elapsed.count() << " seconds" << std::endl;
    return 0;
}