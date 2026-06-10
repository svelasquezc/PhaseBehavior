#include <iostream>
#include <fstream>
#include <chrono>

#include <PhaseBehavior/Utilities/Input.hpp>
#include <PhaseBehavior/Component.hpp>
#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>
#include <PhaseBehavior/PhaseEnvelope.hpp>



int main(){
    using PhaseBehavior::Component;

    auto C1 = Component {"C1",      4599.2,    190.56,   6.1477929425e-3, 16.04,  0.0114};
    auto C2 = Component {"C2",      4872.2,    305.32,   4.8402710552e-3, 30.07,  0.0990};
    auto C3 = Component {"C3",      4251.2,    369.89,   4.5355587810e-3, 44.10, 0.1521};
    auto C4 = Component {"C4",      3796.3,    425.34,   4.3890449440e-3, 58.12, 0.1995};
    auto C5 = Component {"C5",      3368.8,    469.89,   4.3103448281e-3, 72.15,  0.2514};
    auto C6 = Component {"C6",      3012.3,    507.56,   4.2885324643e-3, 86.18, 0.2994};

    PhaseBehavior::Mixture mixture{{C1, 1.0/6.0}, {C2, 1.0/6.0}, {C3, 1.0/6.0}, {C4, 1.0/6.0}, {C5, 1.0/6.0}, {C6, 1.0/6.0}};

    auto start = std::chrono::system_clock::now();
    auto envelope = PhaseBehavior::PhaseEnvelope<PhaseBehavior::EoS::PR::PengRobinson>(mixture);
    auto points = envelope.bruteForce(100, 1000.0, 101.3529, 250);
    auto end = std::chrono::system_clock::now();
    std::ofstream envelopeFile;
    envelopeFile.open("C1toC6-envelope.csv");
    envelopeFile << "Pressure,BubbleTemperature,DewTemperature"<<std::endl;
    for (auto [pressure, bubbleTemp, dewTemp] : points){
        envelopeFile << pressure << "," << bubbleTemp << "," << dewTemp << std::endl;
    }
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    envelopeFile.close();
    std::cout << "File generated in " << elapsed.count() << " seconds" << std::endl;
    return 0;
}