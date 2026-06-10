#include <iostream>
#include <fstream>
#include <chrono>

#include <PhaseBehavior/Utilities/Input.hpp>
#include <PhaseBehavior/Component.hpp>
#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>
#include <PhaseBehavior/PhaseEnvelope.hpp>



int main(){
    auto C1 = PhaseBehavior::Component {"CH4",      666.4,     343.33,   0.0988, 16.043,  0.0104};
    auto C7plus = PhaseBehavior::Component {"C7+",  305.2,     1112,     0.0679, 142.285,  0.4898};

    PhaseBehavior::Mixture mixture{{C1, 0.5}, {C7plus, 0.5}};

    auto start = std::chrono::system_clock::now();
    auto envelope = PhaseBehavior::PhaseEnvelope<PhaseBehavior::EoS::PR::PengRobinson>(mixture);
    auto points = envelope.bruteForce(100, 3000);
    auto end = std::chrono::system_clock::now();
    std::ofstream envelopeFile;
    envelopeFile.open("CH4-C7+-envelope.csv");
    envelopeFile << "Pressure,BubbleTemperature,DewTemperature"<<std::endl;
    for (auto [pressure, bubbleTemp, dewTemp] : points){
        envelopeFile << pressure << "," << bubbleTemp << "," << dewTemp << std::endl;
    }
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    envelopeFile.close();
    std::cout << "File generated in " << elapsed.count() << " seconds" << std::endl;
    return 0;
}