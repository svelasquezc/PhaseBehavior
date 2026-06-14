#include <iostream>
#include <fstream>
#include <chrono>

#include <PhaseBehavior/Utilities/Input.hpp>
#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>
#include <PhaseBehavior/PhaseEnvelope.hpp>



int main(){

    constexpr int numberOfComponents = 10;

    PhaseBehavior::Component CO2   = {"CO2",      7384.3,    304.39,   2.1385799828e-3,  44.01,   0.2667};
    PhaseBehavior::Component C1    = {"CH4",      4594.7,    190.74,   6.1477929425e-3,  16.043,  0.0104};
    PhaseBehavior::Component C2    = {"C2H6",     4871.1,    305.51,   4.8402710552e-3,  30.070,  0.0979};
    PhaseBehavior::Component C3    = {"C3H8",     4247.2,    370.03,   4.5355587810e-3,  44.097,  0.1522};
    PhaseBehavior::Component i_C4  = {"i-C4H10",  3639.7,    408.03,   4.4571224820e-3,  58.123,  0.1852};
    PhaseBehavior::Component n_C4  = {"n-C4H10",  3796.3,    425.34,   4.3890449440e-3,  58.123,  0.1995};
    PhaseBehavior::Component i_C5  = {"i-C5H12",  3381.2,    460.61,   4.2373122060e-3,  72.150,  0.2280};
    PhaseBehavior::Component n_C5  = {"n-C5H12",  3368.8,    469.89,   4.3103448281e-3,  72.150,  0.2514};
    PhaseBehavior::Component n_C6  = {"n-C6H14",  3012.3,    507.56,   4.2885324643e-3,  86.177,  0.2994};
    PhaseBehavior::Component C7    = {"C7+",      2104.3,    617.78,   4.3103448273e-3,  142.285, 0.4898};

    auto envelopeGen = [](auto envelopeName, auto mixture){
        auto start = std::chrono::system_clock::now();
        auto envelope = PhaseBehavior::PhaseEnvelope<PhaseBehavior::EoS::PR::PengRobinson>(mixture);
        auto points = envelope.bruteForce(100, 1000.0, 101.3529, 250);
        auto end = std::chrono::system_clock::now();
        std::ofstream envelopeFile;
        envelopeFile.open(envelopeName);
        envelopeFile << "Pressure,BubbleTemperature,DewTemperature"<<std::endl;
        for (auto [pressure, bubbleTemp, dewTemp] : points){
            envelopeFile << pressure << "," << bubbleTemp << "," << dewTemp << std::endl;
        }
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
        envelopeFile.close();
        std::cout << "File generated in " << elapsed.count() << " seconds" << std::endl;
    };

    if constexpr (numberOfComponents == 6){
        PhaseBehavior::Mixture mixture {{C1, 0.6223}, {C2, 0.1408}, {C3, 0.0835}, {n_C4, 0.0438}, {n_C5, 0.0232}, {n_C6, 0.0864}};
        envelopeGen("SPE3-envelope-6comp.csv", mixture);
    }else if constexpr (numberOfComponents == 8){
        PhaseBehavior::Mixture mixture {{C1, 0.6223}, {C2, 0.1408}, {C3, 0.0835}, {i_C4, 0.0097}, {n_C4, 0.0341}, {i_C5, 0.0084}, {n_C5, 0.0148}, {n_C6, 0.0864}};
        envelopeGen("SPE3-envelope-8comp.csv", mixture);
    }else if constexpr (numberOfComponents == 10){
        PhaseBehavior::Mixture mixture {{CO2, 0.0031}, {C1, 0.6192}, {C2, 0.1408}, {C3, 0.0835}, {i_C4, 0.0097}, {n_C4, 0.0341}, {i_C5, 0.0084}, {n_C5, 0.0148}, {n_C6, 0.0179}, {C7, 0.0685}};

        mixture.interactionCoefficient("CO2", "CH4", 0.105);
        mixture.interactionCoefficient("CH4", "CO2", 0.105);
        mixture.interactionCoefficient("CO2", "C2H6", 0.13);
        mixture.interactionCoefficient("C2H6", "CO2", 0.13);
        mixture.interactionCoefficient("CO2", "C3H8", 0.125);
        mixture.interactionCoefficient("C3H8", "CO2", 0.125);
        mixture.interactionCoefficient("CO2", "i-C4H10", 0.12);
        mixture.interactionCoefficient("i-C4H10", "CO2", 0.12);
        mixture.interactionCoefficient("CO2", "n-C4H10", 0.115);
        mixture.interactionCoefficient("n-C4H10", "CO2", 0.115);
        mixture.interactionCoefficient("CO2", "i-C5H12", 0.115);
        mixture.interactionCoefficient("i-C5H12", "CO2", 0.115);
        mixture.interactionCoefficient("CO2", "n-C5H12", 0.115);
        mixture.interactionCoefficient("n-C5H12", "CO2", 0.115);
        mixture.interactionCoefficient("CO2", "n-C6H14", 0.115);
        mixture.interactionCoefficient("n-C6H14", "CO2", 0.115);
        mixture.interactionCoefficient("CO2", "C7+", 0.115);
        mixture.interactionCoefficient("C7+", "CO2", 0.115);

        envelopeGen("SPE3-envelope.csv", mixture);
    }else{
        std::cout << "Invalid number of components selected. Please select either 6, 8 or 10 components." << std::endl;
    };
    return 0;
}