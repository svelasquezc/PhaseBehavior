#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>
#include <utility>

#include <PhaseBehavior/Utilities/Input.hpp>
#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>
#include <PhaseBehavior/PhaseEnvelope.hpp>

#include <rapidcsv.h>

// Worker implementation: Unrolls the creation and the joins.
// Now accepts an 'offset' so the lambda gets a true global index.
template <typename F, std::size_t... Is>
void run_unrolled_threads_impl(F&& f, std::size_t offset, std::index_sequence<Is...>) {
    // Pack expansion: Unroll thread creation. Pass the offset + local index.
    std::array<std::thread, sizeof...(Is)> threads{ 
        std::thread(f, offset + Is)... 
    };

    // Fold Expression: Unroll the .join() operations.
    (threads[Is].join(), ...);
}

// Batched public interface
template <std::size_t Total, std::size_t BatchSize, typename F>
void run_batched_unrolled_threads(F&& f) {
    // Math to figure out our loop and remainder
    constexpr std::size_t NumBatches = Total / BatchSize;
    constexpr std::size_t Remainder = Total % BatchSize;

    // 1. Execute the full batches using a standard runtime loop
    for (std::size_t b = 0; b < NumBatches; ++b) {
        std::size_t offset = b * BatchSize;
        run_unrolled_threads_impl(std::forward<F>(f), offset, std::make_index_sequence<BatchSize>{});
    }

    // 2. Execute the remainder, if any exist
    if constexpr (Remainder > 0) {
        std::size_t offset = NumBatches * BatchSize;
        run_unrolled_threads_impl(std::forward<F>(f), offset, std::make_index_sequence<Remainder>{});
    }
}


int main(){

    constexpr int numberOfComponents = 6;
    constexpr std::size_t NTHREADS = 16;

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

    
    rapidcsv::Document doc("SPE3_UpscaledProps_th90.csv");

    std::vector<int> time = doc.GetColumn<int>("Time (s)");
    std::vector<double> avgP = doc.GetColumn<double>("Average Pressure (Pa)");
    std::vector<double> xC1 = doc.GetColumn<double>("C1 Composition (1)");
    std::vector<double> xC2 = doc.GetColumn<double>("C2 Composition (1)");
    std::vector<double> xC3 = doc.GetColumn<double>("C3 Composition (1)");
    std::vector<double> xC4 = doc.GetColumn<double>("C4 Composition (1)");
    std::vector<double> xC5 = doc.GetColumn<double>("C5 Composition (1)");
    std::vector<double> xC6 = doc.GetColumn<double>("C6 Composition (1)");

    auto envelopeGen = [=](std::size_t pos){

        PhaseBehavior::Mixture mixture {{C1, xC1[pos]}, {C2, xC2[pos]}, {C3, xC3[pos]}, {n_C4, xC4[pos]}, {n_C5, xC5[pos]}, {n_C6, xC6[pos]}};        
        auto envelope = PhaseBehavior::PhaseEnvelope<PhaseBehavior::EoS::PR::PengRobinson>(mixture);
        auto points = envelope.bruteForce(100, 1000.0, 101.3529, 250);
        std::string envelopeFileName = "SPE3UpscaledPhaseEnvelope-t" + std::to_string(time[pos]) + ".csv";
        std::ofstream envelopeFile;
        envelopeFile.open(envelopeFileName);
        envelopeFile << "Pressure,BubbleTemperature,DewTemperature"<<std::endl;
        for (auto [pressure, bubbleTemp, dewTemp] : points){
            envelopeFile << pressure << "," << bubbleTemp << "," << dewTemp << std::endl;
        }
        envelopeFile.close();
    };

    auto start = std::chrono::system_clock::now();
    
    run_batched_unrolled_threads<78, NTHREADS>(envelopeGen);

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    std::cout << "Files generated in " << elapsed.count() << " seconds" << std::endl;

    return 0;
}