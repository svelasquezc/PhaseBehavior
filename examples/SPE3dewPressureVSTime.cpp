#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>
#include <utility>
#include <functional>

#include <PhaseBehavior/Utilities/Input.hpp>
#include <PhaseBehavior/Utilities/Math.hpp>
#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>
#include <PhaseBehavior/PhaseEnvelope.hpp>

#include <rapidcsv.h>

int main(){

    constexpr double temperature = 350.0; // [K] Example temperature in Kelvin

    auto C1 = PhaseBehavior::Component {"C1",      4599.2,    190.56,   6.1477929425e-3, 16.04,  0.0114};
    auto C2 = PhaseBehavior::Component {"C2",      4872.2,    305.32,   4.8402710552e-3, 30.07,  0.0990};
    auto C3 = PhaseBehavior::Component {"C3",      4251.2,    369.89,   4.5355587810e-3, 44.10, 0.1521};
    auto C4 = PhaseBehavior::Component {"C4",      3796.3,    425.34,   4.3890449440e-3, 58.12, 0.1995};
    auto C5 = PhaseBehavior::Component {"C5",      3368.8,    469.89,   4.3103448281e-3, 72.15,  0.2514};
    auto C6 = PhaseBehavior::Component {"C6",      3012.3,    507.56,   4.2885324643e-3, 86.18, 0.2994};

    std::ofstream outputFile("dew_pressure_results.csv");
    outputFile << "Time,DewPressure"<<std::endl;

    rapidcsv::Document doc("SPE3_UpscaledProps_th90.csv");

    std::vector<int> time = doc.GetColumn<int>("Time (s)");
    std::vector<double> avgP = doc.GetColumn<double>("Average Pressure (Pa)");
    std::vector<double> xC1 = doc.GetColumn<double>("C1 Composition (1)");
    std::vector<double> xC2 = doc.GetColumn<double>("C2 Composition (1)");
    std::vector<double> xC3 = doc.GetColumn<double>("C3 Composition (1)");
    std::vector<double> xC4 = doc.GetColumn<double>("C4 Composition (1)");
    std::vector<double> xC5 = doc.GetColumn<double>("C5 Composition (1)");
    std::vector<double> xC6 = doc.GetColumn<double>("C6 Composition (1)");


    auto objectiveDewPoint = [vaporMolarFraction=1, temperature=temperature](PhaseBehavior::Mixture& mixture, NP_t const& pressure){
                NP_t sum = 0;
                mixture.initializeEquilibriumCoefficients(pressure, temperature); //Wilson Coeff
                for (const auto& mixComp : mixture){
                    const auto Ki = mixComp.equilibriumCoefficient();
                    sum += mixComp.composition()*(Ki-1.0)/(1.0+vaporMolarFraction*(Ki-1.0));
                }
                return sum;
            };

    for(std::size_t i = 0; i < xC1.size(); ++i){
        PhaseBehavior::Mixture mixture {{C1, xC1[i]}, {C2, xC2[i]}, {C3, xC3[i]}, {C4, xC4[i]}, {C5, xC5[i]}, {C6, xC6[i]}};

        auto initialGuessPressure = 1.3e3; //[kPa]

        auto dewConvergence = PhaseBehavior::Math::NewtonRaphson(initialGuessPressure, static_cast<std::function<NP_t(NP_t)>>(std::bind(objectiveDewPoint, mixture, std::placeholders::_1)));

        NP_t dewPressure = 0;
        
        if (dewConvergence){
            dewPressure = dewConvergence.value();
        }else{
            throw std::runtime_error("Failed to converge on dew point");
        }
        outputFile << time[i] << "," << dewPressure << std::endl;
    }

    return 0;

}