#ifndef PHASE_EQUILIBRIUM
#define PHASE_EQUILIBRIUM

#include <algorithm>
#include <iostream>

#include "Utilities/Math.hpp"
#include "Mixture.hpp"

using NP_t = PhaseBehavior::Types::NumericalPrecision;

namespace PhaseBehavior::VaporLiquidEquilibrium {

    void rachfordVLE(Mixture& mixture){

        auto [min, max] = std::minmax_element( mixture.begin(), mixture.end(), [](const auto& lhs, const auto& rhs) {
            return lhs.equilibriumCoefficient() < rhs.equilibriumCoefficient();
        } );

        auto Kimax = max->equilibriumCoefficient();
        auto Kimin = min->equilibriumCoefficient();

        auto lowerLimit = 1.0/(1.0-Kimax + 1e-10);
        auto upperLimit = 1.0/(1.0-Kimin - 1e-10);

        auto objectiveFunction = [&mixture](const NP_t& vaporMolarFraction){
            NP_t sum = 0;
            for (const auto& mixComp : mixture){
                const auto Ki = mixComp.equilibriumCoefficient();
                sum += mixComp.composition()*(Ki-1.0)/(1.0+vaporMolarFraction*(Ki-1.0));
            }
            return sum;
        };

        auto derivativeFunction = [&mixture](const NP_t& vaporMolarFraction){
            NP_t sum = 0;
            for (const auto& mixComp : mixture){
                const auto Ki = mixComp.equilibriumCoefficient();
                sum -= mixComp.composition()*std::pow(Ki-1.0,2)/std::pow(1.0+vaporMolarFraction*(Ki-1.0),2);
            }
            return sum;
        };

        auto divergence = [&lowerLimit, &upperLimit](const NP_t& vaporMolarFraction){
            return !(lowerLimit < vaporMolarFraction && vaporMolarFraction < upperLimit);
        };

        auto hasRoot = objectiveFunction(0) > 0 && objectiveFunction(1) < 0;
        if (hasRoot){
            NP_t initialGuess = (lowerLimit + upperLimit)/2.0;
        
            auto possibleFng = Math::NewtonRaphson<NP_t>(initialGuess,divergence,objectiveFunction, derivativeFunction);

            NP_t vaporMolarFraction;
            if (possibleFng){
                vaporMolarFraction = possibleFng.value();
            }else{
                vaporMolarFraction = Math::bisection<NP_t>(0.0, 1.0, objectiveFunction);
            }

            auto liquidMolarFraction = 1 - vaporMolarFraction;

            mixture.molarFraction("vapor", vaporMolarFraction);
            mixture.molarFraction("liquid", liquidMolarFraction);

            for (auto& mixComp : mixture){
                auto Ki = mixComp.equilibriumCoefficient();
                auto liquidComposition = mixComp.composition()/(1.0 + vaporMolarFraction*(Ki-1));
                auto vaporComposition  = Ki*liquidComposition;
                mixComp.composition("liquid", liquidComposition);
                mixComp.composition("vapor", vaporComposition);
            }
        }
    }

    template<typename EoS>
    void SuccesiveSubstitution(Mixture& mixture, NP_t const& pressure, NP_t const& temperature){

        mixture.initializeEquilibriumCoefficients(pressure, temperature); //Wilson Correlation

        EoS vaporEoS, liquidEoS;

        NP_t fugacitiesSquaredSum = 10;
        while(fugacitiesSquaredSum >= 1e-14){

            rachfordVLE(mixture);

            vaporEoS(mixture, pressure, temperature, "vapor");
            mixture.compressibility("vapor", vaporEoS.selectedCompressibility());
            vaporEoS.fugacities(mixture, "vapor");

            liquidEoS(mixture, pressure, temperature, "liquid");
            mixture.compressibility("liquid", liquidEoS.selectedCompressibility());
            liquidEoS.fugacities(mixture, "liquid");

            for (auto& mixComponent : mixture){
                auto newEquilibriumCoefficient = mixComponent.equilibriumCoefficient()*(mixComponent.fugacity("liquid", pressure)/mixComponent.fugacity("vapor", pressure));
                mixComponent.equilibriumCoefficient(newEquilibriumCoefficient);
            }

            fugacitiesSquaredSum = std::accumulate(mixture.begin(),mixture.end(), static_cast<NP_t>(0.0),
            [&pressure](auto previous, auto second){
                return previous + std::pow(second.fugacity("liquid", pressure)/second.fugacity("vapor", pressure) - 1.0, 2);
            });
        }
    }
}

#endif /* PHASE_EQUILIBRIUM */