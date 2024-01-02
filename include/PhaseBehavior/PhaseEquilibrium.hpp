#ifndef PHASE_EQUILIBRIUM
#define PHASE_EQUILIBRIUM

#include <algorithm>

#include "Utilities/Math.hpp"
#include "Mixture.hpp"

using NP_t = PhaseBehavior::Types::NumericalPrecision;

namespace PhaseBehavior::VaporLiquidEquilibrium {
    void rachfordVLE(Mixture& mixture, const NP_t& pressure, const NP_t& temperature){

        auto [min, max] = minmax_element( mixture.begin(), mixture.end(), [&pressure, &temperature](const auto& lhs, const auto& rhs) {
            return lhs.equilibriumCoefficient(pressure, temperature) < rhs.equilibriumCoefficient(pressure, temperature);
        } );

        auto Kimax = max->equilibriumCoefficient(pressure, temperature);
        auto Kimin = min->equilibriumCoefficient(pressure, temperature);


        auto objectiveFunction = [&mixture, &pressure, &temperature](const NP_t& vaporMolarFraction){
            NP_t sum = 0;
            for (const auto& mixComp : mixture){
                auto Ki = mixComp.equilibriumCoefficient(pressure, temperature);
                sum += mixComp.composition()*(Ki-1)/(1+vaporMolarFraction*(Ki-1));
            }
            return sum;
        };

        auto divergence = [&Kimin, &Kimax](const NP_t& vaporMolarFraction){
            return !((1.0/(1.0-Kimax))< vaporMolarFraction && vaporMolarFraction < (1.0/(1.0-Kimin)));
        };

        auto hasRoot = objectiveFunction(0) > 0 && objectiveFunction(1) < 0;
        if (hasRoot){
            NP_t initialGuess = (1.0/(1.0-Kimax) + 1.0/(1.0-Kimin))/2.0;
        
            auto possibleFng = Math::NewtonRaphson<NP_t>(initialGuess,objectiveFunction,divergence);

            NP_t vaporMolarFraction;
            if (possibleFng){
                vaporMolarFraction = possibleFng.value();
            }else{
                vaporMolarFraction = Math::bisection<NP_t>(1.0/(1.0-Kimax), 1.0/(1.0-Kimin), objectiveFunction);
            }

            auto liquidMolarFraction = 1 - vaporMolarFraction;

            mixture.molarFraction("vapor", vaporMolarFraction);
            mixture.molarFraction("liquid", liquidMolarFraction);

            for (auto& mixComp : mixture){
                auto Ki = mixComp.equilibriumCoefficient(pressure, temperature);
                auto liquidComposition = mixComp.composition()/(1.0 + vaporMolarFraction*(Ki-1));
                auto vaporComposition  = Ki*mixComp.composition()/(1.0 + vaporMolarFraction*(Ki-1));
                mixComp.composition("liquid", liquidComposition);
                mixComp.composition("vapor", vaporComposition);
            }
        }
    }
}

#endif /* PHASE_EQUILIBRIUM */