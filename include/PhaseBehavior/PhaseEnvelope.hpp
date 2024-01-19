#ifndef PHASE_ENVELOPE_HPP
#define PHASE_ENVELOPE_HPP

#include "Utilities/Types.hpp"
#include "Utilities/Math.hpp"
#include "Mixture.hpp"

using NP_t = PhaseBehavior::Types::NumericalPrecision;

namespace PhaseBehavior{
    template<typename EoS>
    NP_t phaseEnvelope(Mixture& mixture, NP_t const& vaporMolarFraction){

        EoS eos;
        NP_t standardConditionsPressure = 14.7; //psia
        NP_t standardConditionsTemperature = 77 + 460; //Rankine
        
        auto objectiveFunction = [&mixture, &vaporMolarFraction, pressure=standardConditionsPressure](NP_t const& temperature){
            NP_t sum = 0;
            mixture.initializeEquilibriumCoefficients(pressure, temperature); //Wilson Coeff
            for (const auto& mixComp : mixture){
                const auto Ki = mixComp.equilibriumCoefficient();
                sum += mixComp.composition()*(Ki-1.0)/(1.0+vaporMolarFraction*(Ki-1.0));
            }
            return sum;
        };

        auto convergenceResult = Math::NewtonRaphson<NP_t>(standardConditionsTemperature, objectiveFunction);

        NP_t initialTemperature = -460;
        if (convergenceResult){
            initialTemperature = convergenceResult.value();
        }

        

        return initialTemperature;
    }
}

#endif /* PHASE_ENVELOPE_HPP */