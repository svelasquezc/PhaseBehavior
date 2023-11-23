#ifndef PENG_ROBINSON_EOS
#define PENG_ROBINSON_EOS

#include <cmath>

#include "GeneralizedCubicEoS.hpp" 
#include "Utilities/Types.hpp"
#include "Component.hpp"

namespace PhaseBehavior::EoS::PR {
    using NP_t = TypesDefinition::NumericalPrecision;

    namespace constants {
        
        constexpr NP_t m1 () {return 1 + std::sqrt(2);}
        constexpr NP_t m2 () {return 1 - std::sqrt(2);}
        constexpr NP_t omegaA () {return 0.457235529;}
        constexpr NP_t omegaB () {return 0.077796074;}
    }

    class PengRobinsonEoS : public GeneralizedCubicEoS<
    constants::m1,
    constants::m2,
    constants::omegaA,
    constants::omegaB,
    PengRobinsonEoS
    >{
        NP_t mParameter(NP_t accentricFactor){
            if (accentricFactor <= 0.49){
                return 0.374640 + 1.54226*accentricFactor 
                - 0.26992*std::pow(accentricFactor,2);
            }
            return 0.379642 + 1.48503*accentricFactor 
            - 0.164423*std::pow(accentricFactor,2) + 0.016666*std::pow(accentricFactor,3);
        }

        NP_t componentAttraction(Component const& component, NP_t const& pressure, NP_t const& temperature){
            return constants::omegaA()*
            std::pow(1+ mParameter(component.accentricFactor())*(1 -
                std::sqrt(component.reducedTemperature(temperature))),2)*(
                component.reducedPressure(pressure)/
                std::pow(component.reducedTemperature(temperature),2));
        }
    };
}

#endif /* PENG_ROBINSON_EOS */