#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include "Types.hpp"

using NP_t = PhaseBehavior::Types::NumericalPrecision;

namespace PhaseBehavior::Constants {
    namespace SI {
        constexpr NP_t universalGasesConstant = 8.314472; // J \cdot mol^{-1} \cdot K^{-1} 
        constexpr NP_t standardConditionsPressure = 101.3529; // kPa
        constexpr NP_t standardConditionsTemperature = 298; //K
    }

    namespace Field{
        constexpr NP_t universalGasesConstant = 10.73; //psia ft3/lb · mol · R    
        constexpr NP_t standardConditionsPressure = 14.7; //psia
        constexpr NP_t standardConditionsTemperature = 520; //R
    }

}

#endif /* CONSTANTS_HPP */