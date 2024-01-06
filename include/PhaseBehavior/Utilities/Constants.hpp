#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include "Types.hpp"

using NP_t = PhaseBehavior::Types::NumericalPrecision;

namespace PhaseBehavior::Constants {
    constexpr NP_t universalGasesConstant = 10.73; //psia ft3/lb · mol · R
    
}

#endif /* CONSTANTS_HPP */