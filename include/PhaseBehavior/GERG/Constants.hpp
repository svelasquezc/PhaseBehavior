#ifndef GERG_CONSTANTS_HPP
#define GERG_CONSTANTS_HPP

#include "../Utilities/Types.hpp"

using NP_t = PhaseBehavior::Types::NumericalPrecision;

namespace PhaseBehavior::EoS::GERG::Constants
{
    constexpr NP_t universalGasesConstant = 8.314472; // J \cdot mol^{-1} \cdot K^{-1} 
    constexpr NP_t jaeschkeSchleyMolarGasConstant = 8.314510; // J \cdot mol^{-1} \cdot K^{-1} 
    constexpr NP_t gasConstantRatio = jaeschkeSchleyMolarGasConstant/universalGasesConstant; // dimensionless
} // namespace PhaseBehavior::EoS::GERG::Constants

#endif /* GERG_CONSTANTS_HPP */