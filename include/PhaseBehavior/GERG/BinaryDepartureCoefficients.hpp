#ifndef BINARYDEPARTURECOEFFICIENTS_HPP
#define BINARYDEPARTURECOEFFICIENTS_HPP

#include <tuple>

#include "../Utilities/Types.hpp"

#include "Types.hpp"

using NP_t = PhaseBehavior::Types::NumericalPrecision;

namespace PhaseBehavior::EoS::GERG::Coefficients::Departure {
    using PolynomialExponentType = GERG::Types::BinaryParameter<std::vector<std::pair<int, NP_t>>>;
}

#endif /* BINARYDEPARTURECOEFFICIENTS_HPP */