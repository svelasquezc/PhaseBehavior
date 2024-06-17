#ifndef PURERESIDUALHELMHOLTZCOEFFICIENTS_HPP
#define PURERESIDUALHELMHOLTZCOEFFICIENTS_HPP

#include <map>
#include <string_view>
#include "../Utilities/Types.hpp"

using NP_t = PhaseBehavior::Types::NumericalPrecision;

namespace PhaseBehavior::EoS::GERG::Coefficients::Residual {
    const static std::map<std::string_view, std::array<NP_t, 7>> residualCoefficient = {
        {"CH4", {19.597508817, -83.959667892, 3.00088, 0.76315, 0.00460, 8.74432, -4.46921}},
        {"N2", {11.083407489, -22.202102428, 2.50031, 0.13732, -0.14660, 0.90066, 0}}
    };
    const static std::map<std::string_view, std::array<NP_t, 4>> exponent = {
        {"CH4", {4.306474465, 0.936220902, 5.577233895, 5.722644361}},
        {"N2", {5.251822620, -5.393067706, 13.788988208, 0}}
    };
}

#endif /* PURERESIDUALHELMHOLTZCOEFFICIENTS_HPP */