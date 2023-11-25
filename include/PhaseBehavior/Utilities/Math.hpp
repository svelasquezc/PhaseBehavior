#ifndef MATH_UTILITIES_HPP
#define MATH_UTILITIES_HPP

#include <vector>
#include <cmath>

namespace PhaseBehavior::Math {

    // Namespace containing mathematical constants
    namespace constants{
        /**
         * @brief Constant for pi.
         *
         * This constant represents the mathematical constant pi.
         */
        template<typename PT>
        constexpr PT pi() { return std::atan(static_cast<PT>(1)) * static_cast<PT>(4); };
    };

    /**
     * @brief Function template to find cubic roots of a cubic equation in standard form.
     *
     * This template function calculates the cubic roots of a cubic equation
     * in the standard form x^3 + ax^2 + bx + c, given coefficients a, b, and c.
     *
     * @tparam PT Type of the coefficients (e.g., long double, double, float).
     * @param a Coefficient of the quadratic term.
     * @param b Coefficient of the linear term.
     * @param c Constant term.
     * @return Vector containing the cubic roots.
     */
    template<typename PT>
    constexpr std::vector<PT> cubicRootsFind(PT a, PT b, PT c) {
        // Calculate intermediate variables for cubic root calculation
        PT Q = (std::pow(a, 2) - 3 * b) / 9.0;
        PT R = (2 * std::pow(a, 3) - 9 * a * b + 27 * c) / 54.0;

        PT Q3 = std::pow(Q, 3.0);
        PT R2 = std::pow(R, 2.0);

        PT M = R2 - Q3;

        // Check conditions for different cases of cubic roots
        if (R2 > Q3) {
            // One root
            PT S = -std::copysign(1.0, R) * std::cbrt(std::abs(R) + std::sqrt(M));
            PT T = S != 0.0 ? Q / S : 0;
            return {S + T - a / 3};
        } else {
            // Three roots
            PT theta = std::acos(R / std::sqrt(Q3));
            return {
                -(2 * std::sqrt(Q) * std::cos(theta / 3)) - a / 3,
                -(2 * std::sqrt(Q) * std::cos((theta - 2 * constants::pi<PT>()) / 3)) - a / 3,
                -(2 * std::sqrt(Q) * std::cos((theta + 2 * constants::pi<PT>()) / 3)) - a / 3
                
            };
        }
    }
};

#endif /* MATH_UTILITIES_HPP */
