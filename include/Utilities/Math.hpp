#ifndef MATH_UTILITIES_HPP
#define MATH_UTILITIES_HPP
#include <vector>
#include <cmath>

namespace Math {

    namespace constants{
        constexpr double pi() { return std::atan(1)*4; };
    };

    template<typename PT>
    constexpr std::vector<PT> cubicRootsFind (PT a, PT b, PT c){
        
        auto Q = (std::pow(a,2) - 3*b)/9;
        auto R = (2*std::pow(a,3) - 9*a*b + 27*c)/54;

        auto Q3 = std::pow(Q,3);
        auto R2 = std::pow(R,2);

        auto M = R2 - Q3;

        if (R2 > Q3){
            auto S = -std::copysign(1.0, R)*std::cbrt(std::abs(R) + std::sqrt(M));
            auto T = S > 0.0 ? Q/S : 0;
            return {S+T-a/3};
        }else {
            auto theta = std::acos(R/std::sqrt(Q3));
            return {
            -(2*std::sqrt(Q)*std::cos(theta/3) - a/3), 
            -(2*std::sqrt(Q)*std::cos((theta+2*constants::pi())/3) - a/3),
            -(2*std::sqrt(Q)*std::cos((theta-2*constants::pi())/3) - a/3)
            };
        }
    }
};
#endif /* MATH_UTILITIES_HPP */