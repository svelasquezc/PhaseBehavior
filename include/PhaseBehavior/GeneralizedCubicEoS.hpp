#ifndef GENERALIZED_CUBIC_EOS_HPP
#define GENERALIZED_CUBIC_EOS_HPP

#include <vector>
#include <cmath>
#include <limits>

#include "Utilities/Math.hpp"
#include "Utilities/Types.hpp"
#include "Mixture.hpp"


namespace PhaseBehavior {

    using NP_t = Types::NumericalPrecision;

    namespace EoS {

        template<
            Types::ConstexprConstant EoSParameter1,
            Types::ConstexprConstant EoSParameter2,
            typename MixingRule
            >
        class GeneralizedCubicEoS {
            private:
            constexpr static NP_t m1 = EoSParameter1();
            constexpr static NP_t m2 = EoSParameter2();
            MixingRule mixingRule;
            public:
            
            std::vector<NP_t> operator() (Mixture const& mixture, NP_t const& pressure, NP_t const& temperature){
                auto [mixtureAttraction, mixtureCovolume] = mixingRule(mixture, pressure, temperature);

                NP_t cuadraticTerm = (m1 + m2 - 1)*mixtureCovolume - 1;
                NP_t linearTerm = mixtureAttraction + m1*m2*std::pow(mixtureCovolume,2) - (m1 + m2)*mixtureCovolume*(mixtureCovolume + 1);
                NP_t constantTerm = -1*(mixtureAttraction*mixtureCovolume + m1*m2*std::pow(mixtureCovolume,2)*(mixtureCovolume + 1));

                auto roots = Math::cubicRootsFind(cuadraticTerm, linearTerm, constantTerm);

                if (roots.size()==1){
                    return {roots[0]};
                }

                NP_t zMin = std::numeric_limits<NP_t>::max();
                NP_t zMax = std::numeric_limits<NP_t>::min();
                
                for (const auto& root : roots){
                    if (root>zMax) zMax = root;
                    if (root<zMin) zMin = root; 
                }

                return {zMax, zMin};
            }

        };
    }
}

#endif /* GENERALIZED_CUBIC_EOS_HPP */