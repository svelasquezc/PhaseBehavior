#ifndef GENERALIZED_CUBIC_EOS_HPP
#define GENERALIZED_CUBIC_EOS_HPP

#include <vector>

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
            public:
            
            std::vector<NP_t> operator() (Mixture const& mixture, NP_t const& pressure, NP_t const& temperature){
                
            }

        };
    }
}

#endif /* GENERALIZED_CUBIC_EOS_HPP */