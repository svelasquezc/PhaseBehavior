#ifndef GENERALIZED_CUBIC_EOS_HPP
#define GENERALIZED_CUBIC_EOS_HPP

#include <vector>

#include "Utilities/Math.hpp"
#include "Utilities/Types.hpp"
#include "Mixture.hpp"


namespace PhaseBehavior {

    using NP_t = TypesDefinition::NumericalPrecision;
    using ConstexprConstantType = decltype(Math::constants::pi<Precision_t>);

    namespace EoS {

        template<
            ConstexprConstantType EoSParameter1,
            ConstexprConstantType EoSParameter2,
            ConstexprConstantType attractionParameter,
            ConstexprConstantType covolumeParameter,
            typename DerivedEoS
            >
        class GeneralizedCubicEoS {
            private:
            constexpr static NP_t m1 = EoSParameter1();
            constexpr static NP_t m2 = EoSParameter2();
            constexpr static NP_t omegaA = attractionParameter();
            constexpr static NP_t omegaB = covolumeParameter();
            DerivedEoS& self() { return static_cast<DerivedEoS&>(*this); };
            public:
            
            std::vector<NP_t> operator() (Mixture const& mixture, NP_t const& pressure, NP_t const& temperature){
                
            }

        };
    }
}

#endif /* GENERALIZED_CUBIC_EOS_HPP */