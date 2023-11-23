#ifndef MIXING_RULES_HPP
#define MIXING_RULES_HPP

#include <tuple>
#include "Utilities/Types.hpp"
#include "Mixture.hpp"


namespace PhaseBehavior::EoS::MixingRules{

    using NP_t = Types::NumericalPrecision;

    template<
    Types::ConstexprConstant attractionParameter,
    Types::ConstexprConstant covolumeParameter,
    typename EoS
    >
    class NonRandomMixingRule{
        
    public:
        std::tuple<NP_t, NP_t> operator() (Mixture const& mixture, NP_t const& pressure, NP_t const& temperature){

        }
    
    };
}

#endif /* MIXING_RULES_HPP */