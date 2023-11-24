#ifndef MIXING_RULES_HPP
#define MIXING_RULES_HPP

#include <tuple>
#include <map>
#include <cmath>
#include <functional>
#include "Utilities/Types.hpp"
#include "Component.hpp"
#include "Mixture.hpp"


namespace PhaseBehavior::EoS::MixingRules{

    using NP_t = Types::NumericalPrecision;

    template<
    Types::ConstexprConstant attractionParameter,
    Types::ConstexprConstant covolumeParameter,
    typename EoS
    >
    class NonRandomMixingRule{
    private:
        using MapType = std::map<std::reference_wrapper<Component const>, NP_t, std::less<const Component>>;

        MapType componentAParameter, componentBParameter;
    public:
        std::tuple<NP_t, NP_t> operator() (Mixture const& mixture, NP_t const& pressure, NP_t const& temperature){

            for(auto const& mixComponent : mixture){

                componentAParameter[mixComponent.pure()] = attractionParameter()*
                EoS::componentAttraction(mixComponent.pure(), pressure, temperature);

                componentBParameter[mixComponent.pure()] = covolumeParameter()*
                mixComponent.pure().reducedPressure(pressure)/mixComponent.pure().reducedTemperature(temperature);

            }

            NP_t mixtureAttraction = 0;
            NP_t mixtureCovolume = 0;

            for(auto const& mixComponentI : mixture){
                for(auto const& mixComponentJ: mixture){

                    auto binaryCoefficient = (1 - mixture.interactionCoefficient(mixComponentI, mixComponentJ))*
                    std::sqrt(componentAParameter[mixComponentI.pure()]*componentAParameter[mixComponentJ.pure()]);

                    mixtureAttraction += mixComponentI.composition()*mixComponentJ.composition()*binaryCoefficient;

                }

                mixtureCovolume += mixComponentI.composition()*componentBParameter[mixComponentI.pure()];
            }

            return {mixtureAttraction, mixtureCovolume};
        }
    
    };
}

#endif /* MIXING_RULES_HPP */