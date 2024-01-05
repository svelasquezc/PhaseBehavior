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
    class RandomMixingRule{
    private:
        
        using MapType = std::map<std::reference_wrapper<Component const>, NP_t, std::less<const Component>>;
        using BinaryAttractionType = std::map<std::pair<Component const&, Component const&>, NP_t>;

        MapType componentAttraction_, componentCovolume_;

        BinaryAttractionType binaryAttraction_;

    public:
        std::tuple<NP_t, NP_t> operator() (Mixture const& mixture, NP_t const& pressure, NP_t const& temperature, std::string const& phaseName = "global"){

            for(auto const& mixComponent : mixture){

                componentAttraction_[mixComponent.pure()] = attractionParameter()*
                EoS::componentAttraction(mixComponent.pure(), pressure, temperature);

                componentCovolume_[mixComponent.pure()] = covolumeParameter()*
                mixComponent.pure().reducedPressure(pressure)/mixComponent.pure().reducedTemperature(temperature);

            }

            NP_t mixtureAttraction = 0;
            NP_t mixtureCovolume = 0;

            for(auto const& mixComponentI : mixture){
                for(auto const& mixComponentJ: mixture){

                    auto binaryAttractionValue = (1 - mixture.interactionCoefficient(mixComponentI, mixComponentJ))*
                    std::sqrt(componentAttraction_[mixComponentI.pure()]*componentAttraction_[mixComponentJ.pure()]);

                    binaryAttraction_[{mixComponentI.pure(), mixComponentJ.pure()}] = binaryAttractionValue;

                    mixtureAttraction += mixComponentI.composition(phaseName)*mixComponentJ.composition(phaseName)*binaryAttractionValue;

                }

                mixtureCovolume += mixComponentI.composition(phaseName)*componentCovolume_[mixComponentI.pure()];
            }

            return {mixtureAttraction, mixtureCovolume};
        }

        NP_t componentAttraction(Component const& component) const{
            return componentAttraction_.at(component);
        }

        NP_t componentCovolume(Component const& component) const{
            return componentCovolume_.at(component);
        }

        NP_t binaryAttraction(Component const& componentI, Component const& componentJ) const{
            return binaryAttraction_.at({componentI, componentJ});
        }
    };
}

#endif /* MIXING_RULES_HPP */