#ifndef MIXTURE_HPP
#define MIXTURE_HPP

#include <vector>
#include <memory>
#include <type_traits>
#include "Component.hpp"

class Mixture {
private: 

    using NP_t = TypesDefinition::NumericalPrecision;

    class MixtureComponent {
        std::shared_ptr<Component> const pureComponent;
        NP_t const molarComposition_;
    };

    unsigned int mutable numberOfComponents_;

    NP_t attractionParameter_;
    NP_t covolumeParameter_;

    std::vector<MixtureComponent> components;

public:
    template<
        typename... MixtureComponents,
        typename = std::enable_if_t<
            std::conjunction_v<std::is_same<MixtureComponent, MixtureComponents>...>
            >
        >
    Mixture(MixtureComponents&&... args): 
    numberOfComponents_(sizeof...(args)),
    components(std::forward<Args>(args)...){}
    
};

#endif /* MIXTURE_HPP */