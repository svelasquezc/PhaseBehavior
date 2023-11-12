#ifndef MIXTURE_HPP
#define MIXTURE_HPP

#include <vector>
#include <memory>
#include "Component.hpp"

class Mixture {
private: 

    using NP_t = TypesDefinition::NumericalPrecision;

    class MixtureComponent {
        std::shared_ptr<Component> pureComponent;
        NP_t molarComposition_;
    };

    mutable unsigned int numberOfComponents_;

    NP_t attractionParameter_;
    NP_t covolumeParameter_;

    std::vector<MixtureComponent> components;

public:
    Mixture(){}
    
};

#endif /* MIXTURE_HPP */