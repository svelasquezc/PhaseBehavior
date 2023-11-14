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
    private:
        std::shared_ptr<Component> const pureComponent_;
        NP_t const molarComposition_;

        const NP_t composition() const {
            return molarComposition_;
        }

        friend class Mixture;
    public:
        MixtureComponent(Component&& pureComponent, NP_t&& molarComposition):
        pureComponent_(std::make_shared<Component>(std::forward<Component>(pureComponent))),
        molarComposition_(std::move(molarComposition)){}

        MixtureComponent(MixtureComponent& rhs) = default;
        MixtureComponent(MixtureComponent&& rhs) = default;


    };

    unsigned int mutable numberOfComponents_;

    NP_t attractionParameter_;
    NP_t covolumeParameter_;

    std::vector<MixtureComponent> components_;

public:
    template<
        typename... MixtureComponents,
        typename = std::enable_if_t<
            std::conjunction_v<std::is_same<MixtureComponent, MixtureComponents>...>
            >
        >
    Mixture(MixtureComponents&&... args): 
    numberOfComponents_(sizeof...(args)),
    components_({std::forward<Args>(args)...}){}
    
    template<typename Pair,
    typename = std::enable_if_t<std::is_constructible<MixtureComponent, Pair>::type>
    >
    Mixture(std::vector<Pair>& inputComponents): 
    numberOfComponents_(inputComponents.size()),
    components_(inputComponents){}
};

#endif /* MIXTURE_HPP */