#ifndef MIXTURE_HPP
#define MIXTURE_HPP

#include <vector>
#include <memory>
#include <type_traits>
#include <utility>
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
        template <typename PairType>
        MixtureComponent(PairType&& componentWithComposition): 
        pureComponent_ (std::make_shared<Component>(std::move(std::get<Component>(componentWithComposition)))),
        molarComposition_ (std::move(std::get<NP_t>(componentWithComposition)))
        {}

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

    Mixture() = delete;
    Mixture(Mixture&) = delete;
    Mixture(Mixture&&) = delete;

    template< 
        typename Container, 
        typename = typename std::enable_if_t<
        std::is_constructible<MixtureComponent, typename Container::value_type>::value
        >
        >
    Mixture(Container& inputComponents): 
    numberOfComponents_(inputComponents.size())
    {
        std::copy(inputComponents.begin(), inputComponents.end(), 
          std::back_inserter(components_));
    }

    template<
        typename... MixtureComponents,
        typename = typename std::enable_if_t<
            std::conjunction_v<std::is_constructible<MixtureComponent, MixtureComponents>...>
            >
        >
    Mixture(MixtureComponents&&... args): 
    numberOfComponents_(sizeof...(args)),
    components_({std::forward<MixtureComponents>(args)...}){}
};

#endif /* MIXTURE_HPP */