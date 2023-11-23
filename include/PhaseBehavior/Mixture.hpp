#ifndef MIXTURE_HPP
#define MIXTURE_HPP

#include <vector>
#include <memory>
#include <type_traits>
#include <utility>
#include <numeric>
#include <cmath>
#include <cassert>
#include <limits>
#include <map>
#include <algorithm>
#include <exception>
#include "Component.hpp"

namespace PhaseBehavior {

    class Mixture {
    private: 

        using NP_t = Types::NumericalPrecision;

        class MixtureComponent {
        private:
            std::shared_ptr<Component> mutable pureComponent_;
            NP_t mutable molarComposition_;

        public:

            MixtureComponent(MixtureComponent& rhs) = default;
            MixtureComponent(MixtureComponent&& rhs) = default;

            template <typename PairType,
            typename = typename std::enable_if_t<!std::is_same_v<PairType, MixtureComponent>>>
            explicit MixtureComponent(PairType&& componentWithComposition)
            {
                pureComponent_ = std::make_shared<Component>(std::forward<Component>(std::get<0>(componentWithComposition)));
                molarComposition_  = std::forward<NP_t>(std::get<1>(componentWithComposition));
            }

            MixtureComponent(Component&& pureComponent, NP_t&& molarComposition):
            pureComponent_(std::make_shared<Component>(std::forward<Component>(pureComponent))),
            molarComposition_(std::move(molarComposition)){}

            void composition(NP_t const& molarComposition){
                molarComposition_ = molarComposition;
            }

            const NP_t composition() const {
                return molarComposition_;
            }

            const decltype(*pureComponent_) pure () const{
                return *pureComponent_;
            }

            bool operator< (MixtureComponent const& rhs) const {
                return pure().name().compare(rhs.pure().name()) < 0;
            }
        };

        unsigned int mutable numberOfComponents_;

        NP_t attractionParameter_;
        NP_t covolumeParameter_;

        std::vector<MixtureComponent> components_;

        using InteractionCoefficientsType = std::map<std::pair<MixtureComponent const&,MixtureComponent const&>, NP_t>;

        InteractionCoefficientsType interactionCoefficients_;

        using MixtureIterator = decltype(components_)::iterator;
        using ConstMixtureIterator = decltype(components_)::const_iterator;

    public:

        Mixture() = delete;
        Mixture(Mixture&) = delete;
        Mixture(Mixture&&) = delete;

        template< 
            typename ContainerType,
            typename = typename std::enable_if_t<
            std::is_constructible_v<MixtureComponent, typename ContainerType::value_type>
            >
            >
        Mixture(ContainerType& inputComponents): 
        numberOfComponents_(inputComponents.size())
        {
            assert(numberOfComponents_ > 1 && "There should be more than one component in a mixture");

            for(auto&& pair : inputComponents){
                //TypeChecker<decltype(pair)> tc;
                components_.emplace_back(std::forward<Component>(pair.first), std::forward<NP_t>(pair.second));
            }

            assert(std::accumulate(components_.begin(),components_.end(), 0,
            [](auto previous, auto second){
                return previous + second.composition();
            }) - 1 <= std::numeric_limits<NP_t>::epsilon() && "The sum of the compositions must be equal to 1");

        }

        template<
            typename... MixtureComponents,
            typename = typename std::enable_if_t<
                std::conjunction_v<std::is_constructible<MixtureComponent, MixtureComponents>...>
                >
            >
        Mixture(MixtureComponents&&... args): 
        numberOfComponents_(sizeof...(args))
        {
            static_assert(sizeof...(args) > 1, "There should be more than one component in a mixture");  
            //TypeChecker<decltype(args)...> tc;
            (components_.emplace_back(
                std::forward<std::common_type_t<std::decay_t<MixtureComponents>...>>(args)
                ),...);

            assert(std::accumulate(components_.begin(),components_.end(), 0,
            [](auto previous, auto second){
                return previous + second.composition();
            }) - 1 <= std::numeric_limits<NP_t>::epsilon() && "The sum of the compositions must be equal to 1");
        }

        ConstMixtureIterator operator[] (std::string const& name){
            return std::find_if(begin(), end(), [&name](auto& component){
                return component.pure().name() == name;
            });
        } 

        auto operator[] (std::size_t componentIndex) -> decltype(components_[componentIndex]){
            return components_[componentIndex];
        }

        ConstMixtureIterator begin() const {return components_.begin();}
        ConstMixtureIterator end() const {return components_.end();}

        void interactionCoefficient(std::string const& componentName1, std::string const& componentName2, NP_t value){
            auto mixComp1 = (*this)[componentName1];
            auto mixComp2 = (*this)[componentName2];
            assert(mixComp1!=this->end() && (componentName1 + " is not a component of the mixture").c_str());
            assert(mixComp2!=this->end() && (componentName2 + " is not a component of the mixture").c_str());
            interactionCoefficients_[{*mixComp1, *mixComp2}] = value;
        }

        void interactionCoefficient(MixtureComponent const& component1, MixtureComponent const& component2, NP_t value){
            interactionCoefficients_[{component1, component2}] = value;
        }

        NP_t interactionCoefficient(MixtureComponent const& component1, MixtureComponent const& component2){
            try{
                return interactionCoefficients_[{component1, component2}];
            }catch(std::out_of_range e){
                return 0.0;
            }
        }

        NP_t interactionCoefficient(std::string const& componentName1, std::string const& componentName2){
            auto mixComp1 = (*this)[componentName1];
            auto mixComp2 = (*this)[componentName2];
            assert(mixComp1!=this->end() && (componentName1 + " is not a component of the mixture").c_str());
            assert(mixComp2!=this->end() && (componentName2 + " is not a component of the mixture").c_str());
            try{
                return interactionCoefficients_[{*mixComp1, *mixComp2}];
            }catch(std::out_of_range e){
                return 0.0;
            }
        }
    };
}

#endif /* MIXTURE_HPP */