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
#include <optional>
#include "Component.hpp"

namespace PhaseBehavior {

    class Mixture {
    private: 

        using NP_t = Types::NumericalPrecision;

        class MixtureComponent {
        private:
            std::shared_ptr<Component> mutable pureComponent_;
            std::map<std::string, NP_t> molarComposition_;
            std::map<std::string, NP_t> fugacityCoefficient_;
            NP_t equilibriumCoefficient_;
        public:

            MixtureComponent(MixtureComponent& rhs) = default;
            MixtureComponent(MixtureComponent&& rhs) = default;

            MixtureComponent& operator= (MixtureComponent const&) = default;
            MixtureComponent& operator= (MixtureComponent&&) = default;

            template <typename PairType,
            typename = typename std::enable_if_t<!std::is_same_v<PairType, MixtureComponent>>>
            explicit MixtureComponent(PairType&& componentWithComposition)
            {
                pureComponent_ = std::make_shared<Component>(std::forward<Component>(std::get<0>(componentWithComposition)));
                molarComposition_["global"]  = std::forward<NP_t>(std::get<1>(componentWithComposition));
            }

            MixtureComponent(Component&& pureComponent, NP_t&& molarComposition):
            pureComponent_(std::make_shared<Component>(std::forward<Component>(pureComponent)))
            {
                molarComposition_["global"] = std::move(molarComposition);
            }

            void composition(NP_t const& molarComposition){
                molarComposition_["global"] = molarComposition;
            }

            const NP_t composition() const noexcept{
                return molarComposition_.at("global");
            }

            const NP_t composition(std::string const& phaseName) const noexcept {
                return molarComposition_.at(phaseName);
            }

            void composition (std::string const& phaseName,  NP_t const& composition){
                molarComposition_[phaseName] = composition;
            }

            const NP_t fugacityCoefficient(std::string const& phaseName) const noexcept{
                return fugacityCoefficient_.at(phaseName);
            }

            void fugacityCoefficient (std::string const& phaseName,  NP_t const& fugacityCoefficient){
                fugacityCoefficient_[phaseName] = fugacityCoefficient;
            }

            const NP_t fugacity(std::string const& phaseName, NP_t const& absolutePressure) const noexcept{
                return composition(phaseName)*fugacityCoefficient(phaseName)*absolutePressure;
            }

            const decltype(*pureComponent_) pure () const{
                return *pureComponent_;
            }

            NP_t equilibriumCoefficient() const {
                return equilibriumCoefficient_;
            }

            void equilibriumCoefficient(NP_t const& equilibriumCoefficient){
                equilibriumCoefficient_ = equilibriumCoefficient;
            }

            bool operator< (MixtureComponent const& rhs) const {
                return pure().name().compare(rhs.pure().name()) < 0;
            }

            bool operator== (MixtureComponent const& rhs) const {
                return pure().name().compare(rhs.pure().name()) == 0;
            }

        };

        unsigned int mutable numberOfComponents_;

        std::vector<MixtureComponent> components_;

        using InteractionCoefficientsType = std::map<std::pair<MixtureComponent const&,MixtureComponent const&>, NP_t>;

        InteractionCoefficientsType interactionCoefficients_;

        using MixtureIterator = decltype(components_)::iterator;
        using ConstMixtureIterator = decltype(components_)::const_iterator;

        std::map<std::string, NP_t> phaseMolarFraction_;
        std::map<std::string, NP_t> phaseCompressibility_;

    public:

        Mixture() = delete;
        Mixture(Mixture const&) = default;
        Mixture(Mixture&&) = default;

        Mixture& operator= (Mixture const&) = default;
        Mixture& operator= (Mixture&&) = default;

        template< 
            typename ContainerType,
            typename = typename std::enable_if_t<
            std::is_constructible_v<MixtureComponent, typename ContainerType::value_type>
            >
            >
        Mixture(ContainerType& inputComponents): 
        numberOfComponents_(inputComponents.size())
        {
            assert(numberOfComponents_ >= 1 && "There should be one or more components in a mixture");

            for(auto&& pair : inputComponents){
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
            static_assert(sizeof...(args) >= 1, "There should be one or more components in a mixture");  
            
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

        auto size() const {return components_.size();}

        MixtureIterator begin() {return components_.begin();}
        MixtureIterator end() {return components_.end();}

        ConstMixtureIterator begin() const {return components_.begin();}
        ConstMixtureIterator end() const {return components_.end();}

        void interactionCoefficient(std::string const& componentName1, std::string const& componentName2, NP_t value){
            auto mixComp1 = (*this)[componentName1];
            auto mixComp2 = (*this)[componentName2];
            assert(mixComp1!=this->end() && (componentName1 + " is not a component of the mixture").c_str());
            assert(mixComp2!=this->end() && (componentName2 + " is not a component of the mixture").c_str());
            interactionCoefficients_[{*mixComp1, *mixComp2}] =  value;
        }

        void interactionCoefficient(MixtureComponent const& component1, MixtureComponent const& component2, NP_t value){
            interactionCoefficients_[{component1, component2}] = value;
        }

        NP_t interactionCoefficient(MixtureComponent const& component1, MixtureComponent const& component2) const{
            try{
                return interactionCoefficients_.at({component1, component2});
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
                return interactionCoefficients_.at({*mixComp1, *mixComp2});
            }catch(std::out_of_range e){
                return 0.0;
            }
        }

        NP_t molarFraction(std::string const& phaseName) const {
            return phaseMolarFraction_.at(phaseName);
        }

        void molarFraction(std::string const& phaseName, NP_t const& molarFraction) {
            phaseMolarFraction_[phaseName] = molarFraction;
        }

        NP_t compressibility(std::string const& phaseName) const {
            return phaseCompressibility_.at(phaseName);
        }

        void compressibility(std::string const& phaseName, NP_t const& compressibility) {
            phaseCompressibility_[phaseName] = compressibility;
        }

        void initializeEquilibriumCoefficients(NP_t const& pressure, NP_t const& temperature){
            for (auto& mixComp : components_){
                mixComp.equilibriumCoefficient(mixComp.pure().equilibriumCoefficient(pressure, temperature));
            }
        }

        NP_t pseudoCriticalPressure(std::string composititionType = "global") const {
            return std::accumulate(components_.begin(), components_.end(), static_cast<NP_t>(0), [&composititionType](auto previous, auto& element){
                return previous + element.composition(composititionType)*element.pure().criticalPressure();
            });
        }

        NP_t pseudoCriticalTemperature(std::string composititionType = "global") const {
            return std::accumulate(components_.begin(), components_.end(), static_cast<NP_t>(0), [&composititionType](auto previous, auto& element){
                return previous + element.composition(composititionType)*element.pure().criticalTemperature();
            });
        }

        NP_t pseudoCriticalVolume(std::string composititionType = "global") const {
            return std::accumulate(components_.begin(), components_.end(), static_cast<NP_t>(0), [&composititionType](auto previous, auto& element){
                return previous + element.composition(composititionType)*element.pure().criticalVolume()*element.pure().molarWeight();
            });
        }
    };
}

#endif /* MIXTURE_HPP */