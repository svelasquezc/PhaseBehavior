#ifndef MIXTURE_HPP
#define MIXTURE_HPP

#include <array>
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

    enum class PhaseName : std::size_t { global = 0, vapor = 1, liquid = 2 };

    static PhaseName phaseNameFromString(std::string const& phaseName);

    static constexpr std::size_t phaseIndex(PhaseName phaseName) noexcept {
        return static_cast<std::size_t>(phaseName);
    }

    static std::size_t phaseIndex(std::string const& phaseName) noexcept {
        return phaseIndex(phaseNameFromString(phaseName));
    }

    static PhaseName phaseNameFromString(std::string const& phaseName){
        if (phaseName == "vapor") return PhaseName::vapor;
        if (phaseName == "liquid") return PhaseName::liquid;
        return PhaseName::global;
    }

    class Mixture {
    private: 

        using NP_t = Types::NumericalPrecision;

        class MixtureComponent {
        private:
            std::shared_ptr<Component> mutable pureComponent_;
            std::array<NP_t, 3> molarComposition_{};
            std::array<NP_t, 3> fugacityCoefficient_{};
            NP_t equilibriumCoefficient_;
        public:

            MixtureComponent(MixtureComponent const& rhs) = default;
            MixtureComponent(MixtureComponent&& rhs) = default;

            MixtureComponent& operator= (MixtureComponent const&) = default;
            MixtureComponent& operator= (MixtureComponent&&) = default;

            template <typename PairType,
            typename = typename std::enable_if_t<
            !std::is_same_v<std::decay_t<PairType>, MixtureComponent>
            &&
            !std::is_same_v<std::decay_t<PairType>, Mixture>
            >>
            explicit MixtureComponent(PairType&& componentWithComposition)
            {

                pureComponent_ = std::make_shared<Component>(std::forward<Component>(std::get<0>(componentWithComposition)));
                molarComposition_[phaseIndex(PhaseName::global)] = std::forward<NP_t>(std::get<1>(componentWithComposition));
            }

            MixtureComponent(Component&& pureComponent, NP_t&& molarComposition):
            pureComponent_(std::make_shared<Component>(std::forward<Component>(pureComponent)))
            {
                molarComposition_[phaseIndex(PhaseName::global)] = std::move(molarComposition);
            }

            void composition(NP_t const& molarComposition){
                molarComposition_[phaseIndex(PhaseName::global)] = molarComposition;
            }

            const NP_t composition() const noexcept{
                return molarComposition_[phaseIndex(PhaseName::global)];
            }

            const NP_t composition(PhaseName phaseName) const noexcept {
                return molarComposition_[phaseIndex(phaseName)];
            }

            const NP_t composition(std::string const& phaseName) const noexcept {
                return composition(phaseNameFromString(phaseName));
            }

            void composition(PhaseName phaseName, NP_t const& composition){
                molarComposition_[phaseIndex(phaseName)] = composition;
            }

            void composition(std::string const& phaseName, NP_t const& compositionValue){
                composition(phaseNameFromString(phaseName), compositionValue);
            }

            const NP_t fugacityCoefficient(PhaseName phaseName) const noexcept{
                return fugacityCoefficient_[phaseIndex(phaseName)];
            }

            const NP_t fugacityCoefficient(std::string const& phaseName) const noexcept{
                return fugacityCoefficient(phaseNameFromString(phaseName));
            }

            void fugacityCoefficient(PhaseName phaseName, NP_t const& fugacityCoefficient){
                fugacityCoefficient_[phaseIndex(phaseName)] = fugacityCoefficient;
            }

            void fugacityCoefficient(std::string const& phaseName, NP_t const& fugacityCoefficientValue){
                fugacityCoefficient(phaseNameFromString(phaseName), fugacityCoefficientValue);
            }

            const NP_t fugacity(PhaseName phaseName, NP_t const& absolutePressure) const noexcept{
                return composition(phaseName)*fugacityCoefficient(phaseName)*absolutePressure;
            }

            const NP_t fugacity(std::string const& phaseName, NP_t const& absolutePressure) const noexcept{
                return fugacity(phaseNameFromString(phaseName), absolutePressure);
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

            auto operator->() const {
                return pureComponent_;
            }

        };

        unsigned int mutable numberOfComponents_;

        std::vector<MixtureComponent> components_;

        using InteractionCoefficientsType = std::map<std::pair<std::string,std::string>, NP_t>;

        InteractionCoefficientsType interactionCoefficients_;

        using MixtureIterator = decltype(components_)::iterator;
        using ConstMixtureIterator = decltype(components_)::const_iterator;

        std::array<NP_t, 3> phaseMolarFraction_{};
        std::array<NP_t, 3> phaseCompressibility_{};

    public:

        Mixture() = delete;
        Mixture(Mixture const&) = default;

        Mixture(Mixture&&) = default;

        Mixture& operator= (Mixture const&) = default;
        Mixture& operator= (Mixture&&) = default;

        void checkConsistency() const {
            assert(std::accumulate(components_.begin(),components_.end(), 0,
            [](auto previous, auto second){
                return previous + second.composition();
            }) - 1 <= std::numeric_limits<NP_t>::epsilon() && "The sum of the compositions must be equal to 1");
        }

        template< 
            typename ContainerType,
            typename = typename std::enable_if_t<
            std::is_constructible_v<MixtureComponent, typename ContainerType::value_type> 
            && 
            !std::is_same_v<std::decay_t<ContainerType>, Mixture>
            >
            >
        Mixture(ContainerType& inputComponents): 
        numberOfComponents_(inputComponents.size())
        {
            assert(numberOfComponents_ >= 1 && "There should be one or more components in a mixture");

            for(auto&& pair : inputComponents){
                components_.emplace_back(std::forward<Component>(std::get<Component>(pair)), std::forward<NP_t>(std::get<NP_t>(pair)));
            }

            checkConsistency();

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

            checkConsistency();
        }

        Mixture(std::initializer_list<std::pair<Component, NP_t>>&& compList){

            assert(compList.size() >= 1 && "There should be one or more components in a mixture");

            for (auto mixComp : compList){
                components_.emplace_back(mixComp);
            }

            checkConsistency();
        }

        MixtureIterator operator[] (std::string const& name){
            return std::find_if(begin(), end(), [&name](auto& component){
                return component.pure().name() == name;
            });
        } 

        MixtureComponent& operator[] (std::size_t componentIndex){
            return components_[componentIndex];
        }

        const MixtureComponent& operator[] (std::size_t componentIndex) const{
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
            interactionCoefficients_[{componentName1, componentName2}] =  value;
        }

        void interactionCoefficient(MixtureComponent const& component1, MixtureComponent const& component2, NP_t value){
            interactionCoefficients_[{component1.pure().name(), component2.pure().name()}] = value;
        }

        NP_t interactionCoefficient(MixtureComponent const& component1, MixtureComponent const& component2) const{
            try{
                return interactionCoefficients_.at({component1.pure().name(), component2.pure().name()});
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
                return interactionCoefficients_.at({componentName1, componentName2});
            }catch(std::out_of_range e){
                return 0.0;
            }
        }

        NP_t molarFraction(PhaseName phaseName) const {
            return phaseMolarFraction_[phaseIndex(phaseName)];
        }

        NP_t molarFraction(std::string const& phaseName) const {
            return molarFraction(phaseNameFromString(phaseName));
        }

        void molarFraction(PhaseName phaseName, NP_t const& molarFractionValue) {
            phaseMolarFraction_[phaseIndex(phaseName)] = molarFractionValue;
        }

        void molarFraction(std::string const& phaseName, NP_t const& molarFractionValue) {
            molarFraction(phaseNameFromString(phaseName), molarFractionValue);
        }

        NP_t compressibility(PhaseName phaseName) const {
            return phaseCompressibility_[phaseIndex(phaseName)];
        }

        NP_t compressibility(std::string const& phaseName) const {
            return compressibility(phaseNameFromString(phaseName));
        }

        void compressibility(PhaseName phaseName, NP_t const& compressibilityValue) {
            phaseCompressibility_[phaseIndex(phaseName)] = compressibilityValue;
        }

        void compressibility(std::string const& phaseName, NP_t const& compressibilityValue) {
            compressibility(phaseNameFromString(phaseName), compressibilityValue);
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

        NP_t pseudoCriticalDensity(std::string composititionType = "global") const {
            return std::accumulate(components_.begin(), components_.end(), static_cast<NP_t>(0), [&composititionType](auto previous, auto& element){
                return previous + element.composition(composititionType)*element.pure().criticalDensity()*element.pure().molarWeight();
            });
        }
    };
}

#endif /* MIXTURE_HPP */