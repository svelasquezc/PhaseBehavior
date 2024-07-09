#ifndef GERG2008_HPP
#define GERG2008_HPP
#include "../Utilities/Types.hpp"
#include "../Utilities/Math.hpp"
#include "../Component.hpp"
#include "Constants.hpp"
#include "../Mixture.hpp"

#include "PureFluidHelmholtzCoefficients.hpp"
#include "BinaryReducingParameters.hpp"
#include "BinaryDepartureCoefficients.hpp"

using NP_t = PhaseBehavior::Types::NumericalPrecision;

namespace PhaseBehavior::EoS::GERG{

    using PhaseBehavior::Component;
    using PhaseBehavior::Mixture;

    class GERG2008{

        private:

        enum class WRT {
            None,
            reducedMixtureDensity,
            reducedTemperatureInverse
        };

        NP_t reducingMixtureDensity_, reducingMixtureTemperature_;
        NP_t idealReducedHelmholtzEnergy_ = 0, residualReducedHelmholtzEnergy_ = 0;
        NP_t compressibilityFactor_ = 1;

        inline NP_t reducingDensity(Mixture const& mixture, std::string compositionType = "global"){

            using Coefficients::ReducingParameters::Density::beta, Coefficients::ReducingParameters::Density::gamma;

            reducingMixtureDensity_ = 0;

            for(auto const& component: mixture){
                reducingMixtureDensity_ += std::pow(component.composition(compositionType),2)/component->molarCriticalDensity();
            }

            for(unsigned int i=0; i<mixture.size()-1; ++i){
                for (unsigned int j=i+1; j<mixture.size(); ++j){
                    reducingMixtureDensity_ += 2*mixture[i].composition(compositionType)*mixture[j].composition(compositionType)*
                    beta.at({mixture[i]->name(), mixture[j]->name()})*gamma.at({mixture[i]->name(), mixture[j]->name()})*
                    (mixture[i].composition(compositionType) + mixture[j].composition(compositionType))/
                    (std::pow(beta.at({mixture[i]->name(), mixture[j]->name()}),2)*mixture[i].composition(compositionType) + mixture[j].composition(compositionType))*
                    (1.0/8.0)*std::pow(1/std::cbrt(mixture[i]->molarCriticalDensity()) + 1/std::cbrt(mixture[j]->molarCriticalDensity()),3);
                }
            }

            reducingMixtureDensity_ = 1.0/reducingMixtureDensity_;

            return reducingMixtureDensity_;
        }

        inline NP_t reducingTemperature(Mixture const& mixture, std::string compositionType = "global"){

            using Coefficients::ReducingParameters::Temperature::beta, Coefficients::ReducingParameters::Temperature::gamma;

            reducingMixtureTemperature_ = 0;

            for(auto const& component: mixture){
                reducingMixtureTemperature_ += std::pow(component.composition(compositionType),2)*component->criticalTemperature();
            }

            for(unsigned int i=0; i<mixture.size()-1; ++i){
                for (unsigned int j=i+1; j<mixture.size(); ++j){
                    reducingMixtureTemperature_ += 2*mixture[i].composition(compositionType)*mixture[j].composition(compositionType)*
                    beta.at({mixture[i]->name(), mixture[j]->name()})*gamma.at({mixture[i]->name(), mixture[j]->name()})*
                    (mixture[i].composition(compositionType) + mixture[j].composition(compositionType))/
                    (std::pow(beta.at({mixture[i]->name(), mixture[j]->name()}),2)*mixture[i].composition(compositionType) + mixture[j].composition(compositionType))*
                    std::sqrt(mixture[i]->criticalTemperature()*mixture[j]->criticalTemperature());
                }
            }

            return reducingMixtureTemperature_;
        }

        inline NP_t reducedTemperatureInverse(NP_t const& temperature){
            return reducingMixtureTemperature_/temperature;
        }


        inline NP_t reducedMixtureDensity(NP_t const& density) const{
            return density/reducingMixtureDensity_;
        }

        inline NP_t const idealReducedHelmholtzFreeEnergy(Component const& component, NP_t const& mixtureDensity, NP_t const& temperature) const {

            using Coefficients::PureFluid::Ideal::pureFluid;

            auto const& coefficients = pureFluid.at(component.name());

            return std::log(mixtureDensity/component.molarCriticalDensity()) + Constants::gasConstantRatio * (
                coefficients.polynomialCoefficient[0] + 
                coefficients.polynomialCoefficient[1]*component.criticalTemperature()/temperature + 
                coefficients.polynomialCoefficient[2]*std::log(component.criticalTemperature()/temperature) +
                coefficients.polynomialCoefficient[3]*std::log(std::abs(std::sinh(coefficients.hyperbolicCoefficient[0]*component.criticalTemperature()/temperature))) + 
                coefficients.polynomialCoefficient[5]*std::log(std::abs(std::sinh(coefficients.hyperbolicCoefficient[2]*component.criticalTemperature()/temperature))) - 
                coefficients.polynomialCoefficient[4]*std::log(std::abs(std::sinh(coefficients.hyperbolicCoefficient[1]*component.criticalTemperature()/temperature))) -
                coefficients.polynomialCoefficient[6]*std::log(std::abs(std::sinh(coefficients.hyperbolicCoefficient[3]*component.criticalTemperature()/temperature)))  
                );
        }

        inline NP_t idealReducedHelmholtzFreeEnergy(Mixture const& mixture, NP_t const& mixtureDensity, NP_t const& temperature, std::string compositionType = "global"){
            idealReducedHelmholtzEnergy_ = 0;
            for (auto const& component : mixture){
                idealReducedHelmholtzEnergy_ += component.composition(compositionType)*idealReducedHelmholtzFreeEnergy(component.pure(), mixtureDensity, temperature) +
                              std::log(component.composition(compositionType));
            }
            return idealReducedHelmholtzEnergy_;
        }

        inline NP_t const residualReducedHelmholtzFreeEnergy(Component const& component, NP_t const& reducedMixtureDensity, NP_t const& reducedTemperatureInverse) const {
            using namespace Coefficients::PureFluid::Residual;

            NP_t pureFluidResidualHelmholtzFreeEnergy = 0;

            auto const& coefficients = pureFluid.at(component.name());

            auto const kPol = coefficients.polynomialExponent.size();

            for (unsigned int term = 0; term < kPol; ++term){

                auto const [d_oi, t_oi] = coefficients.polynomialExponent[term];

                pureFluidResidualHelmholtzFreeEnergy += 
                    coefficients.residualCoefficient[term]*std::pow(reducedMixtureDensity, d_oi)*std::pow(reducedTemperatureInverse, t_oi);
            }

            auto const kExp = coefficients.exponentialExponent.size();

            for (unsigned int term = 0; term < kExp; ++term){

                auto const [d_oi, t_oi, c_oi] = coefficients.exponentialExponent[term];

                pureFluidResidualHelmholtzFreeEnergy += 
                    coefficients.residualCoefficient[kPol + term]*
                    std::pow(reducedMixtureDensity, d_oi)*std::pow(reducedTemperatureInverse, t_oi)*
                    std::exp(-std::pow(reducedMixtureDensity, c_oi));
            }

            return pureFluidResidualHelmholtzFreeEnergy;
        };

        inline NP_t residualDeparture(Component const& component1, Component const& component2, NP_t const& reducedMixtureDensity, NP_t const& reducedTemperatureInverse) const{
            using namespace Coefficients::Departure;

            NP_t residualDeparture = 0;

            auto const F = fParameter(component1.name(),component2.name());

            if (F != 0){

                auto const& coefficients = binaryDeparture.at({component1.name(), component2.name()});

                auto const kPol = coefficients.polynomialExponent.size();

                for (unsigned int term = 0; term < kPol; ++term){

                    auto const [d_ij, t_ij] = coefficients.polynomialExponent[term];

                    residualDeparture += 
                        coefficients.coefficient[term]*std::pow(reducedMixtureDensity, d_ij)*std::pow(reducedTemperatureInverse, t_ij);
                }

                auto const kExp = coefficients.exponentialExponent.size();

                for (unsigned int term = 0; term < kExp; ++term){

                    auto const [d_ij, t_ij, eta_ij, varepsilon_ij, beta_ij, gamma_ij] = coefficients.exponentialExponent[term];

                    residualDeparture += 
                        coefficients.coefficient[kPol + term]*std::pow(reducedMixtureDensity, d_ij)*
                        std::pow(reducedTemperatureInverse, t_ij)*std::exp(
                            -eta_ij*std::pow(reducedMixtureDensity - varepsilon_ij, 2) - beta_ij*(reducedMixtureDensity - gamma_ij)
                        );
                }
            }

            return F*residualDeparture;
        }

        inline NP_t residualReducedHelmholtzFreeEnergy(Mixture const& mixture, NP_t const& reducedMixtureDensity, NP_t const& reducedTemperatureInverse, std::string compositionType = "global"){
            residualReducedHelmholtzEnergy_ = 0;
            for (auto const& component : mixture){
                residualReducedHelmholtzEnergy_ += component.composition(compositionType)*
                                                    residualReducedHelmholtzFreeEnergy(component.pure(), reducedMixtureDensity, reducedTemperatureInverse);
            }

            for(unsigned int i=0; i<mixture.size()-1; ++i){
                for (unsigned int j=i+1; j<mixture.size(); ++j){
                     residualReducedHelmholtzEnergy_ += mixture[i].composition(compositionType)*mixture[j].composition(compositionType)*
                                                        residualDeparture(mixture[i].pure(), mixture[j].pure(), reducedMixtureDensity, reducedTemperatureInverse);
                }
            }

            return residualReducedHelmholtzEnergy_;
        }

        // WRT::None means first order derivative
        // This templated function accounts for derivatives with respect to delta and tau 
        template <WRT wrt1, WRT wrt2 = WRT::None>
        inline NP_t idealReducedHelmholtzFreeEnergyDerivative(Mixture const& mixture, NP_t const& mixtureDensity, NP_t const& temperature, std::string compositionType = "global"){
            // First-Order Derivatives
            if constexpr (wrt2 == WRT::None){
                if constexpr (wrt1 == WRT::reducedMixtureDensity){
                    return reducingMixtureDensity_/mixtureDensity;
                }

                if constexpr (wrt1 == WRT::reducedTemperatureInverse){

                    using Coefficients::PureFluid::Ideal::pureFluid; 

                    NP_t idealReducedHelmholtzFreeEnergyDerivativeValue = 0;

                    for (auto const& component : mixture){
                    
                        auto const& coefficients = pureFluid.at(component->name());

                        idealReducedHelmholtzFreeEnergyDerivativeValue += 
                            component.composition(compositionType)*component->criticalTemperature()/reducingMixtureTemperature_*
                            idealReducedHelmholtzFreeEnergyDerivative<WRT::reducedTemperatureInverse>(component.pure(), mixtureDensity, temperature);
                    }
                }
            }

            // Second-Order Derivatives $\partial \delta^2$, $\partial \delta \partial \tau$
            if constexpr (wrt2 == WRT::reducedMixtureDensity){
                if constexpr (wrt1 == wrt2) {
                    return -std::pow(reducingMixtureDensity_/mixtureDensity,2);
                }

                if constexpr (wrt1 == WRT::reducedTemperatureInverse) return 0;
            }

            if constexpr (wrt2 == WRT::reducedTemperatureInverse && wrt1 == wrt2){
                using Coefficients::PureFluid::Ideal::pureFluid; 

                NP_t idealReducedHelmholtzFreeEnergyDerivativeValue = 0;

                for (auto const& component : mixture){
                
                    auto const& coefficients = pureFluid.at(component->name());

                    idealReducedHelmholtzFreeEnergyDerivativeValue += 
                        component.composition(compositionType)*std::pow(component->criticalTemperature()/reducingMixtureTemperature_,2)*
                        idealReducedHelmholtzFreeEnergyDerivative<WRT::reducedTemperatureInverse, WRT::reducedTemperatureInverse>(component.pure(), mixtureDensity, temperature);
                }
            }
        }

        // WRT::None means first order derivative
        // This templated function accounts for derivatives with respect to delta and tau 
        template <WRT wrt1, WRT wrt2 = WRT::None>
        inline NP_t idealReducedHelmholtzFreeEnergyDerivative(Component const& component, NP_t const& mixtureDensity, NP_t const& temperature){
            // First-Order Derivatives
            if constexpr (wrt2 == WRT::None){
                if constexpr (wrt1 == WRT::reducedMixtureDensity){
                    return component.molarCriticalDensity()/mixtureDensity;
                }

                if constexpr (wrt1 == WRT::reducedTemperatureInverse){

                    using Coefficients::PureFluid::Ideal::pureFluid; 

                    auto const& coefficients = pureFluid.at(component.name());

                    return Constants::gasConstantRatio*(
                            coefficients.polynomialCoefficient[1] + 
                            coefficients.polynomialCoefficient[2]*(temperature/component.criticalTemperature()) +
                            coefficients.polynomialCoefficient[3]*coefficients.hyperbolicCoefficient[0]/std::tanh(coefficients.hyperbolicCoefficient[0]*component.criticalTemperature()/temperature) + 
                            coefficients.polynomialCoefficient[5]*coefficients.hyperbolicCoefficient[2]/std::tanh(coefficients.hyperbolicCoefficient[2]*component.criticalTemperature()/temperature) - 
                            coefficients.polynomialCoefficient[4]*coefficients.hyperbolicCoefficient[1]*std::tanh(coefficients.hyperbolicCoefficient[1]*component.criticalTemperature()/temperature) -
                            coefficients.polynomialCoefficient[6]*coefficients.hyperbolicCoefficient[3]*std::tanh(coefficients.hyperbolicCoefficient[3]*component.criticalTemperature()/temperature)  
                        );
                }
            }

            // Second-Order Derivatives $\partial \delta^2$, $\partial \delta \partial \tau$
            if constexpr (wrt2 == WRT::reducedMixtureDensity){
                if constexpr (wrt1 == wrt2) {
                    return -std::pow(component.molarCriticalDensity()/mixtureDensity,2);
                }

                if constexpr (wrt1 == WRT::reducedTemperatureInverse) return 0;
            }

            if constexpr (wrt2 == WRT::reducedTemperatureInverse && wrt1 == wrt2){
                using Coefficients::PureFluid::Ideal::pureFluid; 

                auto const& coefficients = pureFluid.at(component.name());

                return 
                    Constants::gasConstantRatio*(
                        -coefficients.polynomialCoefficient[2]*std::pow(temperature/component.criticalTemperature(),2) -
                        coefficients.polynomialCoefficient[3]*std::pow(coefficients.hyperbolicCoefficient[0]/std::sinh(coefficients.hyperbolicCoefficient[0]*component.criticalTemperature()/temperature),2) - 
                        coefficients.polynomialCoefficient[5]*std::pow(coefficients.hyperbolicCoefficient[2]/std::sinh(coefficients.hyperbolicCoefficient[2]*component.criticalTemperature()/temperature),2) - 
                        coefficients.polynomialCoefficient[4]*std::pow(coefficients.hyperbolicCoefficient[1]/std::cosh(coefficients.hyperbolicCoefficient[1]*component.criticalTemperature()/temperature),2) -
                        coefficients.polynomialCoefficient[6]*std::pow(coefficients.hyperbolicCoefficient[3]/std::cosh(coefficients.hyperbolicCoefficient[3]*component.criticalTemperature()/temperature),2)  
                    );
                
            }
        }

        // WRT::None means first order derivative
        // This templated function accounts for derivatives with respect to delta and tau 
        template <WRT wrt1, WRT wrt2 = WRT::None>
        inline NP_t residualReducedHelmholtzFreeEnergyDerivative(Mixture const& mixture, NP_t const& reducedMixtureDensity, NP_t const& reducedTemperatureInverse, std::string compositionType = "global"){
            NP_t residualReducedHelmholtzEnergyDerivative_ = 0;
            for (auto const& component : mixture){
                residualReducedHelmholtzEnergyDerivative_ += component.composition(compositionType)*
                                                    residualReducedHelmholtzFreeEnergyDerivative<wrt1, wrt2>(component.pure(), reducedMixtureDensity, reducedTemperatureInverse);
            }

            for(unsigned int i=0; i<mixture.size()-1; ++i){
                for (unsigned int j=i+1; j<mixture.size(); ++j){
                     residualReducedHelmholtzEnergyDerivative_ += mixture[i].composition(compositionType)*mixture[j].composition(compositionType)*
                                                        residualDepartureDerivative<wrt1, wrt2>(mixture[i].pure(), mixture[j].pure(), reducedMixtureDensity, reducedTemperatureInverse);
                }
            }

            return residualReducedHelmholtzEnergyDerivative_;
        }

        template <WRT wrt1, WRT wrt2 = WRT::None>
        inline NP_t residualReducedHelmholtzFreeEnergyDerivative(Component const& component, NP_t const& reducedMixtureDensity, NP_t const& reducedTemperatureInverse) const {
            using namespace Coefficients::PureFluid::Residual;

            NP_t residualReducedHelmholtzFreeEnergyDerivative_ = 0;

            auto const& coefficients = pureFluid.at(component.name());

            auto const kPol = coefficients.polynomialExponent.size();

            for (unsigned int term = 0; term < kPol; ++term){

                auto const [d_oi, t_oi] = coefficients.polynomialExponent[term];

                // First-Order w.r.t delta
                if constexpr (wrt1 == WRT::reducedMixtureDensity && wrt2 == WRT::None){
                    residualReducedHelmholtzFreeEnergyDerivative_ += 
                        coefficients.residualCoefficient[term]*d_oi*std::pow(reducedMixtureDensity, d_oi-1)*
                        std::pow(reducedTemperatureInverse, t_oi);
                }
                // Second-order w.r.t delta^2
                if constexpr (wrt1 == WRT::reducedMixtureDensity && wrt2 == wrt1){
                    residualReducedHelmholtzFreeEnergyDerivative_ += 
                        coefficients.residualCoefficient[term]*d_oi*(d_oi-1)*std::pow(reducedMixtureDensity, d_oi-2)*
                        std::pow(reducedTemperatureInverse, t_oi);
                }
                //Second-order Mixed Derivative \partial delta \partial tau
                if constexpr ((wrt1 == WRT::reducedMixtureDensity && wrt2 == WRT::reducedTemperatureInverse) || (wrt2 == WRT::reducedMixtureDensity && wrt1 == WRT::reducedTemperatureInverse)){
                    residualReducedHelmholtzFreeEnergyDerivative_ += 
                        coefficients.residualCoefficient[term]*d_oi*t_oi*std::pow(reducedMixtureDensity, d_oi-1)*
                        std::pow(reducedTemperatureInverse, t_oi-1);
                }
                // First-Order w.r.t tau
                if constexpr (wrt1 == WRT::reducedTemperatureInverse && wrt2 == WRT::None){
                    residualReducedHelmholtzFreeEnergyDerivative_ += 
                        coefficients.residualCoefficient[term]*t_oi*std::pow(reducedMixtureDensity, d_oi)*
                        std::pow(reducedTemperatureInverse, t_oi - 1);
                }
                // Second-Order w.r.t tau^2
                if constexpr (wrt1 == WRT::reducedTemperatureInverse && wrt2 == wrt1){
                    residualReducedHelmholtzFreeEnergyDerivative_ += 
                        coefficients.residualCoefficient[term]*t_oi*std::pow(reducedMixtureDensity, d_oi)*
                        std::pow(reducedTemperatureInverse, t_oi - 1);
                }
            }

            auto const kExp = coefficients.exponentialExponent.size();

            for (unsigned int term = 0; term < kExp; ++term){

                auto const [d_oi, t_oi, c_oi] = coefficients.exponentialExponent[term];
                // First-Order w.r.t delta
                if constexpr (wrt1 == WRT::reducedMixtureDensity && wrt2 == WRT::None){
                    residualReducedHelmholtzFreeEnergyDerivative_ += 
                        coefficients.residualCoefficient[kPol + term]*
                        std::pow(reducedMixtureDensity, d_oi-1)*(d_oi - c_oi*std::pow(reducedMixtureDensity, c_oi))*
                        std::pow(reducedTemperatureInverse, t_oi)*std::exp(-std::pow(reducedMixtureDensity, c_oi));
                }
                // Second-order w.r.t delta^2
                if constexpr (wrt1 == WRT::reducedMixtureDensity && wrt2 == WRT::reducedMixtureDensity){
                        residualReducedHelmholtzFreeEnergyDerivative_ += 
                            coefficients.residualCoefficient[kPol + term]*std::pow(reducedMixtureDensity, d_oi-2)*
                            ((d_oi - c_oi*std::pow(reducedMixtureDensity, c_oi))*(d_oi - 1 - c_oi*std::pow(reducedMixtureDensity, c_oi)) - std::pow(c_oi,2)*std::pow(reducedMixtureDensity, c_oi))*
                            std::pow(reducedTemperatureInverse, t_oi)*std::exp(-std::pow(reducedMixtureDensity, c_oi));
                }
                //Second-order Mixed Derivative \partial delta \partial tau
                if constexpr ((wrt1 == WRT::reducedMixtureDensity && wrt2 == WRT::reducedTemperatureInverse) || (wrt2 == WRT::reducedMixtureDensity && wrt1 == WRT::reducedTemperatureInverse)){
                    residualReducedHelmholtzFreeEnergyDerivative_ += 
                        coefficients.residualCoefficient[kPol + term]*t_oi*
                        std::pow(reducedMixtureDensity, d_oi-1)*(d_oi - c_oi*std::pow(reducedMixtureDensity, c_oi))*
                        std::pow(reducedTemperatureInverse, t_oi-1)*std::exp(-std::pow(reducedMixtureDensity, c_oi));
                }
                // First-Order w.r.t tau
                if constexpr (wrt1 == WRT::reducedTemperatureInverse && wrt2 == WRT::None){
                    residualReducedHelmholtzFreeEnergyDerivative_ += 
                        coefficients.residualCoefficient[kPol + term]*t_oi*
                        std::pow(reducedMixtureDensity, d_oi)*std::pow(reducedTemperatureInverse, t_oi-1)*
                        std::exp(-std::pow(reducedMixtureDensity, c_oi));
                }
                // Second-Order w.r.t tau^2
                if constexpr (wrt1 == WRT::reducedTemperatureInverse && wrt2 == wrt1){
                    residualReducedHelmholtzFreeEnergyDerivative_ += 
                        coefficients.residualCoefficient[kPol + term]*t_oi*(t_oi - 1)*
                        std::pow(reducedMixtureDensity, d_oi)*std::pow(reducedTemperatureInverse, t_oi-2)*
                        std::exp(-std::pow(reducedMixtureDensity, c_oi));
                }
            }

            return residualReducedHelmholtzFreeEnergyDerivative_;
        }

        template <WRT wrt1, WRT wrt2 = WRT::None>
        inline NP_t residualDepartureDerivative(Component const& component1, Component const& component2, NP_t const& reducedMixtureDensity, NP_t const& reducedTemperatureInverse) const{

            NP_t residualDepartureDerivative_ = 0;

            auto const F = Coefficients::Departure::fParameter(component1.name(),component2.name());

            if (F != 0){

                auto const& coefficients = Coefficients::Departure::binaryDeparture.at({component1.name(), component2.name()});

                auto const kPol = coefficients.polynomialExponent.size();

                for (unsigned int term = 0; term < kPol; ++term){

                    auto const [d_ij, t_ij] = coefficients.polynomialExponent[term];

                    // First-Order w.r.t delta
                    if constexpr (wrt1 == WRT::reducedMixtureDensity && wrt2 == WRT::None){
                        residualDepartureDerivative_ += 
                            coefficients.coefficient[term]*d_ij*std::pow(reducedMixtureDensity, d_ij-1)*std::pow(reducedTemperatureInverse, t_ij);
                    }
                    // Second-order w.r.t delta^2
                    if constexpr (wrt1 == WRT::reducedMixtureDensity && wrt2 == WRT::reducedMixtureDensity){
                        residualDepartureDerivative_ += 
                            coefficients.coefficient[term]*d_ij*(d_ij-1)*std::pow(reducedMixtureDensity, d_ij-2)*std::pow(reducedTemperatureInverse, t_ij);
                    }
                    //Second-order Mixed Derivative \partial delta \partial tau
                    if constexpr ((wrt1 == WRT::reducedMixtureDensity && wrt2 == WRT::reducedTemperatureInverse) || (wrt2 == WRT::reducedMixtureDensity && wrt1 == WRT::reducedTemperatureInverse)){
                        residualDepartureDerivative_ += 
                            coefficients.coefficient[term]*d_ij*t_ij*std::pow(reducedMixtureDensity, d_ij-1)*std::pow(reducedTemperatureInverse, t_ij-1);
                    }
                    // First-Order w.r.t tau
                    if constexpr (wrt1 == WRT::reducedTemperatureInverse && wrt2 == WRT::None){
                        residualDepartureDerivative_ += 
                            coefficients.coefficient[term]*t_ij*std::pow(reducedMixtureDensity, d_ij)*std::pow(reducedTemperatureInverse, t_ij-1);
                    }
                    // Second-Order w.r.t tau^2
                    if constexpr (wrt1 == WRT::reducedTemperatureInverse && wrt2 == wrt1){
                        residualDepartureDerivative_ += 
                            coefficients.coefficient[term]*t_ij*(t_ij-1)*std::pow(reducedMixtureDensity, d_ij)*std::pow(reducedTemperatureInverse, t_ij-2);
                    }
                }

                auto const kExp = coefficients.exponentialExponent.size();

                for (unsigned int term = 0; term < kExp; ++term){

                    auto const [d_ij, t_ij, eta_ij, varepsilon_ij, beta_ij, gamma_ij] = coefficients.exponentialExponent[term];

                     // First-Order w.r.t delta
                    if constexpr (wrt1 == WRT::reducedMixtureDensity && wrt2 == WRT::None){
                        residualDepartureDerivative_ += 
                            coefficients.coefficient[kPol + term]*std::pow(reducedMixtureDensity, d_ij)*
                            std::pow(reducedTemperatureInverse, t_ij)*std::exp(
                                -eta_ij*std::pow(reducedMixtureDensity - varepsilon_ij, 2) - beta_ij*(reducedMixtureDensity - gamma_ij)
                            )*(d_ij/reducedMixtureDensity -2*eta_ij*(reducedMixtureDensity - varepsilon_ij) - beta_ij);
                    }
                    // Second-order w.r.t delta^2
                    if constexpr (wrt1 == WRT::reducedMixtureDensity && wrt2 == WRT::reducedMixtureDensity){
                        residualDepartureDerivative_ += 
                            coefficients.coefficient[kPol + term]*std::pow(reducedMixtureDensity, d_ij)*
                            std::pow(reducedTemperatureInverse, t_ij)*std::exp(
                                -eta_ij*std::pow(reducedMixtureDensity - varepsilon_ij, 2) - beta_ij*(reducedMixtureDensity - gamma_ij)
                            )*(std::pow(d_ij/reducedMixtureDensity -2*eta_ij*(reducedMixtureDensity - varepsilon_ij) - beta_ij,2) - 
                            d_ij/std::pow(reducedMixtureDensity, 2) - 2*eta_ij);
                    }
                    //Second-order Mixed Derivative \partial delta \partial tau
                    if constexpr ((wrt1 == WRT::reducedMixtureDensity && wrt2 == WRT::reducedTemperatureInverse) || (wrt2 == WRT::reducedMixtureDensity && wrt1 == WRT::reducedTemperatureInverse)){
                        residualDepartureDerivative_ += 
                            coefficients.coefficient[kPol + term]*t_ij*std::pow(reducedMixtureDensity, d_ij)*
                            std::pow(reducedTemperatureInverse, t_ij-1)*std::exp(
                                -eta_ij*std::pow(reducedMixtureDensity - varepsilon_ij, 2) - beta_ij*(reducedMixtureDensity - gamma_ij)
                            )*(d_ij/reducedMixtureDensity -2*eta_ij*(reducedMixtureDensity - varepsilon_ij) - beta_ij);
                    }
                    // First-Order w.r.t tau
                    if constexpr (wrt1 == WRT::reducedTemperatureInverse && wrt2 == WRT::None){
                        residualDepartureDerivative_ += 
                           coefficients.coefficient[kPol + term]*t_ij*std::pow(reducedMixtureDensity, d_ij)*
                           std::pow(reducedTemperatureInverse, t_ij-1)*std::exp(
                               -eta_ij*std::pow(reducedMixtureDensity - varepsilon_ij, 2) - beta_ij*(reducedMixtureDensity - gamma_ij)
                           );
                    }
                    // Second-Order w.r.t tau^2
                    if constexpr (wrt1 == WRT::reducedTemperatureInverse && wrt2 == wrt1){
                        residualDepartureDerivative_ += 
                           coefficients.coefficient[kPol + term]*t_ij*(t_ij-1)*std::pow(reducedMixtureDensity, d_ij)*
                           std::pow(reducedTemperatureInverse, t_ij-2)*std::exp(
                               -eta_ij*std::pow(reducedMixtureDensity - varepsilon_ij, 2) - beta_ij*(reducedMixtureDensity - gamma_ij)
                           );
                    }
                }
            }

            return F*residualDepartureDerivative_;
        }

        inline std::tuple<NP_t, NP_t> pressure(Mixture const& mixture, NP_t const& mixtureDensity, NP_t const& temperature){

            NP_t const delta = reducedMixtureDensity(mixtureDensity);
            NP_t const tau = reducedTemperatureInverse(temperature);

            auto const partialHelmholtzPartialDelta = residualReducedHelmholtzFreeEnergyDerivative<WRT::reducedMixtureDensity>(mixture, delta, tau);
            auto const partial2HelmholtzPartialDelta2 = residualReducedHelmholtzFreeEnergyDerivative<WRT::reducedMixtureDensity, WRT::reducedMixtureDensity>(mixture, delta, tau);

            compressibilityFactor_ = 1 + delta*partialHelmholtzPartialDelta;

            auto newPressure = mixtureDensity*Constants::universalGasesConstant*temperature*compressibilityFactor_;

            auto pressureDerivativeVSDensity = temperature*Constants::universalGasesConstant*(
                1 + 2*delta*partialHelmholtzPartialDelta + std::pow(delta, 2)*partial2HelmholtzPartialDelta2
            );

            return {newPressure, pressureDerivativeVSDensity};

        }

        public:

        constexpr inline std::string_view name(){
            return "GERG2008";
        }

        inline NP_t operator()(Mixture const& mixture, NP_t const& specifiedPressure, NP_t const& temperature){

            NP_t pseudoCriticalDensity = mixture.pseudoCriticalDensity();
            NP_t pseudoCriticalTemperature = mixture.pseudoCriticalTemperature();

            reducingDensity(mixture);
            reducingTemperature(mixture);

            NP_t idealDensity = specifiedPressure/(Constants::universalGasesConstant*temperature); // Vapor density estimation based in ideal gas

            NP_t logVolume = -std::log(idealDensity);
            NP_t logSpecifiedPressure = std::log(specifiedPressure);

            NP_t tolerance = 1e-7;

            bool error = false;

            //Translation from AGA 8 NIST code
            for (unsigned int it = 0; it < 50; ++it){
                NP_t density = std::exp(-logVolume);
                auto [calculatedPressure, pressureDerivative] = pressure(mixture, density, temperature);
                NP_t pressureDerivativeAgainstLogVolume = -density * pressureDerivative;
                NP_t logVolumeDifference = (std::log(calculatedPressure) - logSpecifiedPressure)*calculatedPressure/pressureDerivativeAgainstLogVolume;
                logVolume -= logVolumeDifference;
                if(std::abs(logVolumeDifference) < tolerance){
                    if (pressureDerivative > 0){
                        return std::exp(-logVolume);
                    }
                }
            }
            std::cerr << "Error in GERG2008 at pressure: " << specifiedPressure << ". Ideal gas density returned." << std::endl;
            return idealDensity;
        }

        inline NP_t selectedCompressibility() const {
            return compressibilityFactor_;
        }

    };

}
#endif /* GERG2008_HPP */