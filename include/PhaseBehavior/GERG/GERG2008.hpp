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
            NP_t reducingMixtureDensity_, reducingMixtureTemperature_;


            inline NP_t reducingDensity(Mixture const& mixture, std::string compositionType = "global"){

                using Coefficients::ReducingParameters::Density::beta, Coefficients::ReducingParameters::Density::gamma;

                reducingMixtureDensity_ = 0;

                for(auto const& component: mixture){
                    reducingMixtureDensity_ += std::pow(component.composition(compositionType),2)/component->criticalDensity();
                }

                for(unsigned int i=0; i<mixture.size()-1; ++i){
                    for (unsigned int j=i+1; j<mixture.size(); ++j){
                        reducingMixtureDensity_ += 2*mixture[i].composition(compositionType)*mixture[j].composition(compositionType)*
                        beta.at({mixture[i]->name(), mixture[j]->name()})*gamma.at({mixture[i]->name(), mixture[j]->name()})*
                        (mixture[i].composition(compositionType) + mixture[j].composition(compositionType))/
                        (std::pow(beta.at({mixture[i]->name(), mixture[j]->name()}),2)*mixture[i].composition(compositionType) + mixture[j].composition(compositionType))*
                        (1.0/8.0)*std::pow(1/std::cbrt(mixture[i]->criticalDensity()) + 1/std::cbrt(mixture[j]->criticalDensity()),3);
                    }
                }
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
            }

            inline NP_t reducedTemperatureInverse(NP_t const& temperature){
                return reducingMixtureTemperature_/temperature;
            }


            inline NP_t reducedMixtureDensity(NP_t const& density) const{
                return density*reducingMixtureDensity_;
            }

        public:

        inline NP_t idealReducedHelmholtzFreeEnergy(Component const& component, NP_t const& mixtureDensity, NP_t const& temperature) const {

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
            NP_t idealMixtureReducedFreeEnergy = 0;
            for (auto& component : mixture){
                idealMixtureReducedFreeEnergy += component.composition(compositionType)*idealReducedHelmholtzFreeEnergy(component.pure(), mixtureDensity, temperature) +
                              std::log(component.composition(compositionType));
            }
            return idealMixtureReducedFreeEnergy;
        }

        inline NP_t residualHelmholtzFreeEnergy(Component const& component, NP_t const& reducedMixtureDensity, NP_t const& reducedTemperatureInverse) const {
            using namespace Coefficients::PureFluid::Residual;

            NP_t pureFluidResidualHelmholtzFreeEnergy = 0;

            auto const& coefficients = pureFluid.at(component.name());

            const auto kPol = coefficients.polynomialExponent.size();

            for (unsigned int term = 0; term < kPol; ++term){
                pureFluidResidualHelmholtzFreeEnergy += 
                    coefficients.residualCoefficient[term]*
                    std::pow(reducedMixtureDensity, std::get<reducedDensityExponent>(coefficients.polynomialExponent[term]))*
                    std::pow(reducedTemperatureInverse, std::get<reducedTemperatureExponent>(coefficients.polynomialExponent[term]));
            }

            const auto kExp = coefficients.exponentialExponent.size();

            for (unsigned int term = 0; term < kExp; ++term){
                pureFluidResidualHelmholtzFreeEnergy += 
                    coefficients.residualCoefficient[kPol + term]*
                    std::pow(reducedMixtureDensity, std::get<reducedDensityExponent>(coefficients.exponentialExponent[term]))*
                    std::pow(reducedTemperatureInverse, std::get<reducedTemperatureExponent>(coefficients.exponentialExponent[term]))*
                    std::exp(-std::pow(reducedMixtureDensity, std::get<exponentialReducedDensityExponent>(coefficients.exponentialExponent[term])));
            }
        };

        inline NP_t residualDeparture(Component const& component1, Component const& component2, NP_t const& reducedMixtureDensity, NP_t const& reducedTemperatureInverse){
            using namespace Coefficients::Departure;

            NP_t residualDeparture = 0;

            auto F = fParameter(component1.name(),component2.name());

            if (F != 0){

                auto const& coefficients = binaryDeparture.at({component1.name(), component2.name()});

                const auto kPol = coefficients.polynomialExponent.size();

                for (unsigned int term = 0; term < kPol; ++term){
                    residualDeparture += 
                        coefficients.coefficient[term]*
                        std::pow(reducedMixtureDensity, std::get<reducedDensityExponent>(coefficients.polynomialExponent[term]))*
                        std::pow(reducedTemperatureInverse, std::get<reducedTemperatureExponent>(coefficients.polynomialExponent[term]));
                }

                const auto kExp = coefficients.exponentialExponent.size();

                for (unsigned int term = 0; term < kExp; ++term){
                    residualDeparture += 
                        coefficients.coefficient[kPol + term]*
                        std::pow(reducedMixtureDensity, std::get<reducedDensityExponent>(coefficients.exponentialExponent[term]))*
                        std::pow(reducedTemperatureInverse, std::get<reducedTemperatureExponent>(coefficients.exponentialExponent[term]))*
                        std::exp(
                            -std::get<eta>(coefficients.exponentialExponent[term])*std::pow(reducedMixtureDensity - std::get<varepsilon>(coefficients.exponentialExponent[term]) , 2) - 
                            std::get<beta>(coefficients.exponentialExponent[term])*(reducedMixtureDensity - std::get<gamma>(coefficients.exponentialExponent[term]))
                        );
                }
            }

            return residualDeparture*F;
        }
    };

}
#endif /* GERG2008_HPP */