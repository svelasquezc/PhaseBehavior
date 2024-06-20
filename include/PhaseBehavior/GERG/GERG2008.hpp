#ifndef GERG2008_HPP
#define GERG2008_HPP
#include "../Utilities/Types.hpp"
#include "../Utilities/Math.hpp"
#include "../Component.hpp"
#include "Constants.hpp"
#include "../Mixture.hpp"

#include "PureFluidHelmholtzCoefficients.hpp"

using NP_t = PhaseBehavior::Types::NumericalPrecision;

namespace PhaseBehavior::EoS::GERG{

    using PhaseBehavior::Component;
    using PhaseBehavior::Mixture;

    class GERG2008{

        public:

        inline NP_t idealReducedHelmholtzFreeEnergy(Component const& component, NP_t const& mixtureDensity, NP_t const& temperature) const {

            using Coefficients::PureFluid::Ideal::idealCoefficient;
            using Coefficients::PureFluid::Ideal::hyperbolicCoefficient;

            return std::log(mixtureDensity/component.molarCriticalDensity()) + Constants::gasConstantRatio * (
                idealCoefficient.at(component.name())[0] + 
                idealCoefficient.at(component.name())[1]*component.criticalTemperature()/temperature + 
                idealCoefficient.at(component.name())[2]*std::log(component.criticalTemperature()/temperature) +
                idealCoefficient.at(component.name())[3]*std::log(std::abs(std::sinh(hyperbolicCoefficient.at(component.name())[0]*component.criticalTemperature()/temperature))) + 
                idealCoefficient.at(component.name())[5]*std::log(std::abs(std::sinh(hyperbolicCoefficient.at(component.name())[2]*component.criticalTemperature()/temperature))) - 
                idealCoefficient.at(component.name())[4]*std::log(std::abs(std::sinh(hyperbolicCoefficient.at(component.name())[1]*component.criticalTemperature()/temperature))) -
                idealCoefficient.at(component.name())[6]*std::log(std::abs(std::sinh(hyperbolicCoefficient.at(component.name())[3]*component.criticalTemperature()/temperature)))  
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

        inline NP_t residualHelmholtzFreeEnergy(Component const& component, NP_t const& mixtureReducedDensity, NP_t const& inverseReducedTemperature) const {
            using namespace Coefficients::PureFluid::Residual;

            NP_t pureFluidResidualHelmholtzFreeEnergy = 0;

            const auto kPol = polynomialExponent.at(component.name()).size();

            for (unsigned int term = 0; term < kPol; ++term){
                pureFluidResidualHelmholtzFreeEnergy += 
                    residualCoefficient.at(component.name())[term]*
                    std::pow(mixtureReducedDensity, std::get<reducedDensityExponent>(polynomialExponent.at(component.name())[term]))*
                    std::pow(inverseReducedTemperature, std::get<reducedTemperatureExponent>(polynomialExponent.at(component.name())[term]));
            }

            const auto kExp = exponentialExponent.at(component.name()).size();

            for (unsigned int term = 0; term < kExp; ++term){
                pureFluidResidualHelmholtzFreeEnergy += 
                    residualCoefficient.at(component.name())[kPol + term]*
                    std::pow(mixtureReducedDensity, std::get<reducedDensityExponent>(exponentialExponent.at(component.name())[term]))*
                    std::pow(inverseReducedTemperature, std::get<reducedTemperatureExponent>(exponentialExponent.at(component.name())[term]))*
                    std::exp(-std::pow(mixtureReducedDensity, std::get<exponentialReducedDensityExponent>(exponentialExponent.at(component.name())[term])));
            }
        };

    };

}
#endif /* GERG2008_HPP */