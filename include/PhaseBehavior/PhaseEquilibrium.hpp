#ifndef PHASE_EQUILIBRIUM_HPP
#define PHASE_EQUILIBRIUM_HPP

#include <algorithm>
#include <tuple>

#include "Utilities/Math.hpp"
#include "Mixture.hpp"

using NP_t = PhaseBehavior::Types::NumericalPrecision;

namespace PhaseBehavior::VaporLiquidEquilibrium {

    void rachfordVLE(Mixture& mixture){

        auto [min, max] = std::minmax_element( mixture.begin(), mixture.end(), [](const auto& lhs, const auto& rhs) {
            return lhs.equilibriumCoefficient() < rhs.equilibriumCoefficient();
        } );

        auto Kimax = max->equilibriumCoefficient();
        auto Kimin = min->equilibriumCoefficient();

        auto lowerLimit = 1.0/(1.0-Kimax + 1e-10);
        auto upperLimit = 1.0/(1.0-Kimin - 1e-10);

        auto objectiveFunction = [&mixture](const NP_t& vaporMolarFraction){
            NP_t sum = 0;
            for (const auto& mixComp : mixture){
                const auto Ki = mixComp.equilibriumCoefficient();
                sum += mixComp.composition()*(Ki-1.0)/(1.0+vaporMolarFraction*(Ki-1.0));
            }
            return sum;
        };

        auto derivativeFunction = [&mixture](const NP_t& vaporMolarFraction){
            NP_t sum = 0;
            for (const auto& mixComp : mixture){
                const auto Ki = mixComp.equilibriumCoefficient();
                sum -= mixComp.composition()*std::pow(Ki-1.0,2)/std::pow(1.0+vaporMolarFraction*(Ki-1.0),2);
            }
            return sum;
        };

        auto divergence = [&lowerLimit, &upperLimit](const NP_t& vaporMolarFraction){
            return !(lowerLimit < vaporMolarFraction && vaporMolarFraction < upperLimit);
        };

        auto hasRoot = objectiveFunction(0) > 0 && objectiveFunction(1) < 0;
        if (hasRoot){
            NP_t initialGuess = (lowerLimit + upperLimit)/2.0;
        
            auto possibleFng = Math::NewtonRaphson<NP_t>(initialGuess, divergence, objectiveFunction, derivativeFunction);

            NP_t vaporMolarFraction;
            if (possibleFng){
                vaporMolarFraction = possibleFng.value();
            }else{
                vaporMolarFraction = Math::bisection<NP_t>(0.0, 1.0, objectiveFunction);
            }

            auto liquidMolarFraction = 1 - vaporMolarFraction;

            mixture.molarFraction("vapor", vaporMolarFraction);
            mixture.molarFraction("liquid", liquidMolarFraction);

            for (auto& mixComp : mixture){
                auto Ki = mixComp.equilibriumCoefficient();
                auto liquidComposition = mixComp.composition()/(1.0 + vaporMolarFraction*(Ki-1));
                auto vaporComposition  = Ki*liquidComposition;
                mixComp.composition("liquid", liquidComposition);
                mixComp.composition("vapor", vaporComposition);
            }
        }
    }

    template<typename EoS>
    std::tuple<EoS, EoS> succesiveSubstitution(Mixture& mixture, NP_t const& pressure, NP_t const& temperature, bool restartKi=true){

        if(restartKi) mixture.initializeEquilibriumCoefficients(pressure, temperature); //Wilson Correlation

        EoS vaporEoS, liquidEoS;

        NP_t fugacitiesSquaredSum = 10;
        while(fugacitiesSquaredSum >= 1e-14){

            rachfordVLE(mixture);

            vaporEoS(mixture, pressure, temperature, "vapor");
            mixture.compressibility("vapor", vaporEoS.selectedCompressibility());
            vaporEoS.fugacities(mixture, "vapor");

            liquidEoS(mixture, pressure, temperature, "liquid");
            mixture.compressibility("liquid", liquidEoS.selectedCompressibility());
            liquidEoS.fugacities(mixture, "liquid");

            for (auto& mixComponent : mixture){
                auto newEquilibriumCoefficient = mixComponent.equilibriumCoefficient()*(mixComponent.fugacity("liquid", pressure)/mixComponent.fugacity("vapor", pressure));
                mixComponent.equilibriumCoefficient(newEquilibriumCoefficient);
            }

            fugacitiesSquaredSum = std::accumulate(mixture.begin(),mixture.end(), static_cast<NP_t>(0.0),
            [&pressure](auto previous, auto second){
                return previous + std::pow(second.fugacity("liquid", pressure)/second.fugacity("vapor", pressure) - 1.0, 2);
            });
        }

        return {vaporEoS, liquidEoS};
    }

    enum class PhaseStabilityTestResult {
        TrivialSolution,
        LessThanOne,
        GreaterThanOne
    };

    enum class PhaseStabilityResult {
        Stable,
        Unstable
    };

    template<typename EoS>
    PhaseStabilityResult phaseStability(Mixture& mixture, NP_t const& pressure, NP_t const& temperature) noexcept{

        constexpr auto liquidMolesCalculation = [](NP_t composition, NP_t equilibriumCoefficient){
            return composition/equilibriumCoefficient;
        };

        constexpr auto vaporMolesCalculation = [](NP_t composition, NP_t equilibriumCoefficient){
            return composition*equilibriumCoefficient;
        };

        constexpr auto liquidCorrectionFunction = [](NP_t globalFugacity, NP_t liquidFugacity, NP_t totalMoles){
            return (liquidFugacity/globalFugacity)*totalMoles;
        };

        constexpr auto vaporCorrectionFunction = [](NP_t globalFugacity, NP_t vaporFugacity, NP_t totalMoles){
            return (globalFugacity/vaporFugacity)/totalMoles;
        };

        EoS eos;

        eos(mixture, pressure, temperature);
        mixture.compressibility("global", eos.selectedCompressibility());
        eos.fugacities(mixture, "global");

        auto test = [&mixture, &pressure, &temperature](std::string const& phaseName, EoS& eos, auto const& molesCalculation, auto const& correctionFunction){

            mixture.initializeEquilibriumCoefficients(pressure, temperature); //Wilson Correlation

            NP_t correctionSquaredValues = 10; 

            while(correctionSquaredValues >= 1e-10){
                auto totalMoles = std::accumulate(mixture.begin(), mixture.end(), static_cast<NP_t>(0),
                [&phaseName, &molesCalculation](auto previous, auto const& element){
                    return previous + molesCalculation(element.composition(), element.equilibriumCoefficient());
                });

                for (auto& mixComp : mixture){
                    auto moles = molesCalculation(mixComp.composition(), mixComp.equilibriumCoefficient());
                    mixComp.composition(phaseName, moles/totalMoles);
                }

                eos(mixture, pressure, temperature, phaseName);
                mixture.compressibility(phaseName, eos.selectedCompressibility());
                eos.fugacities(mixture, phaseName);

                correctionSquaredValues = 0;
                for (auto& mixComp : mixture){
                    auto correction = correctionFunction(mixComp.fugacity("global", pressure), mixComp.fugacity(phaseName, pressure), totalMoles);
                    auto newEquilibriumCoefficient = mixComp.equilibriumCoefficient()*correction;
                    mixComp.equilibriumCoefficient(newEquilibriumCoefficient);
                    correctionSquaredValues += std::pow(correction - 1, 2.0);
                }

                auto sumOfLogsOfKiSquared = std::accumulate(mixture.begin(), mixture.end(), static_cast<NP_t>(0),
                [](auto previous, auto const& element){
                    return previous + std::pow(std::log(element.equilibriumCoefficient()),2.0);
                });

                if (sumOfLogsOfKiSquared < 1e-4) return PhaseStabilityTestResult::TrivialSolution;
            }
            auto totalMoles = std::accumulate(mixture.begin(), mixture.end(), static_cast<NP_t>(0),
                [&phaseName, &molesCalculation](auto previous, auto const& element){
                    return previous + molesCalculation(element.composition(), element.equilibriumCoefficient());
                });

            return totalMoles > 1 ? PhaseStabilityTestResult::GreaterThanOne : PhaseStabilityTestResult::LessThanOne;
        };

        auto vaporTestResult = test("vapor", eos, vaporMolesCalculation, vaporCorrectionFunction);
        if(vaporTestResult == PhaseStabilityTestResult::GreaterThanOne) return PhaseStabilityResult::Unstable;
        auto liquidTestResult = test("liquid", eos, liquidMolesCalculation, liquidCorrectionFunction);
        if (liquidTestResult == PhaseStabilityTestResult::GreaterThanOne) return PhaseStabilityResult::Unstable;
        return PhaseStabilityResult::Stable;
    }

}

#endif /* PHASE_EQUILIBRIUM_HPP */