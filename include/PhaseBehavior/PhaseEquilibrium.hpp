#ifndef PHASE_EQUILIBRIUM_HPP
#define PHASE_EQUILIBRIUM_HPP

#include <algorithm>
#include <tuple>
#include <vector>

#include <Eigen/Dense>

#include "Utilities/Math.hpp"
#include "Mixture.hpp"
#include "Phase.hpp"

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

            NP_t vaporMolarFraction = 0;
            if (possibleFng){
                vaporMolarFraction = possibleFng.value();
            }else{
                auto possibleBisectionValue = Math::bisection<NP_t>(0.0, 1.0, objectiveFunction);
                if (possibleBisectionValue) vaporMolarFraction = possibleBisectionValue.value();
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
    std::tuple<EoS, EoS> succesiveSubstitution(Mixture& mixture, NP_t const& pressure, NP_t const& temperature, bool restartKi=true, bool useNewton=false, NP_t switchValue=1e-14){

        static constexpr NP_t machineEpsilon = std::sqrt(std::numeric_limits<NP_t>::epsilon());

        if(restartKi) mixture.initializeEquilibriumCoefficients(pressure, temperature); //Wilson Correlation

        EoS vaporEoS, liquidEoS;
        if (!useNewton) switchValue=1e-14;


        NP_t fugacitiesSquaredSum = 10;
        while(fugacitiesSquaredSum >= switchValue){

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

        if (useNewton){
            Eigen::MatrixXd jacobian(mixture.size(), mixture.size());
            Eigen::VectorXd residual, solutionDelta;

            residual.resize(mixture.size());
            solutionDelta.resize(mixture.size());
            jacobian.setZero();
            residual.setZero();
            solutionDelta.setZero();

            //**************************************************** NEWTON PROCEDURE ************************************************************
            while (fugacitiesSquaredSum >= 1e-14){

                //************************************************ RESIDUAL CALCULATION ********************************************************
                rachfordVLE(mixture);

                vaporEoS(mixture, pressure, temperature, "vapor");
                mixture.compressibility("vapor", vaporEoS.selectedCompressibility());
                vaporEoS.fugacities(mixture, "vapor");

                liquidEoS(mixture, pressure, temperature, "liquid");
                mixture.compressibility("liquid", liquidEoS.selectedCompressibility());
                liquidEoS.fugacities(mixture, "liquid");

                for (std::size_t i = 0; i < mixture.size(); ++i){

                    residual[i] = (std::log(mixture[i].composition("vapor")) + std::log(mixture[i].fugacityCoefficient("vapor"))) -
                                            (std::log(mixture[i].composition("liquid")) + std::log(mixture[i].fugacityCoefficient("liquid")));
                }

                //******************************************************************************************************************************

                //************************************************ JACOBIAN CALCULATION ********************************************************
                //************************************************ Component Derivatives *******************************************************
                for (std::size_t j = 0; j < mixture.size(); ++j){

                    auto originalki = mixture[j].equilibriumCoefficient();
                    auto scaledEpsilon = std::log(originalki)*machineEpsilon;
                    auto newKi = std::exp(std::log(originalki) + scaledEpsilon);
                    mixture[j].equilibriumCoefficient(newKi);

                    rachfordVLE(mixture);

                    vaporEoS(mixture, pressure, temperature, "vapor");
                    mixture.compressibility("vapor", vaporEoS.selectedCompressibility());
                    vaporEoS.fugacities(mixture, "vapor");

                    liquidEoS(mixture, pressure, temperature, "liquid");
                    mixture.compressibility("liquid", liquidEoS.selectedCompressibility());
                    liquidEoS.fugacities(mixture, "liquid");
                    
                    for (std::size_t i = 0; i < mixture.size(); ++i){
                        auto modifiedResidual = (std::log(mixture[i].composition("vapor")) + std::log(mixture[i].fugacityCoefficient("vapor"))) -
                                            (std::log(mixture[i].composition("liquid")) + std::log(mixture[i].fugacityCoefficient("liquid")));
                        jacobian(i,j) = (modifiedResidual - residual[i])/scaledEpsilon;
                    }
                    mixture[j].equilibriumCoefficient(originalki);
                }
                //******************************************************************************************************************************
                
                solutionDelta = jacobian.colPivHouseholderQr().solve(-residual);

                for (std::size_t j = 0; j < mixture.size(); ++j){
                    auto originalki = mixture[j].equilibriumCoefficient();
                    auto newKi = std::exp(std::log(originalki) + solutionDelta[j]);
                    mixture[j].equilibriumCoefficient(newKi);
                }

                rachfordVLE(mixture);

                vaporEoS(mixture, pressure, temperature, "vapor");
                mixture.compressibility("vapor", vaporEoS.selectedCompressibility());
                vaporEoS.fugacities(mixture, "vapor");

                liquidEoS(mixture, pressure, temperature, "liquid");
                mixture.compressibility("liquid", liquidEoS.selectedCompressibility());
                liquidEoS.fugacities(mixture, "liquid");

                fugacitiesSquaredSum = std::accumulate(mixture.begin(),mixture.end(), static_cast<NP_t>(0.0),
                [&pressure](auto previous, auto second){
                    return previous + std::pow(second.fugacity("liquid", pressure)/second.fugacity("vapor", pressure) - 1.0, 2);
                });
            }
        }

        return {vaporEoS, liquidEoS};
    }

    enum class PhaseStabilityTestResult {
        TrivialSolution,
        LessThanOne,
        GreaterThanOne
    };

    enum class PhaseStabilityResult {
        Stable = true,
        Unstable = false
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

    template<typename EoS>
    auto isothermalTwoPhaseFlash(Mixture& mixture, NP_t const& pressure, NP_t const& temperature){
        using Phase::phaseIdentification;
        using Phase = Phase::FluidPhase;
        if(phaseStability<EoS>(mixture, pressure, temperature) == PhaseStabilityResult::Unstable){
            auto [vaporEoS, liquidEoS] = succesiveSubstitution<EoS>(mixture, pressure, temperature, false, true, 1e-4);
            return PhaseStabilityResult::Unstable, std::vector<std::shared_ptr<Phase>>{
                phaseIdentification<EoS>(mixture, mixture.compressibility("vapor"),pressure, temperature, vaporEoS, "vapor"),
                phaseIdentification<EoS>(mixture, mixture.compressibility("liquid"), pressure, temperature, liquidEoS, "liquid")
            };
        }else{
            EoS eos;
            eos(mixture, pressure, temperature);
            mixture.compressibility("global", eos.selectedCompressibility());
            return PhaseStabilityResult::Stable, std::vector<std::shared_ptr<Phase>>{
                phaseIdentification<EoS>(mixture, mixture.compressibility("global"), pressure, temperature, eos)
            };
        }
    }
}

#endif /* PHASE_EQUILIBRIUM_HPP */