#ifndef PHASE_ENVELOPE_HPP
#define PHASE_ENVELOPE_HPP

#include <tuple>
#include <vector>

#include <Eigen/Sparse>

#include "Utilities/Types.hpp"
#include "Utilities/Math.hpp"
#include "Mixture.hpp"
#include "PhaseEquilibrium.hpp"

using NP_t = PhaseBehavior::Types::NumericalPrecision;
using PhaseBehavior::Types::FixedQueue;

namespace PhaseBehavior{

    template<typename EoS>
    class PhaseEnvelope{
    private:
        EoS vaporEoS, liquidEoS;
        Mixture& mixture_;
        NP_t vaporMolarFraction_;

        Eigen::SparseMatrix<double> jacobian_;
        std::vector<Eigen::Triplet<double>> nonZeros_;
        Eigen::VectorXd residualVector_, solutionDelta_;

        Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double, int>> solver_;
        std::size_t size_ = 0;
        NP_t currentPressure_ = 0;
        NP_t currentTemperature_ = -460;
        FixedQueue<Eigen::VectorXd, 2> solutions_;
        static constexpr NP_t machineEpsilon = std::sqrt(std::numeric_limits<NP_t>::epsilon());

        struct EnvelopePoint {
            NP_t pressure;
            NP_t bubbleTemperature;
            NP_t dewTemperature;
        };

        using EnvelopePoints = std::vector<EnvelopePoint>;

    public:

        PhaseEnvelope(Mixture& mixture): mixture_(mixture){

            size_ = mixture.size() + 1;

            jacobian_.resize(size_, size_);
            jacobian_.reserve(3*mixture.size() + 1);
            nonZeros_.reserve(3*mixture.size() + 1);
            residualVector_.resize(size_);
            solutionDelta_.resize(size_);
    
            jacobian_.setZero();
            residualVector_.setZero();
            solutionDelta_.setZero();
        }

        std::pair<NP_t, NP_t> nextPoint(NP_t const& vaporMolarFraction){

            solutionDelta_.array() = 10; 

            if (!solutions_.empty()){
                if(solutions_.size()==4){

                }
            }


            //**************************************************** NEWTON PROCEDURE ************************************************************
            while (solutionDelta_.squaredNorm() > 1e-10){

                //************************************************ RESIDUAL CALCULATION ********************************************************
                for (std::size_t i = 0; i < mixture_.size(); ++i){
                    auto ki = mixture_[i].equilibriumCoefficient();
                    auto xi = mixture_[i].composition()/(1+vaporMolarFraction*(ki - 1));
                    auto yi = ki*mixture_[i].composition()/(1+vaporMolarFraction*(ki - 1));
                    mixture_[i].composition("vapor", yi);
                    mixture_[i].composition("liquid", xi);
                }

                vaporEoS(mixture_,currentPressure_, currentTemperature_, "vapor");
                vaporEoS.fugacities(mixture_,currentPressure_, currentTemperature_, "vapor");
                liquidEoS(mixture_,currentPressure_, currentTemperature_, "liquid");
                liquidEoS.fugacities(mixture_,currentPressure_, currentTemperature_, "liquid");

                for (std::size_t i = 0; i < mixture_.size(); ++i){

                    residualVector_[i] = std::log(mixture_[i].equilibriumCoefficient()) + std::log(mixture_[i].fugacityCoefficient("vapor")) -
                                            std::log(mixture_[i].fugacityCoefficient("liquid"));
                }

                residualVector_[size_-1] = std::accumulate(mixture_.begin(), mixture_.end(), static_cast<NP_t>(0), 
                [F = vaporMolarFraction](auto previous, auto const& mixComp){
                    auto ki = mixComp.equilibriumCoefficient();
                    auto zi = mixComp.composition();
                    return previous + (ki-1)*zi/(1+F*(ki-1));
                });

                //******************************************************************************************************************************

                //************************************************ JACOBIAN CALCULATION ********************************************************
                //************************************************ Component Derivatives *******************************************************
                for (std::size_t i=0; i<mixture_.size(); ++i){
                    const auto originalki = mixture_[i].equilibriumCoefficient();
                    const auto originalyi = mixture_[i].composition("vapor");
                    const auto originalxi = mixture_[i].composition("liquid");

                    auto scaledEpsilon = machineEpsilon*originalki;
                    auto ki = originalki+scaledEpsilon;
                    mixture_[i].equilibriumCoefficient(ki);

                    auto xi = mixture_[i].composition()/(1+vaporMolarFraction*(ki - 1));
                    auto yi = ki*mixture_[i].composition()/(1+vaporMolarFraction*(ki - 1));
                    mixture_[i].composition("vapor", yi);
                    mixture_[i].composition("liquid", xi);

                    vaporEoS(mixture_,currentPressure_, currentTemperature_, "vapor");
                    vaporEoS.fugacities(mixture_,currentPressure_, currentTemperature_, "vapor");
                    liquidEoS(mixture_,currentPressure_, currentTemperature_, "liquid");
                    liquidEoS.fugacities(mixture_,currentPressure_, currentTemperature_, "liquid");

                    auto residual = std::log(mixture_[i].equilibriumCoefficient()) + std::log(mixture_[i].fugacityCoefficient("vapor")) -
                                            std::log(mixture_[i].fugacityCoefficient("liquid"));
                    
                    nonZeros_.push_back({i,i,(residual - residualVector_[i])/scaledEpsilon});

                    nonZeros_.push_back({size_ - 1, i, (std::accumulate(mixture_.begin(), mixture_.end(), static_cast<NP_t>(0), 
                                            [F = vaporMolarFraction](auto previous, auto const& mixComp){
                                                auto ki = mixComp.equilibriumCoefficient();
                                                auto zi = mixComp.composition();
                                                return previous + (ki-1)*zi/(1+F*(ki-1));
                                            }) - residualVector_[size_-1])/scaledEpsilon});

                    mixture_[i].equilibriumCoefficient(originalki);
                    mixture_[i].composition("vapor", originalyi);
                    mixture_[i].composition("liquid", originalxi);
                }
                //************************************************ End of Component Derivatives ************************************************


                //************************************************ Temperature Derivatives *****************************************************
                auto originalTemperature = currentTemperature_;
                auto scaledEpsilon = machineEpsilon*originalTemperature;

                currentTemperature_ = originalTemperature + scaledEpsilon;

                vaporEoS(mixture_,currentPressure_, currentTemperature_, "vapor");
                vaporEoS.fugacities(mixture_,currentPressure_, currentTemperature_, "vapor");
                liquidEoS(mixture_,currentPressure_, currentTemperature_, "liquid");
                liquidEoS.fugacities(mixture_,currentPressure_, currentTemperature_, "liquid");

                for (std::size_t i = 0; i < mixture_.size(); ++i){

                    auto residual = std::log(mixture_[i].equilibriumCoefficient()) + std::log(mixture_[i].fugacityCoefficient("vapor")) -
                                            std::log(mixture_[i].fugacityCoefficient("liquid"));

                    nonZeros_.push_back({i, size_ - 1, (residual - residualVector_[i])/scaledEpsilon});
                }

                nonZeros_.push_back({size_ - 1, size_ - 1, (std::accumulate(mixture_.begin(), mixture_.end(), static_cast<NP_t>(0), 
                                            [F = vaporMolarFraction](auto previous, auto const& mixComp){
                                                auto ki = mixComp.equilibriumCoefficient();
                                                auto zi = mixComp.composition();
                                                return previous + (ki-1)*zi/(1+F*(ki-1));
                                            }) - residualVector_[size_-1])/scaledEpsilon});

                currentTemperature_ = originalTemperature;

                vaporEoS(mixture_,currentPressure_, currentTemperature_, "vapor");
                vaporEoS.fugacities(mixture_,currentPressure_, currentTemperature_, "vapor");
                liquidEoS(mixture_,currentPressure_, currentTemperature_, "liquid");
                liquidEoS.fugacities(mixture_,currentPressure_, currentTemperature_, "liquid");
                //************************************************ End of Temperature Derivatives *********************************************

                jacobian_.setFromTriplets(nonZeros_.begin(), nonZeros_.end());
                //Initialize solver with matrix already defined
                solver_.preconditioner().setDroptol(0.00000001);
                solver_.compute(jacobian_);

                solutionDelta_ = solver_.solve(-residualVector_);

                for (std::size_t i=0; i<mixture_.size();++i){
                    auto ki = mixture_[i].equilibriumCoefficient() + solutionDelta_[i];
                    mixture_[i].equilibriumCoefficient(ki);
                }
                currentTemperature_ += solutionDelta_[size_-1];
            }

            return {0,0};
        }

        EnvelopePoints bruteForce(NP_t const maxPressure, std::size_t const& numberOfPoints){

            using VaporLiquidEquilibrium::phaseStability;
            using VaporLiquidEquilibrium::PhaseStabilityResult;

            NP_t standardConditionsPressure = 14.7; //psia
            NP_t standardConditionsTemperature = 77 + 460; //Rankine
            
            currentPressure_ = standardConditionsPressure;

            auto deltaP = (maxPressure - currentPressure_)/numberOfPoints;

            EnvelopePoints envelope;

            auto objectiveBubblePoint = [&mixture = mixture_, vaporMolarFraction=0, pressure=currentPressure_](NP_t const& temperature){
                NP_t sum = 0;
                mixture.initializeEquilibriumCoefficients(pressure, temperature); //Wilson Coeff
                for (const auto& mixComp : mixture){
                    const auto Ki = mixComp.equilibriumCoefficient();
                    sum += mixComp.composition()*(Ki-1.0)/(1.0+vaporMolarFraction*(Ki-1.0));
                }
                return sum;
            };

            auto objectiveDewPoint = [&mixture = mixture_, vaporMolarFraction=1, pressure=currentPressure_](NP_t const& temperature){
                NP_t sum = 0;
                mixture.initializeEquilibriumCoefficients(pressure, temperature); //Wilson Coeff
                for (const auto& mixComp : mixture){
                    const auto Ki = mixComp.equilibriumCoefficient();
                    sum += mixComp.composition()*(Ki-1.0)/(1.0+vaporMolarFraction*(Ki-1.0));
                }
                return sum;
            };

            auto bubbleConvergence = Math::NewtonRaphson<NP_t>(standardConditionsTemperature, objectiveBubblePoint);
            auto dewConvergence = Math::NewtonRaphson<NP_t>(standardConditionsTemperature, objectiveDewPoint);
            NP_t dewTemperature = 0;
            NP_t bubbleTemperature = 0;

            if (bubbleConvergence){
                bubbleTemperature = bubbleConvergence.value();
            }
            if (dewConvergence){
                dewTemperature = dewConvergence.value();
            }            

            while (std::abs(bubbleTemperature - dewTemperature) > 1){
                auto average = (bubbleTemperature + dewTemperature)/2;
                bubbleTemperature = average;
                dewTemperature = average;
                while(phaseStability<EoS>(mixture_, currentPressure_, bubbleTemperature) == PhaseStabilityResult::Unstable) bubbleTemperature -= 1;
                while(phaseStability<EoS>(mixture_, currentPressure_, dewTemperature) == PhaseStabilityResult::Unstable) dewTemperature += 1;
                envelope.push_back({currentPressure_, bubbleTemperature, dewTemperature});
                currentPressure_ += deltaP;
            }
            return envelope;
        }
    };

}

#endif /* PHASE_ENVELOPE_HPP */