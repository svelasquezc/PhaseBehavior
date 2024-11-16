#ifndef GENERALIZED_CUBIC_EOS_HPP
#define GENERALIZED_CUBIC_EOS_HPP

#include <vector>
#include <cmath>
#include <limits>

#include "Utilities/Math.hpp"
#include "Utilities/Types.hpp"
#include "Mixture.hpp"

using NP_t = PhaseBehavior::Types::NumericalPrecision;

namespace PhaseBehavior::EoS {

    template<
        Types::ConstexprConstant EoSParameter1,
        Types::ConstexprConstant EoSParameter2,
        typename MixingRule
        >
    class GeneralizedCubicEoS {
        private:
        constexpr static NP_t m1 = EoSParameter1();
        constexpr static NP_t m2 = EoSParameter2();
        MixingRule mixingRule;
        NP_t selectedCompressibility_ = 0;
        NP_t mixtureAttraction_, mixtureCovolume_;
        public:
        
        NP_t operator() (Component const& component, NP_t const& pressure, NP_t const& temperature){
            auto [componentAttraction, componentCovolume] = mixingRule(component, pressure, temperature);

            NP_t cuadraticTerm = (m1 + m2 - 1)*componentCovolume - 1;
            NP_t linearTerm = componentAttraction + m1*m2*std::pow(componentCovolume,2) - (m1 + m2)*componentCovolume*(componentCovolume + 1);
            NP_t constantTerm = -1*(componentAttraction*componentCovolume + m1*m2*std::pow(componentCovolume,2)*(componentCovolume + 1));

            auto roots = Math::cubicRootsFind(cuadraticTerm, linearTerm, constantTerm);

            if (roots.size()==1){
                return {roots[0]};
            }

            auto [min, max] = std::minmax_element(roots.begin(), roots.end());
            auto zMin = *min, zMax = *max;

            // Root selection process
            if(zMin < componentCovolume) {
                return zMax;
            }else{
                NP_t gibbsEnergyDerivativeWrtT = (zMax - zMin) + std::log((zMin - componentCovolume)/(zMax - componentCovolume)) 
                                            - (componentAttraction/(componentCovolume*(m2-m1)))*std::log(
                                                ((zMin + m1*componentCovolume)*(zMax + m2*componentCovolume))/((zMin + m2*componentCovolume)*(zMax + m1*componentCovolume))
                                                );
                return gibbsEnergyDerivativeWrtT > 0 ? zMin : zMax;
            }
        }

        std::vector<NP_t> operator() (Mixture const& mixture, NP_t const& pressure, NP_t const& temperature, std::string phaseName = "global"){
            auto [mixtureAttraction, mixtureCovolume] = mixingRule(mixture, pressure, temperature, phaseName);
            mixtureAttraction_ = mixtureAttraction;
            mixtureCovolume_ = mixtureCovolume;

            NP_t cuadraticTerm = (m1 + m2 - 1)*mixtureCovolume - 1;
            NP_t linearTerm = mixtureAttraction + m1*m2*std::pow(mixtureCovolume,2) - (m1 + m2)*mixtureCovolume*(mixtureCovolume + 1);
            NP_t constantTerm = -1*(mixtureAttraction*mixtureCovolume + m1*m2*std::pow(mixtureCovolume,2)*(mixtureCovolume + 1));

            auto roots = Math::cubicRootsFind(cuadraticTerm, linearTerm, constantTerm);

            if (roots.size()==1){
                selectedCompressibility_ = roots[0];
                return {roots[0]};
            }

            auto [min, max] = std::minmax_element(roots.begin(), roots.end());
            auto zMin = *min, zMax = *max;

            // Root selection process
            if(zMin < mixtureCovolume) {
                selectedCompressibility_ = zMax;
            }else{
                NP_t gibbsEnergyDerivativeWrtT = (zMax - zMin) + std::log((zMin - mixtureCovolume)/(zMax - mixtureCovolume)) 
                                            - (mixtureAttraction/(mixtureCovolume*(m2-m1)))*std::log(
                                                ((zMin + m1*mixtureCovolume)*(zMax + m2*mixtureCovolume))/((zMin + m2*mixtureCovolume)*(zMax + m1*mixtureCovolume))
                                                );
                selectedCompressibility_ = gibbsEnergyDerivativeWrtT > 0 ? zMin : zMax;
            }
            return {zMax, zMin};
        }

        NP_t selectedCompressibility() const {
            return selectedCompressibility_;
        }

        NP_t fugacity(Component const& component, NP_t const& compressibilityFactor, NP_t const& pressure, NP_t const& temperature){
            auto const Z = compressibilityFactor;
            auto [componentAttraction, componentCovolume] = mixingRule(component, pressure, temperature);
            auto A = componentAttraction;
            auto B = componentCovolume;

            return pressure*std::exp(Z - 1 - std::log(Z - B) - (A/((m1-m2)*B))*std::log((Z + m1*B)/(Z + m2*B)));
        }

        void fugacities(Mixture& mixture, std::string const& phaseName){
            const auto Z = selectedCompressibility_;
            const auto A = mixtureAttraction_;
            const auto B = mixtureCovolume_;
            for (auto& mixtureComponent : mixture){
                const auto B_i = mixingRule.componentCovolume(mixtureComponent.pure());
                NP_t sumOfBinaryAttraction = 0.0;
                for (const auto& relatedComponent : mixture){
                    sumOfBinaryAttraction += 
                    mixingRule.binaryAttraction(mixtureComponent.pure(), relatedComponent.pure())*relatedComponent.composition(phaseName); 
                }
                auto naturalLogOfFugacityCoeff = 
                    -std::log(Z - B) + 
                        (A/((m1-m2)*B))*((2.0*sumOfBinaryAttraction/A) - (B_i/B))*std::log((Z + m2*B)/(Z+ m1*B)) +
                            B_i*(Z-1.0)/B;
                mixtureComponent.fugacityCoefficient(phaseName, std::exp(naturalLogOfFugacityCoeff));
            }
        }

        NP_t volumeShift(Mixture const& mixture, std::string const& phaseName) const{
            NP_t volumeShift = 0;
            for (auto const& mixtureComponent : mixture){
                const auto B_i = mixingRule.componentShiftCovolume(mixtureComponent.pure());
                volumeShift += mixtureComponent.composition(phaseName)*mixtureComponent.pure().shift()*B_i;
            }

            return volumeShift;
        }

        NP_t mixtureCovolume() const {
            return mixtureCovolume_;
        }
    };
}

#endif /* GENERALIZED_CUBIC_EOS_HPP */