#ifndef GENERALIZED_CUBIC_EOS_HPP
#define GENERALIZED_CUBIC_EOS_HPP

#include <vector>
#include <cmath>
#include <limits>

#include "Utilities/Math.hpp"
#include "Utilities/Types.hpp"
#include "Mixture.hpp"


namespace PhaseBehavior {

    using NP_t = Types::NumericalPrecision;

    namespace EoS {

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
            public:
            
            std::vector<NP_t> operator() (Mixture const& mixture, NP_t const& pressure, NP_t const& temperature, std::string phaseName = "global"){
                auto [mixtureAttraction, mixtureCovolume] = mixingRule(mixture, pressure, temperature, phaseName);

                NP_t cuadraticTerm = (m1 + m2 - 1)*mixtureCovolume - 1;
                NP_t linearTerm = mixtureAttraction + m1*m2*std::pow(mixtureCovolume,2) - (m1 + m2)*mixtureCovolume*(mixtureCovolume + 1);
                NP_t constantTerm = -1*(mixtureAttraction*mixtureCovolume + m1*m2*std::pow(mixtureCovolume,2)*(mixtureCovolume + 1));

                auto roots = Math::cubicRootsFind(cuadraticTerm, linearTerm, constantTerm);

                if (roots.size()==1){
                    selectedCompressibility_ = roots[0];
                    return {roots[0]};
                }

                NP_t zMin = std::numeric_limits<NP_t>::max();
                NP_t zMax = std::numeric_limits<NP_t>::min();
                
                for (const auto& root : roots){
                    if (root>zMax) zMax = root;
                    if (root<zMin) zMin = root; 
                }

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
        };
    }
}

#endif /* GENERALIZED_CUBIC_EOS_HPP */