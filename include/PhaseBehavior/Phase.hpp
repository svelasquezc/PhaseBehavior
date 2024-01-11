#ifndef PHASE_HPP
#define PHASE_HPP

#include <algorithm>
#include <string>

#include "Utilities/Types.hpp"
#include "Utilities/Constants.hpp"
#include "Mixture.hpp"

using NP_t = PhaseBehavior::Types::NumericalPrecision;

namespace PhaseBehavior::Phase{
    class FluidPhase {
    protected:
        NP_t molecularWeight_;
        NP_t compressibility_;
        NP_t density_;
        NP_t viscosity_;
        NP_t molarVolume_;
    public:
        void molecularWeight(Mixture const& mixture, std::string const& phaseName){
            molecularWeight_ = std::accumulate(mixture.begin(), mixture.end(),
             static_cast<NP_t>(0.0), [](auto previous, auto element){
                return previous + element.composition(phaseName)*element.pure().molecularWeight();
             });
        }

        NP_t density(NP_t const& compressibility, NP_t const& pressure, NP_t const& temperature){
            compressibility_ = compressibility;
            density_ = pressure*(molecularWeight_/compressibility_)/(Constants::universalGasesConstant*temperature);
            return density_;
        }

        virtual NP_t viscosity(Mixture const& mixture, NP_t const& compressibility, NP_t const& pressure, NP_t const& temperature) = 0;

        NP_t molarVolume(NP_t const& compressibility, NP_t const& pressure, NP_t const& temperature, NP_t const& shift = 0){
            molarVolume_ = Constants::universalGasesConstant*temperature*compressibility/pressure - shift;
            return molarVolume_;
        }
        NP_t density() {
            density_ = molecularWeight_/molarVolume_;
            return density_;
        }
        NP_t molecularWeight(){
            return molecularWeight_;
        }
    };

    namespace Traits{

        template<typename DerivedPhase>
        class phase_traits : public FluidPhase{
        private:
            constexpr static DerivedPhase& self(){ return static_cast<DerivedPhase>(*this); }
        public:
            NP_t viscosity(Mixture const& mixture, NP_t const& compressibility, NP_t const& pressure, NP_t const& temperature) override {
                return self().viscosity(mixture, compressibility, pressure, temperature);
            }
        };
    };

    class VaporLikePhase : public Traits::phase_traits<VaporLikePhase> {
        public:
        NP_t viscosity(Mixture const& mixture, NP_t const& compressibility, NP_t const& pressure, NP_t const& temperature){
            NP_t x_v = 3.5 + 986/temperature + 0.01*molecularWeight_;
            NP_t y_v = 2.4 - 0.2*x_v;
            NP_t k_v = (9.4+0.02*molecularWeight_)*std::pow(temperature, 1.5)/(209+  temperature + 19*molecularWeight_);
            viscosity_ = 1e-4*k_v*std::exp(x_v*std::pow(density_/62.4, y_v));
            return viscosity_;
        }
    };

    class LiquidLikePhase : public Traits::phase_traits<LiquidLikePhase> {

        public:
        
        NP_t viscosity(Mixture const& mixture, NP_t const& compressibility, NP_t const& pressure, NP_t const& temperature){
            auto componentReferenceViscosity = [&temperature](Component const& component){
                auto viscosityParameter = 5.4402*std::pow(component.criticalTemperature(), 1.0/6.0)/(
                                            std::sqrt(component.molarWeight())*std::pow(component.criticalPressure(),2.0/3.0));
                auto reducedTemperature = component.reducedTemperature(temperature);
                if (reducedTemperature > 1.5){
                    return 17.78*1e-5*std::pow(4.58*reducedTemperature - 1.67, 0.625)/viscosityParameter;
                }

                return 34*1e-5*std::pow(reducedTemperature,0.94)/viscosityParameter;
            };
            auto reducedDensity = density()*mixture.pseudoCriticalVolume()/molecularWeight_;

            NP_t numerator = 0;
            NP_t denominator = 0;

            for (const auto& mixComp : mixture){
                numerator += mixComp.composition("liquid")*componentReferenceViscosity(mixComp.pure())*std::sqrt(mixComp.pure().molarWeight());
                denominator += mixComp.composition("liquid")*std::sqrt(mixComp.pure().molarWeight());
            }

            auto mixtureReferenceViscosity = numerator/denominator;

            auto mixtureViscosityParameter = 5.4402*std::pow(mixture.pseudoCriticalTemperature(),1.0/6.0)/(std::sqrt(molecularWeight_)*std::pow(mixture.pseudoCriticalPressure(),2.0/3.0));

            // LBC Correlation
            viscosity_ = mixtureReferenceViscosity + (std::pow(0.1023 + 0.023364*reducedDensity + 0.058533*std::pow(reducedDensity,2) - 
                                                    0.040758*std::pow(reducedDensity,3) + 0.0093324*std::pow(reducedDensity,4), 4) - 1e-4)/mixtureViscosityParameter;

            return viscosity_;
        }
    };
};


#endif /* PHASE_HPP */