#ifndef PHASE_HPP
#define PHASE_HPP

#include <algorithm>
#include <string>
#include <memory>

#include "Utilities/Types.hpp"
#include "Utilities/Constants.hpp"
#include "Mixture.hpp"

using NP_t = PhaseBehavior::Types::NumericalPrecision;

namespace PhaseBehavior::Phase{

    class FluidPhase {
    protected:
        std::string phaseName_;
        NP_t molecularWeight_;
        NP_t compressibility_;
        NP_t density_;
        NP_t viscosity_;
        NP_t molarVolume_;
    public:
        void molecularWeight(Mixture const& mixture, std::string const& compositionType){
            molecularWeight_ = std::accumulate(mixture.begin(), mixture.end(),
             static_cast<NP_t>(0.0), [phaseName_=this->phaseName_, &compositionType](auto previous, auto& element){
                return previous + element.composition(compositionType)*element.pure().molarWeight();
             });
        }

        NP_t density(NP_t const& compressibility, NP_t const& pressure, NP_t const& temperature){
            compressibility_ = compressibility;
            density_ = pressure*(molecularWeight_/compressibility_)/(Constants::universalGasesConstant*temperature);
            return density_;
        }

        virtual NP_t viscosity(Mixture const& mixture, NP_t const& compressibility, NP_t const& pressure, NP_t const& temperature) = 0;

        template<typename EoS>
        NP_t molarVolume(Mixture const& mixture, NP_t const& compressibility, NP_t const& pressure, NP_t const& temperature, EoS const& eos, std::string compositionType = ""){
            if(compositionType == "") compositionType = phaseName_;
            molarVolume_ = Constants::universalGasesConstant*(temperature*compressibility/pressure - eos.volumeShift(mixture, compositionType));
            return molarVolume_;
        }

        
        NP_t molarVolume(NP_t const& compressibility, NP_t const& pressure, NP_t const& temperature){
            molarVolume_ = Constants::universalGasesConstant*temperature*compressibility/pressure;
            return molarVolume_;
        }

        NP_t density(){
            density_ = molecularWeight_/molarVolume_;
            return density_;
        }
        NP_t molecularWeight() const{
            return molecularWeight_;
        }

        NP_t viscosity() const{
            return viscosity_;
        }
    };

    namespace Traits{

        template<typename DerivedPhase>
        class phase_traits : public FluidPhase{
        private:
            constexpr DerivedPhase& self(){ return static_cast<DerivedPhase&>(*this); }
        public:

            NP_t viscosity(Mixture const& mixture, NP_t const& compressibility, NP_t const& pressure, NP_t const& temperature) override {
                return self().viscosity(mixture, compressibility, pressure, temperature);
            }
        };
    };

    class VaporLikePhase : public Traits::phase_traits<VaporLikePhase> {
        public:

        VaporLikePhase(Mixture const& mixture, std::string compositionType = ""){
            phaseName_ = "vapor";
            if (compositionType == "") compositionType = phaseName_;
            molecularWeight(mixture, compositionType);
        }

        NP_t viscosity(Mixture const& mixture, NP_t const& compressibility, NP_t const& pressure, NP_t const& temperature){
            NP_t x_v = 3.5 + 986/temperature + 0.01*molecularWeight_;
            NP_t y_v = 2.4 - 0.2*x_v;
            NP_t k_v = (9.4+0.02*molecularWeight_)*std::pow(temperature, 1.5)/(209 + temperature + 19*molecularWeight_);
            viscosity_ = 1e-4*k_v*std::exp(x_v*std::pow(density_/62.4, y_v));
            return viscosity_;
        }

        NP_t viscosity(){
            return FluidPhase::viscosity();
        }
    };

    class LiquidLikePhase : public Traits::phase_traits<LiquidLikePhase> {

        public:

        LiquidLikePhase(Mixture const& mixture, std::string compositionType = ""){
            phaseName_ = "liquid";
            if (compositionType == "") compositionType = phaseName_;
            molecularWeight(mixture, compositionType);
        }
        
        NP_t viscosity(Mixture const& mixture, NP_t const& compressibility, NP_t const& pressure, NP_t const& temperature){
            auto componentReferenceViscosity = [&temperature](Component const& component){

                NP_t viscosityParameter = 5.4402*std::pow(component.criticalTemperature(), 1.0/6.0)/(
                                            std::sqrt(component.molarWeight())*std::pow(component.criticalPressure(),2.0/3.0));

                NP_t reducedTemperature = component.reducedTemperature(temperature);

                if (reducedTemperature <= 1.5){
                    return static_cast<NP_t>(34.0*1e-5*std::pow(reducedTemperature,0.94)/viscosityParameter);    
                }
                
                return static_cast<NP_t>(17.78*1e-5*std::pow(4.58*reducedTemperature - 1.67, 0.625)/viscosityParameter);
                
            };

            auto reducedDensity = density_*mixture.pseudoCriticalVolume()/molecularWeight_;

            NP_t numerator = 0;
            NP_t denominator = 0;

            for (const auto& mixComp : mixture){
                auto molarWeightPonderation = mixComp.composition(phaseName_)*std::sqrt(mixComp.pure().molarWeight());
                numerator += componentReferenceViscosity(mixComp.pure())*molarWeightPonderation;
                denominator += molarWeightPonderation;
            }
            
            auto mixtureReferenceViscosity = numerator/denominator;

            auto mixtureViscosityParameter = 5.4402*std::pow(mixture.pseudoCriticalTemperature(),1.0/6.0)/(std::sqrt(molecularWeight())*std::pow(mixture.pseudoCriticalPressure(),2.0/3.0));

            // LBC Correlation
            viscosity_ = mixtureReferenceViscosity + (std::pow(0.1023 + 0.023364*reducedDensity + 0.058533*std::pow(reducedDensity,2.0) - 
                                                    0.040758*std::pow(reducedDensity,3.0) + 0.0093324*std::pow(reducedDensity,4.0), 4.0) - 1e-4)/mixtureViscosityParameter;

            return viscosity_;
        }

        NP_t viscosity(){
            return FluidPhase::viscosity();
        }
    };

    template<typename EoS>
    std::unique_ptr<FluidPhase> singlePhaseIdentification(Mixture const& mixture, NP_t const& compressibility, NP_t const& pressure, NP_t const& temperature, EoS const& eos){
        std::unique_ptr<FluidPhase> liquidLike = std::make_unique<LiquidLikePhase>(mixture, "global");
        std::unique_ptr<FluidPhase> vaporLike = std::make_unique<VaporLikePhase>(mixture, "global");

        auto molarVolume = liquidLike->molarVolume(mixture, compressibility, pressure, temperature, eos, "global");
        vaporLike->molarVolume(mixture, compressibility, pressure, temperature, eos, "global");
        auto mixtureCovolume = eos.mixtureCovolume();
        if(molarVolume/mixtureCovolume < 1.75) return liquidLike;
        return vaporLike;
    };
};


#endif /* PHASE_HPP */