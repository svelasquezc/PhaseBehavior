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
        NP_t volume_;
    public:
        void molecularWeight(Mixture const& mixture, std::string const& phaseName){
            molecularWeight_ = std::accumulate(mixture.begin(), mixture.end(),
             static_cast<NP_t>(0.0), [](auto previous, auto element){
                return previous + element.composition(phaseName)*element.pure().molecularWeight();
             });
        }

        void density(NP_t const& compressibility, NP_t const& pressure, NP_t const& temperature){
            compressibility_ = compressibility;
            density_ = pressure*(molecularWeight_/compressibility_)/(Constants::universalGasesConstant*temperature);
        }
        virtual void viscosity() = 0;
        NP_t volume(){
            return volume_;
        }
        NP_t density() {
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
            void viscosity() override {
                self().viscosity();
            }
        };
    };

    class VaporLikePhase : public Traits::phase_traits<VaporLikePhase> {

    };

    class LiquidLikePhase : public Traits::phase_traits<LiquidLikePhase> {

    };
};


#endif /* PHASE_HPP */