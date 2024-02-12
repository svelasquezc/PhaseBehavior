#ifndef CVD_HPP
#define CVD_HPP

#include "../Utilities/Types.hpp"
#include "../Utilities/Constants.hpp"
#include "../Mixture.hpp"
#include "../PhaseEquilibrium.hpp"

using NP_t = Types::NumericalPrecision;

namespace PhaseBehavior::BlackOil {

    template<typename EoS>
    class BlackOilPropertiesEstimation {
        private:
        struct BlackOilPVTRow {
            NP_t pressure;
            NP_t oilVolumetricFactor;
            NP_t gasVolumetricFactor;
            NP_t dissolvedGas;
            NP_t volatilizedOil;
        };
        std::vector<BlackOilPVTRow> PVTTable;

        public:
        BlackOilPropertiesEstimation(Mixture& mixture, 
        NP_t reservoirPressure, NP_t reservoirTemperature, 
        NP_t stage1Pressure, NP_t stage1Temperature,
        NP_t stage2Pressure, NP_t stage2Temperature,
        NP_t stockTankPressure, NP_t stockTankTemperature){

            using VaporLiquidEquilibrium::phaseStability;
            using VaporLiquidEquilibrium::PhaseStabilityResult;

            NP_t dewPressure = reservoirPressure;
            while(phaseStability<EoS>(mixture, dewPressure, bubbleTemperature) == PhaseStabilityResult::Stable) dewPressure-=1; //psia 

            NP_t initialGasEquivalentVolume = 1e6;
            auto totalMoles = Constants::standardConditionsPressure * initialGasEquivalentVolume/
                                (Constants::universalGasesConstant * Constants::standardConditionsTemperature);
            
            


        }
    };
}

#endif /* CVD_HPP */