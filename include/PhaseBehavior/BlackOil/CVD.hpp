#ifndef CVD_HPP
#define CVD_HPP

#include <tuple>

#include "../Utilities/Types.hpp"
#include "../Utilities/Constants.hpp"
#include "../Mixture.hpp"
#include "../PhaseEquilibrium.hpp"
#include "../Phase.hpp"

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

        struct StockTankOilProperties {
            NP_t gasFraction;
            NP_t oilFraction;
            NP_t molecularWeight;
            NP_t density;
        };

        std::vector<BlackOilPVTRow> PVTTable;

        public:
        BlackOilPropertiesEstimation(Mixture mixture, 
        NP_t const& reservoirPressure, NP_t const& reservoirTemperature, 
        NP_t const& stage1Pressure, NP_t const& stage1Temperature,
        NP_t const& stage2Pressure, NP_t const& stage2Temperature,
        NP_t const& stockTankPressure, NP_t const& stockTankTemperature,
        NP_t const& abandonmentPressure){

            using VaporLiquidEquilibrium::phaseStability;
            using VaporLiquidEquilibrium::PhaseStabilityResult;

            EoS eos;

            NP_t dewPressure = reservoirPressure;
            while(phaseStability<EoS>(mixture, dewPressure, reservoirTemperature) == PhaseStabilityResult::Stable) dewPressure-=1; //psia 

            NP_t initialGasEquivalentVolume = 1e6;
            auto totalMoles = Constants::standardConditionsPressure * initialGasEquivalentVolume/
                                (Constants::universalGasesConstant * Constants::standardConditionsTemperature);
            
            eos(mixture, dewPressure, reservoirTemperature);
            mixture.compressibility("global", eos.selectedCompressibility());
            auto gasPhase = Phase::singlePhaseIdentification(mixture, mixture.compressibility("global"),dewPressure, reservoirTemperature, eos);

            auto pvtMolecularWeight = gasPhase->molecularWeight();
            auto pvtDensity = gasPhase->density();

            auto associatedPVTVolume = totalMoles * pvtMolecularWeight/pvtDensity;

            auto [initialGasFraction, initialOilFraction, stockTankMW, stockTankDensity] = oilSeparator(mixture, stage1Pressure, stage1Temperature,
                                                                                            stage2Pressure, stage2Temperature, stockTankPressure, stockTankTemperature, eos);

            auto totalSurfaceGas = initialGasFraction * totalMoles * 379.56;
            auto totalStockTankOil = initialOilFraction * totalMoles * stockTankMW/(stockTankDensity*5.615);

            auto initialGasFormationVolumeFactor = associatedPVTVolume/(totalSurfaceGas*5.615);
            auto initialVolatilizedOil = 1e6*totalStockTankOil/totalSurfaceGas;

            auto pressure = dewPressure;
            while(pressure > abandonmentPressure){
                
            }
        }

        StockTankOilProperties oilSeparator(Mixture const& mixture, NP_t const& stage1Pressure, NP_t const& stage1Temperature,
                                            NP_t const& stage2Pressure, NP_t const& stage2Temperature,
                                            NP_t const& stockTankPressure, NP_t const& stockTankTemperature, EoS& eos){

            auto firstStageMixture = mixture;
            auto secondStageMixture = mixture;
            auto stockTankMixture = mixture;

            VaporLiquidEquilibrium::isothermalTwoPhaseFlash(firstStageMixture, stage1Pressure, stage1Temperature);
            
            for (std::size_t i = 0; i < mixture.size(); ++i){
                secondStageMixture[i].composition("global", firstStageMixture[i].composition("liquid"));
            }

            VaporLiquidEquilibrium::isothermalTwoPhaseFlash(secondStageMixture, stage2Pressure, stage2Temperature);

            for (std::size_t i = 0; i < mixture.size(); ++i){
                stockTankMixture[i].composition("global", secondStageMixture[i].composition("liquid"));
            }

            VaporLiquidEquilibrium::isothermalTwoPhaseFlash(stockTankMixture, stockTankPressure, stockTankTemperature);

            auto oilFraction = firstStageMixture.molarFraction("liquid")*secondStageMixture.molarFraction("liquid")*stockTankMixture.molarFraction("liquid");
            auto gasFraction = 1 - oilFraction;

            auto oilPhase = Phase::singlePhaseIdentification(stockTankMixture, mixture.compressibility("liquid"),stockTankPressure, stockTankTemperature, eos);

            auto molecularWeight = oilPhase->molecularWeight();
            auto density = oilPhase->density();

            return {gasFraction, oilFraction, molecularWeight, density};
        }
    };
}

#endif /* CVD_HPP */