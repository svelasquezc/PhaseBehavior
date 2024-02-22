#ifndef CVD_HPP
#define CVD_HPP

#include <tuple>

#include "../Utilities/Types.hpp"
#include "../Utilities/Constants.hpp"
#include "../Mixture.hpp"
#include "../PhaseEquilibrium.hpp"
#include "../Phase.hpp"

using NP_t = PhaseBehavior::Types::NumericalPrecision;

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

        std::vector<BlackOilPVTRow> PVTTable_;

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
            auto fluid = Phase::singlePhaseIdentification(mixture, mixture.compressibility("global"), dewPressure, reservoirTemperature, eos);

            auto pvtMolecularWeight = fluid->molecularWeight();
            auto pvtDensity = fluid->density();

            auto associatedPVTVolume = totalMoles * pvtMolecularWeight/pvtDensity;

            auto [surfaceGasMolarFraction, stockTankOilMolarFraction, stockTankOilMW, stockTankOilDensity] = separator(mixture, stage1Pressure, stage1Temperature,
                                                                                            stage2Pressure, stage2Temperature, stockTankPressure, stockTankTemperature, eos);

            auto totalSurfaceGas = surfaceGasMolarFraction * totalMoles * 379.56;
            auto totalStockTankOil = stockTankOilMolarFraction * totalMoles * stockTankOilMW/(stockTankOilDensity*5.615);

            auto initialGasVolumetricFactor = associatedPVTVolume/(totalSurfaceGas*5.615);
            auto initialVolatilizedOil = 1e6*totalStockTankOil/totalSurfaceGas;

            PVTTable_.push_back({dewPressure, 0.0, initialGasVolumetricFactor, 0.0, initialVolatilizedOil});

            NP_t reservoirOilMoles  = 0;
            NP_t gasAndExcessGasMoles = totalMoles;
            NP_t reservoirOilVolume = 0;
            NP_t gasAndExcessGasVolume = associatedPVTVolume;
            NP_t reservoirGasVolume = 0;

            auto pressure = dewPressure - 1;


            while(pressure >= abandonmentPressure){
                VaporLiquidEquilibrium::isothermalTwoPhaseFlash<EoS>(mixture, pressure, reservoirTemperature);
                auto gasPhase = Phase::VaporLikePhase(mixture, mixture.compressibility("vapor"), pressure, reservoirTemperature, eos);
                auto oilPhase = Phase::LiquidLikePhase(mixture, mixture.compressibility("liquid"), pressure, reservoirTemperature, eos);
                reservoirOilMoles = totalMoles*(1-mixture.molarFraction("vapor"));
                gasAndExcessGasMoles = totalMoles*mixture.molarFraction("vapor");

                gasAndExcessGasVolume = gasAndExcessGasMoles*gasPhase.molecularWeight()/(gasPhase.density()*5.615);
                reservoirOilVolume = reservoirOilMoles*oilPhase.molecularWeight()/(oilPhase.density()*5.615);
                reservoirGasVolume = associatedPVTVolume - reservoirOilVolume;
                auto excessGasVolume = gasAndExcessGasVolume - reservoirGasVolume;
                auto reservoirGasMoles = reservoirGasVolume*gasPhase.density()*5.615/gasPhase.molecularWeight();
                auto excessGasMoles = excessGasVolume*gasPhase.density()*5.615/gasPhase.molecularWeight();

                auto [surfaceGasInGasMolarFraction, stockTankOilInGasMolarFraction,
                        stockTankOilInGasMolWeight, stockTankOilInGasDensity] = separator(mixture, stage1Pressure, stage1Temperature,
                                                                            stage2Pressure, stage2Temperature, stockTankPressure,
                                                                                stockTankTemperature, eos, "vapor");

                auto standardConditionsGasInGasVolume = surfaceGasInGasMolarFraction*reservoirGasMoles*379.56;
                auto stockTankOilInGasVolume = stockTankOilInGasMolarFraction*reservoirOilMoles*stockTankOilInGasMolWeight/(stockTankOilInGasDensity*5.615);

                auto [surfaceGasInOilMolarFraction, stockTankOilInOilMolarFraction,
                        stockTankOilInOilMolWeight, stockTankOilInOilDensity] = separator(mixture, stage1Pressure, stage1Temperature,
                                                                            stage2Pressure, stage2Temperature, stockTankPressure,
                                                                                stockTankTemperature, eos, "liquid");

                auto standardConditionsGasInOilVolume = surfaceGasInOilMolarFraction*reservoirGasMoles*379.56;
                auto stockTankOilInOilVolume = stockTankOilInOilMolarFraction*reservoirOilMoles*stockTankOilInOilMolWeight/(stockTankOilInOilDensity*5.615);

                auto gasVolumetricFactor = reservoirGasVolume/standardConditionsGasInGasVolume;
                auto oilVolumetricFactor = reservoirOilVolume/stockTankOilInOilVolume;

                auto volatilizedOil = stockTankOilInGasVolume/standardConditionsGasInGasVolume;
                auto dissolvedGas = standardConditionsGasInOilVolume/stockTankOilInOilVolume;
                if (static_cast<int>(pressure) % 50 == 0){
                    PVTTable_.push_back({pressure, oilVolumetricFactor, gasVolumetricFactor, dissolvedGas, volatilizedOil});
                }
                for (auto& mixtureComponent : mixture){
                    mixtureComponent.composition((mixtureComponent.composition()*totalMoles - mixtureComponent.composition("vapor")*excessGasMoles)/totalMoles);
                }
                totalMoles = totalMoles - excessGasMoles;
                pressure -= 1;
            }
        }

        StockTankOilProperties separator(Mixture const& mixture, NP_t const& stage1Pressure, NP_t const& stage1Temperature,
                                            NP_t const& stage2Pressure, NP_t const& stage2Temperature,
                                            NP_t const& stockTankPressure, NP_t const& stockTankTemperature, EoS& eos, std::string separationType = "global"){

            auto firstStageMixture = mixture;
            auto secondStageMixture = mixture;
            auto stockTankMixture = mixture;

            for (std::size_t i = 0; i < mixture.size(); ++i){
                firstStageMixture[i].composition("global", mixture[i].composition(separationType));
            }

            VaporLiquidEquilibrium::isothermalTwoPhaseFlash<EoS>(firstStageMixture, stage1Pressure, stage1Temperature);
            
            for (std::size_t i = 0; i < mixture.size(); ++i){
                secondStageMixture[i].composition("global", firstStageMixture[i].composition("liquid"));
            }

            VaporLiquidEquilibrium::isothermalTwoPhaseFlash<EoS>(secondStageMixture, stage2Pressure, stage2Temperature);

            for (std::size_t i = 0; i < mixture.size(); ++i){
                stockTankMixture[i].composition("global", secondStageMixture[i].composition("liquid"));
            }

            VaporLiquidEquilibrium::isothermalTwoPhaseFlash<EoS>(stockTankMixture, stockTankPressure, stockTankTemperature);

            auto oilFraction = firstStageMixture.molarFraction("liquid")*secondStageMixture.molarFraction("liquid")*stockTankMixture.molarFraction("liquid");
            auto gasFraction = 1 - oilFraction;

            auto oilPhase = Phase::singlePhaseIdentification(stockTankMixture, mixture.compressibility("liquid"),stockTankPressure, stockTankTemperature, eos);

            auto molecularWeight = oilPhase->molecularWeight();
            auto density = oilPhase->density();

            return {gasFraction, oilFraction, molecularWeight, density};
        }

        auto PVTTable(){
            return PVTTable_;
        }
    };
}

#endif /* CVD_HPP */