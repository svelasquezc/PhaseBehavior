#include <iostream>
#include <tuple>
#include <fstream>
#include <sstream>
#include <array>

#include <PhaseBehavior/Component.hpp>
#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/Phase.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>
#include <PhaseBehavior/PhaseEquilibrium.hpp>

using PR = PhaseBehavior::EoS::PR::PengRobinson;
using SRK = PhaseBehavior::EoS::SRK::SoaveRedlichKwong;

struct Composition {
    double x1; // e.g., C1
    double x2; // e.g., C3
    double x3; // e.g., C8
};

// Computes the total number of valid points on the ternary grid
constexpr std::size_t get_num_points(std::size_t N) {
    return (N + 1) * (N + 2) / 2;
}

// Generates the grid at compile time
template <std::size_t N>
constexpr auto generate_ternary_grid() {
    constexpr std::size_t num_points = get_num_points(N);
    std::array<Composition, num_points> grid{};
    
    std::size_t index = 0;
    double step = 1.0 / static_cast<double>(N);
    
    for (std::size_t i = 0; i <= N; ++i) {
        double x1 = i * step;
        
        for (std::size_t j = 0; j <= N - i; ++j) {
            double x2 = j * step;
            double x3 = 1.0 - x1 - x2;
            
            // Mitigate minor floating-point inaccuracies
            if (x3 < 1e-15) x3 = 0.0; 
            if (x3 > 1.0) x3 = 1.0;

            grid[index++] = {x1, x2, x3};
        }
    }
    return grid;
}


int main(){

    std::ofstream file("C1C3C8Densities.csv");
    
    using PhaseBehavior::Component, PhaseBehavior::VaporLiquidEquilibrium::isothermalTwoPhaseFlash;

    auto C1 = Component {"CH4",      4599.2,    190.56,   6.1477929425e-3, 16.04 ,  0.0114};
    auto C3 = Component {"C3H8",     4251.2,    369.89,   4.5355587810e-3, 44.097,  0.1521};
    auto C8 = Component {"n-C8H18",  2497.0,    569.32,   4.2571306948e-3, 114.23,  0.3950};

    //auto casePressure = 9.5e3 /*[kPa]*/;
    //auto casePressure = 5.3e3 /*[kPa]*/;
    auto casePressure = 120e2 /*[kPa]*/;
    auto caseTemperature = 400.0 /*[K]*/;

    auto grid = generate_ternary_grid<200>();

    file << "xC1,xC3,xC8,Density" <<std::endl;

    for(auto [zC1, zC3, zC8] : grid){

        file << zC1 << "," << zC3 << "," << zC8 << ",";

        PhaseBehavior::Mixture C1C3C8{{C1, zC1}, {C3, zC3}, {C8, zC8}};

        auto result = isothermalTwoPhaseFlash<PR>(C1C3C8, casePressure, caseTemperature);

        if (result == PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable){

            auto liquidPhase = PhaseBehavior::Phase::LiquidLikePhase(C1C3C8);
            auto liquidMolarDensity = 1.0 / liquidPhase.molarVolume(C1C3C8.compressibility(PhaseBehavior::PhaseName::liquid), casePressure, caseTemperature);
            auto gasPhase = PhaseBehavior::Phase::VaporLikePhase(C1C3C8);
            auto gasMolarDensity = 1.0 / gasPhase.molarVolume(C1C3C8.compressibility(PhaseBehavior::PhaseName::vapor), casePressure, caseTemperature);

            auto mixtureGasMolarFraction = C1C3C8.molarFraction(PhaseBehavior::PhaseName::vapor);

            auto mixtureAvgDensity = (gasMolarDensity*mixtureGasMolarFraction + liquidMolarDensity*(1.0 - mixtureGasMolarFraction))*C1C3C8.molarWeight();
            file << mixtureAvgDensity << std::endl;
        }else{
            PR prEoS;
            prEoS(C1C3C8, casePressure, caseTemperature);
            C1C3C8.compressibility(PhaseBehavior::PhaseName::global, prEoS.selectedCompressibility());
            auto phase = PhaseBehavior::Phase::singlePhaseIdentification(C1C3C8, C1C3C8.compressibility(PhaseBehavior::PhaseName::global), casePressure, caseTemperature, prEoS);
            file << phase->density() << std::endl;
        }   
    }

    file.close();
    return 0;
}