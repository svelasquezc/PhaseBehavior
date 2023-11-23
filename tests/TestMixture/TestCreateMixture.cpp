#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include <cmath>

#include <PhaseBehavior/Utilities/Types.hpp>
#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/Component.hpp>

template<typename DeducedType>
        class TypeChecker;

using Precision_t = TypesDefinition::NumericalPrecision;
using MixComp = std::pair<Component, Precision_t>;

TEST_CASE("Can create Mixture objects", "[mixture]"){

    SECTION("Mixture object from file"){
        
        std::vector<MixComp> components;

        {
            auto pvtFile = std::ifstream("PVT.csv");
            std::string line;

            // Get rid of the header line
            std::getline(pvtFile, line);
            while (std::getline(pvtFile, line)){
                std::istringstream ss(line);
                std::string name;
                Precision_t composition, Pc, Tc, w, MW, Vc;
                ss >> name >> composition >> Pc >> Tc >> w >> MW >> Vc;
                components.push_back({Component(std::move(name), Pc, Tc, Vc, MW, w ), composition});
            }
            pvtFile.close();
        }
        Mixture mixture(components);
        {
            auto interactionCoefficientsFile = std::ifstream("InteractionCoefficients.csv");
            std::string line;

            std::getline(interactionCoefficientsFile, line);

            for (const auto& mixComponent1 : mixture){
                std::getline(interactionCoefficientsFile, line);
                std::istringstream ss(line);
                std::string name;
                ss >> name;
                for (const auto& mixComponent2 : mixture){
                    Precision_t interactionCoefficientValue = 0.0;
                    ss >> interactionCoefficientValue;
                    if (std::abs(interactionCoefficientValue) > 0.0){
                        mixture.interactionCoefficient(mixComponent1.pure()->name(), mixComponent2.pure()->name(),interactionCoefficientValue);
                    }
                }
            }
        }

        CHECK(Catch::Approx(mixture.interactionCoefficient("CO2", "C1"))==0.105);
        CHECK(Catch::Approx(mixture.interactionCoefficient("CO2", "C2"))==0.13);
        CHECK(Catch::Approx(mixture.interactionCoefficient("CO2", "C3"))==0.125);
        CHECK(Catch::Approx(mixture.interactionCoefficient("CO2", "i-C4"))==0.12);
        CHECK(Catch::Approx(mixture.interactionCoefficient("CO2", "n-C4"))==0.115);
        CHECK(Catch::Approx(mixture.interactionCoefficient("CO2", "C7+"))==0.115);

    }

    SECTION("Mixture object from known components"){
        Component CO2 {"CO2",1071,547.91,0.2667,4401,0.0344};
        Component C1 {"C1",1071,547.91,0.2667,4401,0.0344};
        Component C2 {"C2",1071,547.91,0.2667,4401,0.0344};
        MixComp CO2inMix = {CO2, 0.0031}; 
        MixComp C1inMix = {C1, 0.6192};
        MixComp C2inMix = {C2, 0.3192};
        
        //REQUIRE_THROWS(Mixture(CO2inMix));
        auto m1 = Mixture(CO2inMix, C1inMix);
        auto m2 = Mixture(CO2inMix, C1inMix, C2inMix);
    }

}