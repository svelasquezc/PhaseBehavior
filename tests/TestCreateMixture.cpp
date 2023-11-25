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
#include <PhaseBehavior/Utilities/Input.hpp>
#include <PhaseBehavior/Component.hpp>

template<typename DeducedType>
class TypeChecker;

using Precision_t = PhaseBehavior::Types::NumericalPrecision;
using MixComp = std::pair<PhaseBehavior::Component, Precision_t>;

TEST_CASE("Can create Mixture objects", "[mixture]"){

    SECTION("Mixture object from file"){
        
        auto mixture = PhaseBehavior::Input::createMixtureFromFile("PVT.csv", "InteractionCoefficients.csv");

        CHECK(Catch::Approx(mixture.interactionCoefficient("CO2", "C1"))==0.105);
        CHECK(Catch::Approx(mixture.interactionCoefficient("CO2", "C2"))==0.13);
        CHECK(Catch::Approx(mixture.interactionCoefficient("CO2", "C3"))==0.125);
        CHECK(Catch::Approx(mixture.interactionCoefficient("CO2", "i-C4"))==0.12);
        CHECK(Catch::Approx(mixture.interactionCoefficient("CO2", "n-C4"))==0.115);
        CHECK(Catch::Approx(mixture.interactionCoefficient("CO2", "C7+"))==0.115);
        CHECK(Catch::Approx(mixture.interactionCoefficient("C1", "C2"))==0.0);
        CHECK(Catch::Approx(mixture.interactionCoefficient("n-C6", "i-C4"))==0.0);

    }

    SECTION("Mixture object from known components"){
        PhaseBehavior::Component CO2 {"CO2",1071,547.91,0.2667,4401,0.0344};
        PhaseBehavior::Component C1 {"C1",1071,547.91,0.2667,4401,0.0344};
        PhaseBehavior::Component C2 {"C2",1071,547.91,0.2667,4401,0.0344};
        MixComp CO2inMix = {CO2, 0.0031}; 
        MixComp C1inMix = {C1, 0.6192};
        MixComp C2inMix = {C2, 0.3192};
        
        //REQUIRE_THROWS(Mixture(CO2inMix));
        auto m1 = PhaseBehavior::Mixture(CO2inMix, C1inMix);
        auto m2 = PhaseBehavior::Mixture(CO2inMix, C1inMix, C2inMix);
    }

}