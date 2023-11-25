#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <PhaseBehavior/Utilities/Input.hpp>

#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/MixingRules.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>


using PRMixingRule = PhaseBehavior::EoS::MixingRules::NonRandomMixingRule<
PhaseBehavior::EoS::PR::constants::omegaA, 
PhaseBehavior::EoS::PR::constants::omegaB, 
PhaseBehavior::EoS::PR::PengRobinson
>;

TEST_CASE("Can evaluate the attraction and covolume of a Mixture", "[NonRandomMixingRule]"){
    auto mixture = PhaseBehavior::Input::createMixtureFromFile("PVT.csv", "InteractionCoefficients.csv");

    PRMixingRule mixingRule;

    auto [A1,B1] = mixingRule(mixture, 500 /*psia*/, 50 + 459.67 /*R*/);
    auto [A2,B2] = mixingRule(mixture, 400 /*psia*/, 50 + 459.67 /*R*/);
    auto [A3,B3] = mixingRule(mixture, 300 /*psia*/, 50 + 459.67 /*R*/);
    auto [A4,B4] = mixingRule(mixture, 500 /*psia*/, 100 + 459.67 /*R*/);
    auto [A5,B5] = mixingRule(mixture, 400 /*psia*/, 100 + 459.67 /*R*/);
    auto [A6,B6] = mixingRule(mixture, 300 /*psia*/, 100 + 459.67 /*R*/);

    CHECK(Catch::Approx(A1)==0.4401); CHECK(Catch::Approx(B1)==0.0693);
    CHECK(Catch::Approx(A2)==0.3521); CHECK(Catch::Approx(B2)==0.0554);
    CHECK(Catch::Approx(A3)==0.2641); CHECK(Catch::Approx(B3)==0.0416);
    CHECK(Catch::Approx(A4)==0.3437); CHECK(Catch::Approx(B4)==0.0631);
    CHECK(Catch::Approx(A5)==0.2778); CHECK(Catch::Approx(B5)==0.0505);
    CHECK(Catch::Approx(A6)==0.2084); CHECK(Catch::Approx(B6)==0.0379);
    
}