#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <PhaseBehavior/Utilities/Input.hpp>
#include <PhaseBehavior/Utilities/Math.hpp>

#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/MixingRules.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>


using PRMixingRule = PhaseBehavior::EoS::MixingRules::RandomMixingRule<
PhaseBehavior::EoS::PR::constants::omegaA, 
PhaseBehavior::EoS::PR::constants::omegaB, 
PhaseBehavior::EoS::PR::PengRobinson
>;

TEST_CASE("Can evaluate the attraction and covolume of a Mixture", "[NonRandomMixingRule]"){
    auto mixture = PhaseBehavior::Input::createMixtureFromFile("PVT.csv", "InteractionCoefficients.csv");

    PRMixingRule mixingRule;

    auto [A1,B1] = mixingRule(mixture, 500 /*psia*/, 50 + 460 /*R*/);
    auto [A2,B2] = mixingRule(mixture, 400 /*psia*/, 50 + 460 /*R*/);
    auto [A3,B3] = mixingRule(mixture, 300 /*psia*/, 50 + 460 /*R*/);
    auto [A4,B4] = mixingRule(mixture, 500 /*psia*/, 100 + 460 /*R*/);
    auto [A5,B5] = mixingRule(mixture, 400 /*psia*/, 100 + 460 /*R*/);
    auto [A6,B6] = mixingRule(mixture, 300 /*psia*/, 100 + 460 /*R*/);

    CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(A1,4))==0.4401); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(B1,4))==0.0693);
    CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(A2,4))==0.3521); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(B2,4))==0.0554);
    CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(A3,4))==0.2641); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(B3,4))==0.0416);
    CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(A4,4))==0.3473); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(B4,4))==0.0631);
    CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(A5,4))==0.2778); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(B5,4))==0.0505);
    CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(A6,4))==0.2084); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(B6,4))==0.0379);
    
}