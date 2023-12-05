#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <PhaseBehavior/Utilities/Input.hpp>
#include <PhaseBehavior/Utilities/Math.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>

TEST_CASE("Can create and evaluate an EoS", "[PengRobinson]"){
    auto mixture = PhaseBehavior::Input::createMixtureFromFile("PVT.csv", "InteractionCoefficients.csv");

    PhaseBehavior::EoS::SRK::SoaveRedlichKwong srk;

    auto A1 = srk(mixture, 500 /*psia*/, 50 + 460 /*R*/);
    auto A2 = srk(mixture, 400 /*psia*/, 50 + 460 /*R*/);
    auto A3 = srk(mixture, 300 /*psia*/, 50 + 460 /*R*/);
    auto A4 = srk(mixture, 500 /*psia*/, 100 + 460 /*R*/);
    auto A5 = srk(mixture, 400 /*psia*/, 100 + 460 /*R*/);
    auto A6 = srk(mixture, 300 /*psia*/, 100 + 460 /*R*/);

    CHECK(A1.size() == 1); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(A1[0],4))==0.1634);
    CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(A2[0],4))==0.6240);
    CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(A2[1],4))==0.1390);
    CHECK(A3.size() == 1); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(A3[0],4))==0.7529);
    CHECK(A4.size() == 1); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(A4[0],4))==0.6826);
    CHECK(A5.size() == 1); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(A5[0],4))==0.7614);
    CHECK(A6.size() == 1); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(A6[0],4))==0.8296);
}
