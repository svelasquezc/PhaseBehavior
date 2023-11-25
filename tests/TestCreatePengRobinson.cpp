#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <PhaseBehavior/Utilities/Input.hpp>
#include <PhaseBehavior/Utilities/Math.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>

TEST_CASE("Can create and evaluate an EoS", "[PengRobinson]"){
    auto mixture = PhaseBehavior::Input::createMixtureFromFile("PVT.csv", "InteractionCoefficients.csv");

    PhaseBehavior::EoS::PR::PengRobinson pr;

    auto A1 = pr(mixture, 500 /*psia*/, 50 + 460 /*R*/);
    auto A2 = pr(mixture, 400 /*psia*/, 50 + 460 /*R*/);
    auto A3 = pr(mixture, 300 /*psia*/, 50 + 460 /*R*/);
    auto A4 = pr(mixture, 500 /*psia*/, 100 + 460 /*R*/);
    auto A5 = pr(mixture, 400 /*psia*/, 100 + 460 /*R*/);
    auto A6 = pr(mixture, 300 /*psia*/, 100 + 460 /*R*/);

    CHECK(A1.size() == 1); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(A1[0],5))==0.1480);
    CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(A2[0],5))==0.6059);
    CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(A2[1],5))==0.1267);
    CHECK(A3.size() == 1); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(A3[0],5))==0.7368);
    CHECK(A4.size() == 1); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(A4[0],5))==0.6609);
    CHECK(A5.size() == 1); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(A5[0],5))==0.7420);
    CHECK(A6.size() == 1); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(A6[0],5))==0.8141);
}
