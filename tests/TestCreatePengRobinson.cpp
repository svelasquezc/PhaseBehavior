#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <PhaseBehavior/Utilities/Input.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>

TEST_CASE("Can create and evaluate an EoS", "[PengRobinson]"){
    auto mixture = PhaseBehavior::Input::createMixtureFromFile("PVT.csv", "InteractionCoefficients.csv");

    PhaseBehavior::EoS::PR::PengRobinson pr;

    auto A1 = pr(mixture, 500 /*psia*/, 50 + 459.67 /*R*/);
    auto A2 = pr(mixture, 400 /*psia*/, 50 + 459.67 /*R*/);
    auto A3 = pr(mixture, 300 /*psia*/, 50 + 459.67 /*R*/);
    auto A4 = pr(mixture, 500 /*psia*/, 100 + 459.67 /*R*/);
    auto A5 = pr(mixture, 400 /*psia*/, 100 + 459.67 /*R*/);
    auto A6 = pr(mixture, 300 /*psia*/, 100 + 459.67 /*R*/);

    CHECK(A1.size() == 1); CHECK(Catch::Approx(A1[0])==0.1480);
    CHECK(Catch::Approx(A2[0])==0.6059); CHECK(Catch::Approx(A2[1])==0.1267);
    CHECK(A3.size() == 1); CHECK(Catch::Approx(A3[0])==0.7368);
    CHECK(A4.size() == 1); CHECK(Catch::Approx(A4[0])==0.6609);
    CHECK(A5.size() == 1); CHECK(Catch::Approx(A5[0])==0.7420);
    CHECK(A6.size() == 1); CHECK(Catch::Approx(A6[0])==0.8141);
}
