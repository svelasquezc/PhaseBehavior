#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include<PhaseBehavior/GERG/BinaryReducingParameters.hpp>

TEST_CASE("Can access to the value regardless of key permutation", "[GERG]"){
    using namespace PhaseBehavior::EoS::GERG::Coefficients::ReducingParameters;

    REQUIRE(Density::beta.at({"CH4", "N2"}) == Density::beta.at({"N2", "CH4"}));

}