#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <cmath>

#include<PhaseBehavior/Utilities/Types.hpp>
#include <PhaseBehavior/Utilities/Math.hpp>

TEST_CASE("can find a root of a cubic polynomial", "[NewtonRaphson]"){

    SECTION("Finds only one real root in a polynomial (two imaginary roots)"){

        auto cubicPolynomial = [](const PhaseBehavior::Types::NumericalPrecision& x){
            return std::pow(x,3) -7.8693*std::pow(x,2) + 13.3771*x - 6.5354;             
        };

        auto divergenceCriterium = [](const PhaseBehavior::Types::NumericalPrecision& x){
            return std::abs(x)>10;
        };

        auto root = PhaseBehavior::Math::NewtonRaphson<PhaseBehavior::Types::NumericalPrecision>(5, cubicPolynomial, divergenceCriterium);

        REQUIRE(root == Catch::Approx(5.7357));
    }

}
