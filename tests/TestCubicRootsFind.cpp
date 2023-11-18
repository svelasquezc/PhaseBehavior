#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include<PhaseBehavior/Utilities/Types.hpp>
#include <PhaseBehavior/Utilities/Math.hpp>

TEST_CASE("cubic roots find function can solve cubic polynomials", "[cubicRootsFind]"){

    SECTION("Finds three roots in a polynomial"){
        auto roots = Math::cubicRootsFind<TypesDefinition::NumericalPrecision>(-6,11,-6);

        CHECK(roots.size() == 3);

        CHECK(roots[0] == Catch::Approx(1));
        CHECK(roots[1] == Catch::Approx(2));
        CHECK(roots[2] == Catch::Approx(3));
    }

    SECTION("Finds only one real root in a polynomial (two imaginary roots)"){
        auto roots = Math::cubicRootsFind<TypesDefinition::NumericalPrecision>(-7.8693,13.3771,-6.5354);
        REQUIRE(roots.size() == 1);

        REQUIRE(roots[0] == Catch::Approx(5.7357));
    }

}
