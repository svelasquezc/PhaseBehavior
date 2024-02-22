#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <PhaseBehavior/Utilities/Types.hpp>
#include <PhaseBehavior/Utilities/Input.hpp>
#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/BlackOil/CVD.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>

using NP_t = PhaseBehavior::Types::NumericalPrecision;

struct PVTTableTest{
    NP_t P;
    NP_t Bo;
    NP_t Bg;
    NP_t Rs;
    NP_t Rv;
};


TEST_CASE("Can calculate Black Oil Properties for a Mixture", "[BlackOilEstimation]"){

    SECTION("Using PR at 3500 psia and 250 F"){
        auto mixture = PhaseBehavior::Input::createMixtureFromFile<true>("PVTPR.csv", "InteractionCoefficients.csv");

        auto BOProperties = PhaseBehavior::BlackOil::BlackOilPropertiesEstimation<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 3500, 250 + 460, 300, 62 + 460, 20, 60 + 460, 14.7, 60 + 460, 100);

        
    }

}