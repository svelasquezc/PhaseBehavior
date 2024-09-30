#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <PhaseBehavior/Utilities/Types.hpp>
#include <PhaseBehavior/Utilities/Math.hpp>
#include <PhaseBehavior/Utilities/Input.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>
#include <PhaseBehavior/PhaseEquilibrium.hpp>
#include <PhaseBehavior/Phase.hpp>

TEST_CASE("Can calculate Phase properties using different EoS", "[Phase]"){

    SECTION("Using PR at 500 psi and 50 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTPR.csv", "InteractionCoefficients.csv");

        auto [vaporEoS, liquidEoS, unused] = PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 500 /*psia*/, 50 + 460 /*R*/);

        PhaseBehavior::Phase::VaporLikePhase vaporPhase(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquidPhase(mixture);

        auto fng = mixture.molarFraction("vapor");
        auto Z_liq = mixture.compressibility("liquid");
        auto Z_vap = mixture.compressibility("vapor");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7233);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.molecularWeight(),3))==20.241);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.molecularWeight(),3))==68.141);

        vaporPhase.molarVolume(Z_vap, 500, 50+460);
        liquidPhase.molarVolume(Z_liq, 500, 50+460);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.density(),3))==2.162);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.density(),3))==38.067);

        vaporPhase.viscosity(mixture, Z_vap, 500, 50 + 460);
        liquidPhase.viscosity(mixture, Z_liq, 500, 50 + 460);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.viscosity(),3))==0.011);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.viscosity(),3))==0.218);

        vaporPhase.molarVolume(mixture, Z_vap, 500, 50+460, vaporEoS);
        liquidPhase.molarVolume(mixture, Z_liq, 500, 50+460, liquidEoS);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.density(),3))==2.146);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.density(),3))==38.008);

    }

    SECTION("Using PR at 400 psi and 50 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTPR.csv", "InteractionCoefficients.csv");

        auto [vaporEoS, liquidEoS, unused] = PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 400 /*psia*/, 50 + 460 /*R*/);

        PhaseBehavior::Phase::VaporLikePhase vaporPhase(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquidPhase(mixture);

        auto fng = mixture.molarFraction("vapor");
        auto Z_liq = mixture.compressibility("liquid");
        auto Z_vap = mixture.compressibility("vapor");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7521);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.molecularWeight(),3))==20.550);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.molecularWeight(),3))==72.774);

        vaporPhase.molarVolume (Z_vap, 400, 50+460);
        liquidPhase.molarVolume(Z_liq, 400, 50+460);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.density(),3))==1.706);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.density(),3))==38.745);

        vaporPhase.viscosity (mixture, Z_vap, 400, 50 + 460);
        liquidPhase.viscosity(mixture, Z_liq, 400, 50 + 460);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.viscosity(),3))==0.011);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.viscosity(),3))==0.237);

        vaporPhase.molarVolume (mixture, Z_vap, 400, 50+460, vaporEoS);
        liquidPhase.molarVolume(mixture, Z_liq, 400, 50+460, liquidEoS);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.density(),3))==1.696);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.density(),3))==38.855);

    }

    SECTION("Using PR at 300 psi and 50 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTPR.csv", "InteractionCoefficients.csv");

        auto [vaporEoS, liquidEoS, unused] = PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 300 /*psia*/, 50 + 460 /*R*/);

        PhaseBehavior::Phase::VaporLikePhase vaporPhase(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquidPhase(mixture);

        auto fng = mixture.molarFraction("vapor");
        auto Z_liq = mixture.compressibility("liquid");
        auto Z_vap = mixture.compressibility("vapor");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7831);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.molecularWeight(),3))==20.979);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.molecularWeight(),3))==78.673);

        vaporPhase.molarVolume (Z_vap, 300, 50+460);
        liquidPhase.molarVolume(Z_liq, 300, 50+460);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.density(),3))==1.269);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.density(),3))==39.467);

        vaporPhase.viscosity (mixture, Z_vap, 300, 50 + 460);
        liquidPhase.viscosity(mixture, Z_liq, 300, 50 + 460);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.viscosity(),3))==0.011);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.viscosity(),3))==0.259);

        vaporPhase.molarVolume (mixture, Z_vap, 300, 50+460, vaporEoS);
        liquidPhase.molarVolume(mixture, Z_liq, 300, 50+460, liquidEoS);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.density(),3))==1.263);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.density(),3))==39.797);
    }

    SECTION("Using PR at 500 psi and 100 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTPR.csv", "InteractionCoefficients.csv");

        auto [vaporEoS, liquidEoS, unused] = PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 500 /*psia*/, 100 + 460 /*R*/);

        PhaseBehavior::Phase::VaporLikePhase vaporPhase(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquidPhase(mixture);

        auto fng = mixture.molarFraction("vapor");
        auto Z_liq = mixture.compressibility("liquid");
        auto Z_vap = mixture.compressibility("vapor");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7941);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.molecularWeight(),3))==21.745);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.molecularWeight(),3))==78.814);

        vaporPhase.molarVolume (Z_vap, 500, 100+460);
        liquidPhase.molarVolume(Z_liq, 500, 100+460);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.density(),3))==2.062);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.density(),3))==37.920);

        vaporPhase.viscosity (mixture, Z_vap, 500, 100 + 460);
        liquidPhase.viscosity(mixture, Z_liq, 500, 100 + 460);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.viscosity(),3))==0.012);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.viscosity(),3))==0.200);

        vaporPhase.molarVolume (mixture, Z_vap, 500, 100+460, vaporEoS);
        liquidPhase.molarVolume(mixture, Z_liq, 500, 100+460, liquidEoS);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.density(),3))==2.048);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.density(),3))==38.311);
    }

    SECTION("Using PR at 400 psi and 100 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTPR.csv", "InteractionCoefficients.csv");

        auto [vaporEoS, liquidEoS, unused] = PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 400 /*psia*/, 100 + 460 /*R*/);

        PhaseBehavior::Phase::VaporLikePhase vaporPhase(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquidPhase(mixture);

        auto fng = mixture.molarFraction("vapor");
        auto Z_liq = mixture.compressibility("liquid");
        auto Z_vap = mixture.compressibility("vapor");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.8155);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.molecularWeight(),3))==22.056);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.molecularWeight(),3))==84.067);

        vaporPhase.molarVolume (Z_vap, 400, 100+460);
        liquidPhase.molarVolume(Z_liq, 400, 100+460);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.density(),3))==1.634);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.density(),3))==38.537);

        vaporPhase.viscosity (mixture, Z_vap, 400, 100 + 460);
        liquidPhase.viscosity(mixture, Z_liq, 400, 100 + 460);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.viscosity(),3))==0.012);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.viscosity(),3))==0.215);

        vaporPhase.molarVolume (mixture, Z_vap, 400, 100+460, vaporEoS);
        liquidPhase.molarVolume(mixture, Z_liq, 400, 100+460, liquidEoS);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.density(),3))==1.626);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.density(),3))==39.113);

    }

    SECTION("Using PR at 300 psi and 100 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTPR.csv", "InteractionCoefficients.csv");

        auto [vaporEoS, liquidEoS, unused] = PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 300 /*psia*/, 100 + 460 /*R*/);

        PhaseBehavior::Phase::VaporLikePhase vaporPhase(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquidPhase(mixture);

        auto fng = mixture.molarFraction("vapor");
        auto Z_liq = mixture.compressibility("liquid");
        auto Z_vap = mixture.compressibility("vapor");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.8383);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.molecularWeight(),3))==22.468);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.molecularWeight(),3))==90.669);

        vaporPhase.molarVolume (Z_vap, 300, 100+460);
        liquidPhase.molarVolume(Z_liq, 300, 100+460);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.density(),3))==1.219);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.density(),3))==39.187);

        vaporPhase.viscosity (mixture, Z_vap, 300, 100+460);
        liquidPhase.viscosity(mixture, Z_liq, 300, 100+460);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.viscosity(),3)+0.001)==0.012);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.viscosity(),3))==0.232);

        vaporPhase.molarVolume (mixture, Z_vap, 300, 100+460, vaporEoS);
        liquidPhase.molarVolume(mixture, Z_liq, 300, 100+460, liquidEoS);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.density(),3))==1.214);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.density(),3))==39.998);
    }

    SECTION("Using SRK at 500 psi and 50 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTSRK.csv", "InteractionCoefficients.csv");

        auto [vaporEoS, liquidEoS, unused] = PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 500 /*psia*/, 50 + 460 /*R*/);

        PhaseBehavior::Phase::VaporLikePhase vaporPhase(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquidPhase(mixture);

        auto fng = mixture.molarFraction("vapor");
        auto Z_liq = mixture.compressibility("liquid");
        auto Z_vap = mixture.compressibility("vapor");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7245);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.molecularWeight(),3))==20.207);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.molecularWeight(),3))==68.430);

        vaporPhase.molarVolume(Z_vap, 500, 50+460);
        liquidPhase.molarVolume(Z_liq, 500, 50+460);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.density(),3))==2.108);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.density(),3))==33.853);

        vaporPhase.viscosity(mixture, Z_vap, 500, 50 + 460);
        liquidPhase.viscosity(mixture, Z_liq, 500, 50 + 460);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.viscosity(),3))==0.011);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.viscosity(),3))==0.114);

        vaporPhase.molarVolume(mixture, Z_vap, 500, 50+460, vaporEoS);
        liquidPhase.molarVolume(mixture, Z_liq, 500, 50+460, liquidEoS);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.density(),3))==2.113);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.density(),3))==38.527);
    }

    SECTION("Using SRK at 400 psi and 50 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTSRK.csv", "InteractionCoefficients.csv");

        auto [vaporEoS, liquidEoS, unused] = PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 400 /*psia*/, 50 + 460 /*R*/);

        PhaseBehavior::Phase::VaporLikePhase vaporPhase(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquidPhase(mixture);

        auto fng = mixture.molarFraction("vapor");
        auto Z_liq = mixture.compressibility("liquid");
        auto Z_vap = mixture.compressibility("vapor");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7530);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.molecularWeight(),3))==20.522);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.molecularWeight(),3))==73.039);

        vaporPhase.molarVolume (Z_vap, 400, 50+460);
        liquidPhase.molarVolume(Z_liq, 400, 50+460);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.density(),3))==1.672);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.density(),3))==34.448);

        vaporPhase.viscosity (mixture, Z_vap, 400, 50 + 460);
        liquidPhase.viscosity(mixture, Z_liq, 400, 50 + 460);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.viscosity(),3))==0.011);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.viscosity(),3))==0.121);

        vaporPhase.molarVolume (mixture, Z_vap, 400, 50+460, vaporEoS);
        liquidPhase.molarVolume(mixture, Z_liq, 400, 50+460, liquidEoS);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.density(),3))==1.675);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.density(),3))==39.446);
    }

    SECTION("Using SRK at 300 psi and 50 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTSRK.csv", "InteractionCoefficients.csv");

        auto [vaporEoS, liquidEoS, unused] = PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 300 /*psia*/, 50 + 460 /*R*/);

        PhaseBehavior::Phase::VaporLikePhase vaporPhase(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquidPhase(mixture);

        auto fng = mixture.molarFraction("vapor");
        auto Z_liq = mixture.compressibility("liquid");
        auto Z_vap = mixture.compressibility("vapor");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7836);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.molecularWeight(),3))==20.955);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.molecularWeight(),3))==78.900);

        vaporPhase.molarVolume (Z_vap, 300, 50+460);
        liquidPhase.molarVolume(Z_liq, 300, 50+460);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.density(),3))==1.249);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.density(),3))==35.087);

        vaporPhase.viscosity (mixture, Z_vap, 300, 50 + 460);
        liquidPhase.viscosity(mixture, Z_liq, 300, 50 + 460);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.viscosity(),3))==0.011);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.viscosity(),3))==0.128);

        vaporPhase.molarVolume (mixture, Z_vap, 300, 50+460, vaporEoS);
        liquidPhase.molarVolume(mixture, Z_liq, 300, 50+460, liquidEoS);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.density(),3))==1.251);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.density(),3))==40.477);
    }

    SECTION("Using SRK at 500 psi and 100 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTSRK.csv", "InteractionCoefficients.csv");

        auto [vaporEoS, liquidEoS, unused] = PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 500 /*psia*/, 100 + 460 /*R*/);

        PhaseBehavior::Phase::VaporLikePhase vaporPhase(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquidPhase(mixture);

        auto fng = mixture.molarFraction("vapor");
        auto Z_liq = mixture.compressibility("liquid");
        auto Z_vap = mixture.compressibility("vapor");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7942);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.molecularWeight(),3))==21.705);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.molecularWeight(),3))==78.982);

        vaporPhase.molarVolume (Z_vap, 500, 100+460);
        liquidPhase.molarVolume(Z_liq, 500, 100+460);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.density(),3))==2.014);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.density(),3))==33.692);

        vaporPhase.viscosity (mixture, Z_vap, 500, 100 + 460);
        liquidPhase.viscosity(mixture, Z_liq, 500, 100 + 460);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.viscosity(),3))==0.012);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.viscosity(),3))==0.108);

        vaporPhase.molarVolume (mixture, Z_vap, 500, 100+460, vaporEoS);
        liquidPhase.molarVolume(mixture, Z_liq, 500, 100+460, liquidEoS);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.density(),3))==2.018);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.density(),3))==38.733);
    }

    SECTION("Using SRK at 400 psi and 100 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTSRK.csv", "InteractionCoefficients.csv");

        auto [vaporEoS, liquidEoS, unused] = PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 400 /*psia*/, 100 + 460 /*R*/);

        PhaseBehavior::Phase::VaporLikePhase vaporPhase(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquidPhase(mixture);

        auto fng = mixture.molarFraction("vapor");
        auto Z_liq = mixture.compressibility("liquid");
        auto Z_vap = mixture.compressibility("vapor");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.8156);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.molecularWeight(),3))==22.023);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.molecularWeight(),3))==84.230);

        vaporPhase.molarVolume (Z_vap, 400, 100+460);
        liquidPhase.molarVolume(Z_liq, 400, 100+460);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.density(),3))==1.603);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.density(),3))==34.237);

        vaporPhase.viscosity (mixture, Z_vap, 400, 100 + 460);
        liquidPhase.viscosity(mixture, Z_liq, 400, 100 + 460);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.viscosity(),3))==0.012);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.viscosity(),3))==0.114);

        vaporPhase.molarVolume (mixture, Z_vap, 400, 100+460, vaporEoS);
        liquidPhase.molarVolume(mixture, Z_liq, 400, 100+460, liquidEoS);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.density(),3))==1.606);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.density(),3))==39.610);
    }

    SECTION("Using SRK at 300 psi and 100 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTSRK.csv", "InteractionCoefficients.csv");

        auto [vaporEoS, liquidEoS, unused] = PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 300 /*psia*/, 100 + 460 /*R*/);

        PhaseBehavior::Phase::VaporLikePhase vaporPhase(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquidPhase(mixture);

        auto fng = mixture.molarFraction("vapor");
        auto Z_liq = mixture.compressibility("liquid");
        auto Z_vap = mixture.compressibility("vapor");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.8384);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.molecularWeight(),3))==22.442);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.molecularWeight(),3))==90.819);

        vaporPhase.molarVolume (Z_vap, 300, 100+460);
        liquidPhase.molarVolume(Z_liq, 300, 100+460);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.density(),3))==1.201);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.density(),3))==34.815);

        vaporPhase.viscosity (mixture, Z_vap, 300, 100 + 460);
        liquidPhase.viscosity(mixture, Z_liq, 300, 100 + 460);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.viscosity(),3)+0.001)==0.012);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.viscosity(),3))==0.119);

        vaporPhase.molarVolume (mixture, Z_vap, 300, 100+460, vaporEoS);
        liquidPhase.molarVolume(mixture, Z_liq, 300, 100+460, liquidEoS);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporPhase.density(),3))==1.203);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidPhase.density(),3))==40.587);
    }
}
