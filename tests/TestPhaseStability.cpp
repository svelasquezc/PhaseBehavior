#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <PhaseBehavior/Utilities/Types.hpp>
#include <PhaseBehavior/Utilities/Math.hpp>
#include <PhaseBehavior/Utilities/Input.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>
#include <PhaseBehavior/PhaseEquilibrium.hpp>
#include <PhaseBehavior/Phase.hpp>

TEST_CASE("Can calculate Phase properties using different EoS", "[Phase]"){

    SECTION("Using PR at 3000 psi and 250 F"){
        {
            constexpr bool hasShift = true;

            auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTPR.csv", "InteractionCoefficients.csv");

            auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 3000 /*psia*/, 250 + 460 /*R*/);

            REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Stable);

            for (auto& mixtureComponent : mixture){
                mixtureComponent.composition("vapor", mixtureComponent.composition());
            }

            auto fng = 1;

            auto Z_vap = mixture.compressibility("global");

            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==1.0000);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.7486);

            PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==33.49);
        }
    }

    SECTION("Using PR at 2900 psi and 250 F"){
        {
            constexpr bool hasShift = true;

            auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTPR.csv", "InteractionCoefficients.csv");

            auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 2900 /*psia*/, 250 + 460 /*R*/);

            REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Stable);

            for (auto& mixtureComponent : mixture){
                mixtureComponent.composition("vapor", mixtureComponent.composition());
            }

            auto fng = 1;
            
            auto Z_vap = mixture.compressibility("global");

            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==1.0000);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.7393);

            PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==33.49);
        }
    }

    SECTION("Using PR at 2800 psi and 250 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTPR.csv", "InteractionCoefficients.csv");

        auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 2800 /*psia*/, 250 + 460 /*R*/);

        REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Stable);

        for (auto& mixtureComponent : mixture){
            mixtureComponent.composition("vapor", mixtureComponent.composition());
        }

        auto fng = 1;

        auto Z_vap = mixture.compressibility("global");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==1.0000);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.7303);

        PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==33.49);

    }

    SECTION("Using PR at 2700 psi and 250 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTPR.csv", "InteractionCoefficients.csv");

        auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 2700 /*psia*/, 250 + 460 /*R*/);

        REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable);

        PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 2700 /*psia*/, 250 + 460 /*R*/, false);

        auto fng = mixture.molarFraction("vapor");

        auto Z_vap = mixture.compressibility("vapor");
        auto Z_liq = mixture.compressibility("liquid");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.8055);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.7426);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.6738);

        PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquid(mixture);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==31.37);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquid.molecularWeight(),2))==42.31);
        
    }

    SECTION("Using PR at 2600 psi and 250 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTPR.csv", "InteractionCoefficients.csv");

        auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 2600 /*psia*/, 250 + 460 /*R*/);

        REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable);

        PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 2600 /*psia*/, 250 + 460 /*R*/, false);

        auto fng = mixture.molarFraction("vapor");

        auto Z_vap = mixture.compressibility("vapor");
        auto Z_liq = mixture.compressibility("liquid");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7719);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.7520);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.6539);

        PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquid(mixture);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==30.04);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquid.molecularWeight(),2))==45.19);
        
    }

    SECTION("Using PR at 2500 psi and 250 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTPR.csv", "InteractionCoefficients.csv");

        auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 2500 /*psia*/, 250 + 460 /*R*/);

        REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable);

        PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 2500 /*psia*/, 250 + 460 /*R*/, false);

        auto fng = mixture.molarFraction("vapor");

        auto Z_vap = mixture.compressibility("vapor");
        auto Z_liq = mixture.compressibility("liquid");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7645);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.7594);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.6356);

        PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquid(mixture);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==29.11);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquid.molecularWeight(),2))==47.72);
        
    }

    SECTION("Using PR at 2400 psi and 250 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTPR.csv", "InteractionCoefficients.csv");

        auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 2400 /*psia*/, 250 + 460 /*R*/);

        REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable);

        PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 2400 /*psia*/, 250 + 460 /*R*/, false);

        auto fng = mixture.molarFraction("vapor");

        auto Z_vap = mixture.compressibility("vapor");
        auto Z_liq = mixture.compressibility("liquid");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7652);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.7661);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.6179);

        PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquid(mixture);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==28.40);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquid.molecularWeight(),2))==50.09);
        
    }

    SECTION("Using PR at 2300 psi and 250 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTPR.csv", "InteractionCoefficients.csv");

        auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 2300 /*psia*/, 250 + 460 /*R*/);

        REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable);

        PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 2300 /*psia*/, 250 + 460 /*R*/, false);

        auto fng = mixture.molarFraction("vapor");

        auto Z_vap = mixture.compressibility("vapor");
        auto Z_liq = mixture.compressibility("liquid");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7692);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.7724);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.6003);

        PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquid(mixture);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==27.83);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquid.molecularWeight(),2))==52.38);
        
    }

    SECTION("Using PR at 2200 psi and 250 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTPR.csv", "InteractionCoefficients.csv");

        auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 2200 /*psia*/, 250 + 460 /*R*/);

        REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable);

        PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 2200 /*psia*/, 250 + 460 /*R*/, false);

        auto fng = mixture.molarFraction("vapor");

        auto Z_vap = mixture.compressibility("vapor");
        auto Z_liq = mixture.compressibility("liquid");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7750);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.7786);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.5826);

        PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquid(mixture);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==27.36);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquid.molecularWeight(),2))==54.64);
        
    }

    SECTION("Using PR at 2100 psi and 250 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTPR.csv", "InteractionCoefficients.csv");

        auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 2100 /*psia*/, 250 + 460 /*R*/);

        REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable);

        PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 2100 /*psia*/, 250 + 460 /*R*/, false);

        auto fng = mixture.molarFraction("vapor");

        auto Z_vap = mixture.compressibility("vapor");
        auto Z_liq = mixture.compressibility("liquid");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7816);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.7848);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.5648);

        PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquid(mixture);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==26.96);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquid.molecularWeight(),2))==56.88);
        
    }

    SECTION("Using PR at 2000 psi and 250 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTPR.csv", "InteractionCoefficients.csv");

        auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 2000 /*psia*/, 250 + 460 /*R*/);

        REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable);

        PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 2000 /*psia*/, 250 + 460 /*R*/, false);

        auto fng = mixture.molarFraction("vapor");

        auto Z_vap = mixture.compressibility("vapor");
        auto Z_liq = mixture.compressibility("liquid");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7887);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.7911);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.5466);

        PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquid(mixture);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==26.63);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquid.molecularWeight(),2))==59.13);
        
    }

    SECTION("Using PR at 1900 psi and 250 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTPR.csv", "InteractionCoefficients.csv");

        auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 1900 /*psia*/, 250 + 460 /*R*/);

        REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable);

        PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 1900 /*psia*/, 250 + 460 /*R*/, false);

        auto fng = mixture.molarFraction("vapor");

        auto Z_vap = mixture.compressibility("vapor");
        auto Z_liq = mixture.compressibility("liquid");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7961);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.7975);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.5280);

        PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquid(mixture);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==26.34);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquid.molecularWeight(),2))==61.42);
        
    }

    SECTION("Using PR at 1800 psi and 250 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTPR.csv", "InteractionCoefficients.csv");

        auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 1800 /*psia*/, 250 + 460 /*R*/);

        REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable);

        PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 1800 /*psia*/, 250 + 460 /*R*/, false);

        auto fng = mixture.molarFraction("vapor");

        auto Z_vap = mixture.compressibility("vapor");
        auto Z_liq = mixture.compressibility("liquid");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.8036);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.8042);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.5089);

        PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquid(mixture);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==26.10);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquid.molecularWeight(),2))==63.75);
        
    }

    SECTION("Using PR at 1700 psi and 250 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTPR.csv", "InteractionCoefficients.csv");

        auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 1700 /*psia*/, 250 + 460 /*R*/);

        REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable);

        PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 1700 /*psia*/, 250 + 460 /*R*/, false);

        auto fng = mixture.molarFraction("vapor");

        auto Z_vap = mixture.compressibility("vapor");
        auto Z_liq = mixture.compressibility("liquid");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.8112);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.8111);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.4894);

        PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquid(mixture);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==25.90);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquid.molecularWeight(),2))==66.13);
        
    }

    SECTION("Using PR at 1600 psi and 250 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTPR.csv", "InteractionCoefficients.csv");

        auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 1600 /*psia*/, 250 + 460 /*R*/);

        REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable);

        PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 1600 /*psia*/, 250 + 460 /*R*/, false);

        auto fng = mixture.molarFraction("vapor");

        auto Z_vap = mixture.compressibility("vapor");
        auto Z_liq = mixture.compressibility("liquid");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.8188);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.8183);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.4692);

        PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquid(mixture);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==25.73);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquid.molecularWeight(),2))==68.59);
        
    }

    SECTION("Using SRK at 3000 psi and 250 F"){
        {
            constexpr bool hasShift = true;

            auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTSRK.csv", "InteractionCoefficients.csv");

            auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 3000 /*psia*/, 250 + 460 /*R*/);

            REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Stable);

            for (auto& mixtureComponent : mixture){
                mixtureComponent.composition("vapor", mixtureComponent.composition());
            }

            auto fng = 1;

            auto Z_vap = mixture.compressibility("global");

            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==1.0000);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.8061);

            PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==33.49);
        }
    }

    SECTION("Using SRK at 2900 psi and 250 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTSRK.csv", "InteractionCoefficients.csv");

        auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 2900 /*psia*/, 250 + 460 /*R*/);

        REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable);

        PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 2900 /*psia*/, 250 + 460 /*R*/, false);

        auto fng = mixture.molarFraction("vapor");

        auto Z_vap = mixture.compressibility("vapor");
        auto Z_liq = mixture.compressibility("liquid");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.8719);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.7989);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.7766);

        PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquid(mixture);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==32.95);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquid.molecularWeight(),2))==37.17);       
    }

    SECTION("Using SRK at 2800 psi and 250 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTSRK.csv", "InteractionCoefficients.csv");

        auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 2800 /*psia*/, 250 + 460 /*R*/);

        REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable);

        PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 2800 /*psia*/, 250 + 460 /*R*/, false);

        auto fng = mixture.molarFraction("vapor");

        auto Z_vap = mixture.compressibility("vapor");
        auto Z_liq = mixture.compressibility("liquid");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7196);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.8102);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.7526);

        PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquid(mixture);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==30.57);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquid.molecularWeight(),2))==41.00);
    }

    SECTION("Using SRK at 2700 psi and 250 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTSRK.csv", "InteractionCoefficients.csv");

        auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 2700 /*psia*/, 250 + 460 /*R*/);

        REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable);

        PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 2700 /*psia*/, 250 + 460 /*R*/, false);

        auto fng = mixture.molarFraction("vapor");

        auto Z_vap = mixture.compressibility("vapor");
        auto Z_liq = mixture.compressibility("liquid");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7118);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.8160);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.7337);

        PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquid(mixture);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==29.36);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquid.molecularWeight(),2))==43.71);
    }

    SECTION("Using SRK at 2600 psi and 250 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTSRK.csv", "InteractionCoefficients.csv");

        auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 2600 /*psia*/, 250 + 460 /*R*/);

        REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable);

        PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 2600 /*psia*/, 250 + 460 /*R*/, false);

        auto fng = mixture.molarFraction("vapor");

        auto Z_vap = mixture.compressibility("vapor");
        auto Z_liq = mixture.compressibility("liquid");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7167);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.8205);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.7159);

        PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquid(mixture);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==28.51);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquid.molecularWeight(),2))==46.09);
        
    }

    SECTION("Using SRK at 2500 psi and 250 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTSRK.csv", "InteractionCoefficients.csv");

        auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 2500 /*psia*/, 250 + 460 /*R*/);

        REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable);

        PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 2500 /*psia*/, 250 + 460 /*R*/, false);

        auto fng = mixture.molarFraction("vapor");

        auto Z_vap = mixture.compressibility("vapor");
        auto Z_liq = mixture.compressibility("liquid");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7250);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.8245);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.6983);

        PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquid(mixture);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==27.87);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquid.molecularWeight(),2))==48.34);
        
    }

    SECTION("Using SRK at 2400 psi and 250 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTSRK.csv", "InteractionCoefficients.csv");

        auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 2400 /*psia*/, 250 + 460 /*R*/);

        REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable);

        PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 2400 /*psia*/, 250 + 460 /*R*/, false);

        auto fng = mixture.molarFraction("vapor");

        auto Z_vap = mixture.compressibility("vapor");
        auto Z_liq = mixture.compressibility("liquid");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7345);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.8282);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.6807);

        PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquid(mixture);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==27.34);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquid.molecularWeight(),2))==50.51);
        
    }

    SECTION("Using SRK at 2300 psi and 250 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTSRK.csv", "InteractionCoefficients.csv");

        auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 2300 /*psia*/, 250 + 460 /*R*/);

        REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable);

        PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 2300 /*psia*/, 250 + 460 /*R*/, false);

        auto fng = mixture.molarFraction("vapor");

        auto Z_vap = mixture.compressibility("vapor");
        auto Z_liq = mixture.compressibility("liquid");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7443);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.8319);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.6628);

        PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquid(mixture);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==26.91);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquid.molecularWeight(),2))==52.65);
        
    }

    SECTION("Using SRK at 2200 psi and 250 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTSRK.csv", "InteractionCoefficients.csv");

        auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 2200 /*psia*/, 250 + 460 /*R*/);

        REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable);

        PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 2200 /*psia*/, 250 + 460 /*R*/, false);

        auto fng = mixture.molarFraction("vapor");

        auto Z_vap = mixture.compressibility("vapor");
        auto Z_liq = mixture.compressibility("liquid");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7542);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.8357);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.6445);

        PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquid(mixture);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==26.55);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquid.molecularWeight(),2))==54.79);
        
    }

    SECTION("Using SRK at 2100 psi and 250 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTSRK.csv", "InteractionCoefficients.csv");

        auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 2100 /*psia*/, 250 + 460 /*R*/);

        REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable);

        PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 2100 /*psia*/, 250 + 460 /*R*/, false);

        auto fng = mixture.molarFraction("vapor");

        auto Z_vap = mixture.compressibility("vapor");
        auto Z_liq = mixture.compressibility("liquid");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7639);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.8396);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.6258);

        PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquid(mixture);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==26.25);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquid.molecularWeight(),2))==56.94);
        
    }

    SECTION("Using SRK at 2000 psi and 250 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTSRK.csv", "InteractionCoefficients.csv");

        auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 2000 /*psia*/, 250 + 460 /*R*/);

        REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable);

        PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 2000 /*psia*/, 250 + 460 /*R*/, false);

        auto fng = mixture.molarFraction("vapor");

        auto Z_vap = mixture.compressibility("vapor");
        auto Z_liq = mixture.compressibility("liquid");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7735);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.8437);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.6065);

        PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquid(mixture);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==25.99);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquid.molecularWeight(),2))==59.12);
        
    }

    SECTION("Using SRK at 1900 psi and 250 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTSRK.csv", "InteractionCoefficients.csv");

        auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 1900 /*psia*/, 250 + 460 /*R*/);

        REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable);

        PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 1900 /*psia*/, 250 + 460 /*R*/, false);

        auto fng = mixture.molarFraction("vapor");

        auto Z_vap = mixture.compressibility("vapor");
        auto Z_liq = mixture.compressibility("liquid");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7830);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.8480);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.5867);

        PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquid(mixture);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==25.78);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquid.molecularWeight(),2))==61.34);
        
    }

    SECTION("Using SRK at 1800 psi and 250 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTSRK.csv", "InteractionCoefficients.csv");

        auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 1800 /*psia*/, 250 + 460 /*R*/);

        REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable);

        PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 1800 /*psia*/, 250 + 460 /*R*/, false);

        auto fng = mixture.molarFraction("vapor");

        auto Z_vap = mixture.compressibility("vapor");
        auto Z_liq = mixture.compressibility("liquid");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7922);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.8525);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.5662);

        PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquid(mixture);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==25.59);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquid.molecularWeight(),2))==63.62);
        
    }

    SECTION("Using SRK at 1700 psi and 250 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTSRK.csv", "InteractionCoefficients.csv");

        auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 1700 /*psia*/, 250 + 460 /*R*/);

        REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable);

        PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 1700 /*psia*/, 250 + 460 /*R*/, false);

        auto fng = mixture.molarFraction("vapor");

        auto Z_vap = mixture.compressibility("vapor");
        auto Z_liq = mixture.compressibility("liquid");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.8013);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.8573);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.5451);

        PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquid(mixture);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==25.44);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquid.molecularWeight(),2))==65.96);
    }

    SECTION("Using SRK at 1600 psi and 250 F"){

        constexpr bool hasShift = true;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile<hasShift>("PVTSRK.csv", "InteractionCoefficients.csv");

        auto stabilityResult = PhaseBehavior::VaporLiquidEquilibrium::phaseStability<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 1600 /*psia*/, 250 + 460 /*R*/);

        REQUIRE(stabilityResult==PhaseBehavior::VaporLiquidEquilibrium::PhaseStabilityResult::Unstable);

        PhaseBehavior::VaporLiquidEquilibrium::succesiveSubstitution<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 1600 /*psia*/, 250 + 460 /*R*/, false);

        auto fng = mixture.molarFraction("vapor");

        auto Z_vap = mixture.compressibility("vapor");
        auto Z_liq = mixture.compressibility("liquid");

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.8102);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.8623);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.5233);

        PhaseBehavior::Phase::VaporLikePhase vapor(mixture);
        PhaseBehavior::Phase::LiquidLikePhase liquid(mixture);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vapor.molecularWeight(),2))==25.32);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquid.molecularWeight(),2))==68.39);
        
    }
}