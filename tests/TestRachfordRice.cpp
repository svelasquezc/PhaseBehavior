#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <PhaseBehavior/Utilities/Input.hpp>
#include <PhaseBehavior/Utilities/Math.hpp>

#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/MixingRules.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>
#include <PhaseBehavior/PhaseEquilibrium.hpp>


TEST_CASE("Can obtain the compressibility, A, and B, for a two phase mixture", "[RachfordVLE]"){

    SECTION("Using PengRobinson"){

        using MixingRule = PhaseBehavior::EoS::MixingRules::RandomMixingRule<
                            PhaseBehavior::EoS::PR::constants::omegaA, 
                            PhaseBehavior::EoS::PR::constants::omegaB, 
                            PhaseBehavior::EoS::PR::PengRobinson
                            >;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile("PVT.csv", "InteractionCoefficients.csv");

        MixingRule mixingRule;

        PhaseBehavior::EoS::PR::PengRobinson pr;

        // 500 psia and 50 F
        {
            PhaseBehavior::VaporLiquidEquilibrium::rachfordVLE(mixture, 500 /*psia*/, 50 + 460 /*R*/);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture.molarFraction("vapor"),4))==0.7252);

            auto [vaporA,vaporB] = mixingRule(mixture, 500 /*psia*/, 50 + 460 /*R*/, "vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporA,4))==0.1780); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporB,4))==0.0442);

            pr(mixture, 500 /*psia*/, 50 + 460 /*R*/, "vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(pr.selectedCompressibility(),4))==0.8671);


            auto [liquidA,liquidB] = mixingRule(mixture, 500 /*psia*/, 50 + 460 /*R*/, "liquid");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidA,4))==1.6922); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidB,4))==0.1353);

            pr(mixture, 500 /*psia*/, 50 + 460 /*R*/, "liquid");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(pr.selectedCompressibility(),4))==0.1665);
        }
        // 400 psia and 50 F
        {
            PhaseBehavior::VaporLiquidEquilibrium::rachfordVLE(mixture, 400 /*psia*/, 50 + 460 /*R*/);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture.molarFraction("vapor"),4))==0.7533);
            auto [vaporA,vaporB] = mixingRule(mixture, 400 /*psia*/, 50 + 460 /*R*/, "vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporA,4))==0.1480); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporB,4))==0.0360);
            pr(mixture, 400 /*psia*/, 50 + 460 /*R*/, "vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(pr.selectedCompressibility(),4))==0.8878);
            auto [liquidA,liquidB] = mixingRule(mixture, 400 /*psia*/, 50 + 460 /*R*/, "liquid");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidA,4))==1.5146); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidB,4))==0.1148);
            pr(mixture, 400 /*psia*/, 50 + 460 /*R*/, "liquid");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(pr.selectedCompressibility(),4))==0.1395);
        }
        // 300 psia and 50 F
        {
            PhaseBehavior::VaporLiquidEquilibrium::rachfordVLE(mixture, 300 /*psia*/, 50 + 460 /*R*/);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture.molarFraction("vapor"),4))==0.7840);
            auto [vaporA,vaporB] = mixingRule(mixture, 300 /*psia*/, 50 + 460 /*R*/, "vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporA,4))==0.1165); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporB,4))==0.0275);
            pr(mixture, 300 /*psia*/, 50 + 460 /*R*/, "vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(pr.selectedCompressibility(),4))==0.9104);
            auto [liquidA,liquidB] = mixingRule(mixture, 300 /*psia*/, 50 + 460 /*R*/, "liquid");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidA,4))==1.3010); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidB,4))==0.0925);
            pr(mixture, 300 /*psia*/, 50 + 460 /*R*/, "liquid");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(pr.selectedCompressibility(),4))==0.1108);
        }
        // 500 psia and 100 F
        {
            PhaseBehavior::VaporLiquidEquilibrium::rachfordVLE(mixture, 500 /*psia*/, 100 + 460 /*R*/);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture.molarFraction("vapor"),4))==0.7954);
            auto [vaporA,vaporB] = mixingRule(mixture, 500 /*psia*/, 100 + 460 /*R*/, "vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporA,4))==0.1590); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporB,4))==0.0425);
            pr(mixture, 500 /*psia*/, 100 + 460 /*R*/, "vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(pr.selectedCompressibility(),4))==0.8863);
            auto [liquidA,liquidB] = mixingRule(mixture, 500 /*psia*/, 100 + 460 /*R*/, "liquid");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidA,4))==1.7703); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidB,4))==0.1433);
            pr(mixture, 500 /*psia*/, 100 + 460 /*R*/, "liquid");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(pr.selectedCompressibility(),4))==0.1768);
        }
        // 400 psia and 100 F
        {
            PhaseBehavior::VaporLiquidEquilibrium::rachfordVLE(mixture, 400 /*psia*/, 100 + 460 /*R*/);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture.molarFraction("vapor"),4))==0.8170);
            auto [vaporA,vaporB] = mixingRule(mixture, 400 /*psia*/, 100 + 460 /*R*/, "vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporA,4))==0.1320); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporB,4))==0.0345);
            pr(mixture, 400 /*psia*/, 100 + 460 /*R*/, "vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(pr.selectedCompressibility(),4))==0.9039);
            auto [liquidA,liquidB] = mixingRule(mixture, 400 /*psia*/, 100 + 460 /*R*/, "liquid");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidA,4))==1.5832); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidB,4))==0.1217);
            pr(mixture, 400 /*psia*/, 100 + 460 /*R*/, "liquid");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(pr.selectedCompressibility(),4))==0.1483);
        }
        // 300 psia and 100 F
        {
            PhaseBehavior::VaporLiquidEquilibrium::rachfordVLE(mixture, 300 /*psia*/, 100 + 460 /*R*/);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture.molarFraction("vapor"),4))==0.8400);
            auto [vaporA,vaporB] = mixingRule(mixture, 300 /*psia*/, 100 + 460 /*R*/, "vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporA,4))==0.1034); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporB,4))==0.0264);
            pr(mixture, 300 /*psia*/, 100 + 460 /*R*/, "vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(pr.selectedCompressibility(),4))==0.9234);
            auto [liquidA,liquidB] = mixingRule(mixture, 300 /*psia*/, 100 + 460 /*R*/, "liquid");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidA,4))==1.3559); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidB,4))==0.0980);
            pr(mixture, 300 /*psia*/, 100 + 460 /*R*/, "liquid");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(pr.selectedCompressibility(),4))==0.1178);
        }
    }

    SECTION("Using SoaveRedlichKwong"){

        using MixingRule = PhaseBehavior::EoS::MixingRules::RandomMixingRule<
                            PhaseBehavior::EoS::SRK::constants::omegaA, 
                            PhaseBehavior::EoS::SRK::constants::omegaB, 
                            PhaseBehavior::EoS::SRK::SoaveRedlichKwong
                            >;

        auto mixture = PhaseBehavior::Input::createMixtureFromFile("PVT.csv", "InteractionCoefficients.csv");

        MixingRule mixingRule;

        PhaseBehavior::EoS::SRK::SoaveRedlichKwong srk;

        // 500 psia and 50 F
        {
            PhaseBehavior::VaporLiquidEquilibrium::rachfordVLE(mixture, 500 /*psia*/, 50 + 460 /*R*/);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture.molarFraction("vapor"),4))==0.7252);

            auto [vaporA,vaporB] = mixingRule(mixture, 500 /*psia*/, 50 + 460 /*R*/, "vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporA,4))==0.1613); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporB,4))==0.0493);

            srk(mixture, 500 /*psia*/, 50 + 460 /*R*/, "vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(srk.selectedCompressibility(),4))==0.8865);


            auto [liquidA,liquidB] = mixingRule(mixture, 500 /*psia*/, 50 + 460 /*R*/, "liquid");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidA,4))==1.6592); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidB,4))==0.1507);

            srk(mixture, 500 /*psia*/, 50 + 460 /*R*/, "liquid");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(srk.selectedCompressibility(),4))==0.1875);
        }
        // 400 psia and 50 F
        {
            PhaseBehavior::VaporLiquidEquilibrium::rachfordVLE(mixture, 400 /*psia*/, 50 + 460 /*R*/);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture.molarFraction("vapor"),4))==0.7533);
            auto [vaporA,vaporB] = mixingRule(mixture, 400 /*psia*/, 50 + 460 /*R*/, "vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporA,4))==0.1345); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporB,4))==0.0401);
            srk(mixture, 400 /*psia*/, 50 + 460 /*R*/, "vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(srk.selectedCompressibility(),4))==0.9039);
            auto [liquidA,liquidB] = mixingRule(mixture, 400 /*psia*/, 50 + 460 /*R*/, "liquid");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidA,4))==1.4881); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidB,4))==0.1279);
            srk(mixture, 400 /*psia*/, 50 + 460 /*R*/, "liquid");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(srk.selectedCompressibility(),4))==0.1571);
        }
        // 300 psia and 50 F
        {
            PhaseBehavior::VaporLiquidEquilibrium::rachfordVLE(mixture, 300 /*psia*/, 50 + 460 /*R*/);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture.molarFraction("vapor"),4))==0.7840);
            auto [vaporA,vaporB] = mixingRule(mixture, 300 /*psia*/, 50 + 460 /*R*/, "vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporA,4))==0.1061); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporB,4))==0.0306);
            srk(mixture, 300 /*psia*/, 50 + 460 /*R*/, "vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(srk.selectedCompressibility(),4))==0.9231);
            auto [liquidA,liquidB] = mixingRule(mixture, 300 /*psia*/, 50 + 460 /*R*/, "liquid");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidA,4))==1.2814); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidB,4))==0.1031);
            srk(mixture, 300 /*psia*/, 50 + 460 /*R*/, "liquid");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(srk.selectedCompressibility(),4))==0.1248);
        }
        // 500 psia and 100 F
        {
            PhaseBehavior::VaporLiquidEquilibrium::rachfordVLE(mixture, 500 /*psia*/, 100 + 460 /*R*/);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture.molarFraction("vapor"),4))==0.7954);
            auto [vaporA,vaporB] = mixingRule(mixture, 500 /*psia*/, 100 + 460 /*R*/, "vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporA,4))==0.1432); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporB,4))==0.0473);
            srk(mixture, 500 /*psia*/, 100 + 460 /*R*/, "vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(srk.selectedCompressibility(),4))==0.9047);
            auto [liquidA,liquidB] = mixingRule(mixture, 500 /*psia*/, 100 + 460 /*R*/, "liquid");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidA,4))==1.7338); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidB,4))==0.1596);
            srk(mixture, 500 /*psia*/, 100 + 460 /*R*/, "liquid");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(srk.selectedCompressibility(),4))==0.1992);
        }
        // 400 psia and 100 F
        {
            PhaseBehavior::VaporLiquidEquilibrium::rachfordVLE(mixture, 400 /*psia*/, 100 + 460 /*R*/);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture.molarFraction("vapor"),4))==0.8170);
            auto [vaporA,vaporB] = mixingRule(mixture, 400 /*psia*/, 100 + 460 /*R*/, "vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporA,4))==0.1192); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporB,4))==0.0384);
            srk(mixture, 400 /*psia*/, 100 + 460 /*R*/, "vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(srk.selectedCompressibility(),4))==0.9192);
            auto [liquidA,liquidB] = mixingRule(mixture, 400 /*psia*/, 100 + 460 /*R*/, "liquid");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidA,4))==1.5537); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidB,4))==0.1355);
            srk(mixture, 400 /*psia*/, 100 + 460 /*R*/, "liquid");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(srk.selectedCompressibility(),4))==0.1670);
        }
        // 300 psia and 100 F
        {
            PhaseBehavior::VaporLiquidEquilibrium::rachfordVLE(mixture, 300 /*psia*/, 100 + 460 /*R*/);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture.molarFraction("vapor"),4))==0.8400);
            auto [vaporA,vaporB] = mixingRule(mixture, 300 /*psia*/, 100 + 460 /*R*/, "vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporA,4))==0.0936); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(vaporB,4))==0.0294);
            srk(mixture, 300 /*psia*/, 100 + 460 /*R*/, "vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(srk.selectedCompressibility(),4))==0.9354);
            auto [liquidA,liquidB] = mixingRule(mixture, 300 /*psia*/, 100 + 460 /*R*/, "liquid");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidA,4))==1.3337); CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(liquidB,4))==0.1091);
            srk(mixture, 300 /*psia*/, 100 + 460 /*R*/, "liquid");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(srk.selectedCompressibility(),4))==0.1326);
        }
    }
}