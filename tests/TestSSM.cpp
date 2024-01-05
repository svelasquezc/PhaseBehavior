#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <PhaseBehavior/Utilities/Input.hpp>
#include <PhaseBehavior/Utilities/Math.hpp>

#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/MixingRules.hpp>
#include <PhaseBehavior/EquationsOfState.hpp>
#include <PhaseBehavior/PhaseEquilibrium.hpp>


TEST_CASE("Can obtain the compressibility, A, and B, for a two phase mixture", "[RachfordVLE]"){

    SECTION("Using PR for 500 psia and 50 F"){

        auto mixture = PhaseBehavior::Input::createMixtureFromFile("PVT.csv", "InteractionCoefficients.csv");

        {
            PhaseBehavior::VaporLiquidEquilibrium::SuccesiveSubstitution<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 500 /*psia*/, 50 + 460 /*R*/);

            auto Z_liq = mixture.compressibility("liquid");
            auto Z_vap = mixture.compressibility("vapor");
            auto fng = mixture.molarFraction("vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7233);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.1636);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.8555);

            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["CO2"]->equilibriumCoefficient(),3))==2.049);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C1"]->equilibriumCoefficient(),3))==4.667);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C2"]->equilibriumCoefficient(),3))==0.947);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C3"]->equilibriumCoefficient(),3))==0.297);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["i-C4"]->equilibriumCoefficient(),3))==0.131);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C4"]->equilibriumCoefficient(),3))==0.094);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["i-C5"]->equilibriumCoefficient(),3))==0.041);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C5"]->equilibriumCoefficient(),3))==0.031);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C6"]->equilibriumCoefficient(),3))==0.011);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C7+"]->equilibriumCoefficient(),3))==0.000);
            
        }
    }

    SECTION("Using SRK for 500 psia and 50 F"){

        auto mixture = PhaseBehavior::Input::createMixtureFromFile("PVT.csv", "InteractionCoefficients.csv");
    
        PhaseBehavior::VaporLiquidEquilibrium::SuccesiveSubstitution<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 500 /*psia*/, 50 + 460 /*R*/);

        auto Z_liq = mixture.compressibility("liquid");
        auto Z_vap = mixture.compressibility("vapor");
        auto fng = mixture.molarFraction("vapor");
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7245);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.1847);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.8759);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["CO2"]->equilibriumCoefficient(),3))==1.973);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C1"]->equilibriumCoefficient(),3))==4.803);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C2"]->equilibriumCoefficient(),3))==0.948);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C3"]->equilibriumCoefficient(),3))==0.292);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["i-C4"]->equilibriumCoefficient(),3))==0.127);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C4"]->equilibriumCoefficient(),3))==0.091);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["i-C5"]->equilibriumCoefficient(),3))==0.039);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C5"]->equilibriumCoefficient(),3))==0.029);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C6"]->equilibriumCoefficient(),3))==0.010);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C7+"]->equilibriumCoefficient(),3))==0.000);
        
    }

    SECTION("Using PR for 400 psia and 50 F"){

        auto mixture = PhaseBehavior::Input::createMixtureFromFile("PVT.csv", "InteractionCoefficients.csv");

        {
            PhaseBehavior::VaporLiquidEquilibrium::SuccesiveSubstitution<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 400 /*psia*/, 50 + 460 /*R*/);

            auto Z_liq = mixture.compressibility("liquid");
            auto Z_vap = mixture.compressibility("vapor");
            auto fng = mixture.molarFraction("vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7521);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.1373);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.8805);

            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["CO2"]->equilibriumCoefficient(),3))==2.458);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C1"]->equilibriumCoefficient(),3))==5.759);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C2"]->equilibriumCoefficient(),3))==1.101);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C3"]->equilibriumCoefficient(),3))==0.329);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["i-C4"]->equilibriumCoefficient(),3))==0.140);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C4"]->equilibriumCoefficient(),3))==0.100);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["i-C5"]->equilibriumCoefficient(),3))==0.042);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C5"]->equilibriumCoefficient(),3))==0.032);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C6"]->equilibriumCoefficient(),3))==0.010);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C7+"]->equilibriumCoefficient(),3))==0.000);
            
        }
    }

    SECTION("Using SRK for 400 psia and 50 F"){

        auto mixture = PhaseBehavior::Input::createMixtureFromFile("PVT.csv", "InteractionCoefficients.csv");
    
        PhaseBehavior::VaporLiquidEquilibrium::SuccesiveSubstitution<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 400 /*psia*/, 50 + 460 /*R*/);

        auto Z_liq = mixture.compressibility("liquid");
        auto Z_vap = mixture.compressibility("vapor");
        auto fng = mixture.molarFraction("vapor");
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7530);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.1550);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.8973);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["CO2"]->equilibriumCoefficient(),3))==2.367);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C1"]->equilibriumCoefficient(),3))==5.937);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C2"]->equilibriumCoefficient(),3))==1.105);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C3"]->equilibriumCoefficient(),3))==0.326);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["i-C4"]->equilibriumCoefficient(),3))==0.137);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C4"]->equilibriumCoefficient(),3))==0.097);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["i-C5"]->equilibriumCoefficient(),3))==0.040);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C5"]->equilibriumCoefficient(),3))==0.030);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C6"]->equilibriumCoefficient(),3))==0.010);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C7+"]->equilibriumCoefficient(),3))==0.000);
        
    }

    SECTION("Using PR for 300 psia and 50 F"){

        auto mixture = PhaseBehavior::Input::createMixtureFromFile("PVT.csv", "InteractionCoefficients.csv");

        {
            PhaseBehavior::VaporLiquidEquilibrium::SuccesiveSubstitution<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 300 /*psia*/, 50 + 460 /*R*/);

            auto Z_liq = mixture.compressibility("liquid");
            auto Z_vap = mixture.compressibility("vapor");
            auto fng = mixture.molarFraction("vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7831);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.1093);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.9065);

            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["CO2"]->equilibriumCoefficient(),3))==3.149);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C1"]->equilibriumCoefficient(),3))==7.585);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C2"]->equilibriumCoefficient(),3))==1.367);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C3"]->equilibriumCoefficient(),3))==0.390);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["i-C4"]->equilibriumCoefficient(),3))==0.160);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C4"]->equilibriumCoefficient(),3))==0.113);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["i-C5"]->equilibriumCoefficient(),3))==0.045);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C5"]->equilibriumCoefficient(),3))==0.034);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C6"]->equilibriumCoefficient(),3))==0.011);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C7+"]->equilibriumCoefficient(),3))==0.000);
            
        }
    }

    SECTION("Using SRK for 300 psia and 50 F"){

        auto mixture = PhaseBehavior::Input::createMixtureFromFile("PVT.csv", "InteractionCoefficients.csv");
    
        PhaseBehavior::VaporLiquidEquilibrium::SuccesiveSubstitution<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 300 /*psia*/, 50 + 460 /*R*/);

        auto Z_liq = mixture.compressibility("liquid");
        auto Z_vap = mixture.compressibility("vapor");
        auto fng = mixture.molarFraction("vapor");
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7836);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.1233);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.9195);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["CO2"]->equilibriumCoefficient(),3))==3.031);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C1"]->equilibriumCoefficient(),3))==7.832);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C2"]->equilibriumCoefficient(),3))==1.376);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C3"]->equilibriumCoefficient(),3))==0.388);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["i-C4"]->equilibriumCoefficient(),3))==0.158);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C4"]->equilibriumCoefficient(),3))==0.110);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["i-C5"]->equilibriumCoefficient(),3))==0.044);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C5"]->equilibriumCoefficient(),3))==0.033);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C6"]->equilibriumCoefficient(),3))==0.010);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C7+"]->equilibriumCoefficient(),3))==0.000);
        
    }

    SECTION("Using PR for 300 psia and 100 F"){

        auto mixture = PhaseBehavior::Input::createMixtureFromFile("PVT.csv", "InteractionCoefficients.csv");

        {
            PhaseBehavior::VaporLiquidEquilibrium::SuccesiveSubstitution<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 300 /*psia*/, 100 + 460 /*R*/);

            auto Z_liq = mixture.compressibility("liquid");
            auto Z_vap = mixture.compressibility("vapor");
            auto fng = mixture.molarFraction("vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.8383);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.1155);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.9206);

            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["CO2"]->equilibriumCoefficient(),3))==4.430);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C1"]->equilibriumCoefficient(),3))==9.082);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C2"]->equilibriumCoefficient(),3))==2.060);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C3"]->equilibriumCoefficient(),3))==0.693);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["i-C4"]->equilibriumCoefficient(),3))==0.319);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C4"]->equilibriumCoefficient(),3))==0.235);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["i-C5"]->equilibriumCoefficient(),3))==0.106);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C5"]->equilibriumCoefficient(),3))==0.084);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C6"]->equilibriumCoefficient(),3))==0.031);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C7+"]->equilibriumCoefficient(),3))==0.001);
            
        }
    }

    SECTION("Using SRK for 300 psia and 100 F"){

        auto mixture = PhaseBehavior::Input::createMixtureFromFile("PVT.csv", "InteractionCoefficients.csv");
    
        PhaseBehavior::VaporLiquidEquilibrium::SuccesiveSubstitution<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 300 /*psia*/, 100 + 460 /*R*/);

        auto Z_liq = mixture.compressibility("liquid");
        auto Z_vap = mixture.compressibility("vapor");
        auto fng = mixture.molarFraction("vapor");
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.8384);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.1302);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.9330);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["CO2"]->equilibriumCoefficient(),3))==4.267);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C1"]->equilibriumCoefficient(),3))==9.311);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C2"]->equilibriumCoefficient(),3))==2.075);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C3"]->equilibriumCoefficient(),3))==0.693);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["i-C4"]->equilibriumCoefficient(),3))==0.317);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C4"]->equilibriumCoefficient(),3))==0.233);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["i-C5"]->equilibriumCoefficient(),3))==0.104);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C5"]->equilibriumCoefficient(),3))==0.081);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C6"]->equilibriumCoefficient(),3))==0.029);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C7+"]->equilibriumCoefficient(),3))==0.001);
        
    }

    SECTION("Using PR for 400 psia and 100 F"){

        auto mixture = PhaseBehavior::Input::createMixtureFromFile("PVT.csv", "InteractionCoefficients.csv");

        {
            PhaseBehavior::VaporLiquidEquilibrium::SuccesiveSubstitution<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 400 /*psia*/, 100 + 460 /*R*/);

            auto Z_liq = mixture.compressibility("liquid");
            auto Z_vap = mixture.compressibility("vapor");
            auto fng = mixture.molarFraction("vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.8155);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.1452);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.8985);

            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["CO2"]->equilibriumCoefficient(),3))==3.421);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C1"]->equilibriumCoefficient(),3))==6.864);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C2"]->equilibriumCoefficient(),3))==1.634);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C3"]->equilibriumCoefficient(),3))==0.571);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["i-C4"]->equilibriumCoefficient(),3))==0.270);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C4"]->equilibriumCoefficient(),3))==0.202);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["i-C5"]->equilibriumCoefficient(),3))==0.094);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C5"]->equilibriumCoefficient(),3))==0.074);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C6"]->equilibriumCoefficient(),3))==0.028);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C7+"]->equilibriumCoefficient(),3))==0.001);
            
        }
    }

    SECTION("Using SRK for 400 psia and 100 F"){

        auto mixture = PhaseBehavior::Input::createMixtureFromFile("PVT.csv", "InteractionCoefficients.csv");
    
        PhaseBehavior::VaporLiquidEquilibrium::SuccesiveSubstitution<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 400 /*psia*/, 100 + 460 /*R*/);

        auto Z_liq = mixture.compressibility("liquid");
        auto Z_vap = mixture.compressibility("vapor");
        auto fng = mixture.molarFraction("vapor");
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.8156);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.1638);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.9145);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["CO2"]->equilibriumCoefficient(),3))==3.296);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C1"]->equilibriumCoefficient(),3))==7.030);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C2"]->equilibriumCoefficient(),3))==1.642);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C3"]->equilibriumCoefficient(),3))==0.569);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["i-C4"]->equilibriumCoefficient(),3))==0.268);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C4"]->equilibriumCoefficient(),3))==0.198);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["i-C5"]->equilibriumCoefficient(),3))==0.091);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C5"]->equilibriumCoefficient(),3))==0.072);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C6"]->equilibriumCoefficient(),3))==0.027);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C7+"]->equilibriumCoefficient(),3))==0.001);
        
    }

    SECTION("Using PR for 500 psia and 100 F"){

        auto mixture = PhaseBehavior::Input::createMixtureFromFile("PVT.csv", "InteractionCoefficients.csv");

        {
            PhaseBehavior::VaporLiquidEquilibrium::SuccesiveSubstitution<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 500 /*psia*/, 100 + 460 /*R*/);

            auto Z_liq = mixture.compressibility("liquid");
            auto Z_vap = mixture.compressibility("vapor");
            auto fng = mixture.molarFraction("vapor");
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7941);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.1730);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.8777);

            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["CO2"]->equilibriumCoefficient(),3))==2.819);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C1"]->equilibriumCoefficient(),3))==5.538);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C2"]->equilibriumCoefficient(),3))==1.383);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C3"]->equilibriumCoefficient(),3))==0.502);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["i-C4"]->equilibriumCoefficient(),3))==0.245);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C4"]->equilibriumCoefficient(),3))==0.184);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["i-C5"]->equilibriumCoefficient(),3))==0.088);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C5"]->equilibriumCoefficient(),3))==0.070);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C6"]->equilibriumCoefficient(),3))==0.028);
            CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C7+"]->equilibriumCoefficient(),3))==0.001);
            
        }
    }

    SECTION("Using SRK for 500 psia and 100 F"){

        auto mixture = PhaseBehavior::Input::createMixtureFromFile("PVT.csv", "InteractionCoefficients.csv");
    
        PhaseBehavior::VaporLiquidEquilibrium::SuccesiveSubstitution<PhaseBehavior::EoS::SRK::SoaveRedlichKwong>(mixture, 500 /*psia*/, 100 + 460 /*R*/);

        auto Z_liq = mixture.compressibility("liquid");
        auto Z_vap = mixture.compressibility("vapor");
        auto fng = mixture.molarFraction("vapor");
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(fng,4))==0.7942);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_liq,4))==0.1951);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(Z_vap,4))==0.8969);

        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["CO2"]->equilibriumCoefficient(),3))==2.718);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C1"]->equilibriumCoefficient(),3))==5.666);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C2"]->equilibriumCoefficient(),3))==1.388);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C3"]->equilibriumCoefficient(),3))==0.499);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["i-C4"]->equilibriumCoefficient(),3))==0.241);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C4"]->equilibriumCoefficient(),3))==0.180);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["i-C5"]->equilibriumCoefficient(),3))==0.085);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C5"]->equilibriumCoefficient(),3))==0.068);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["n-C6"]->equilibriumCoefficient(),3))==0.026);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(mixture["C7+"]->equilibriumCoefficient(),3))==0.001);
        
    }
}