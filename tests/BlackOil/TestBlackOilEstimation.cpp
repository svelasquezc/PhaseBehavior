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

    auto mixture = PhaseBehavior::Input::createMixtureFromFile<true>("PVTPR.csv", "InteractionCoefficients.csv");

    auto BOProperties = PhaseBehavior::BlackOil::BlackOilPropertiesEstimation<PhaseBehavior::EoS::PR::PengRobinson>(mixture, 3500, 250 + 460, 300, 62 + 460, 20, 60 + 460, 14.7, 60 + 460, 100);

    auto bomTable = BOProperties.PVTTable();

    std::vector<PVTTableTest> bomTestData;
    {
        auto bomFile = std::ifstream("BOMHandoutData.csv");
        std::string line;
        // Get rid of the header line
        std::getline(bomFile, line);
        while (std::getline(bomFile, line)){
            std::istringstream ss(line);
            PVTTableTest row;
            ss >> row.P >> row.Bo >> row.Bg >> row.Rs >> row.Rv;
            if(row.P <= 2750)
            bomTestData.push_back(row);
        }
        bomFile.close();
    }

    for (std::size_t i = 0; i < bomTable.size(); ++i){
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(bomTable[i].pressure,1))==bomTestData[i].P);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(bomTable[i].oilVolumetricFactor,3))==bomTestData[i].Bo);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(bomTable[i].gasVolumetricFactor,3))==bomTestData[i].Bg);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(bomTable[i].dissolvedGas,3))==bomTestData[i].Rs);
        CHECK(Catch::Approx(PhaseBehavior::Math::roundUp(bomTable[i].volatilizedOil,3))==bomTestData[i].Rv);
    }

}