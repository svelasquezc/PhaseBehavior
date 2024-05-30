#include <iostream>
#include <fstream>
#include <chrono>

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


int main(){
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
}