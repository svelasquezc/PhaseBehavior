#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>

#include <PhaseBehavior/Utilities/Types.hpp>
#include <PhaseBehavior/Mixture.hpp>
#include <PhaseBehavior/Component.hpp>

using Precision_t = TypesDefinition::NumericalPrecision;

TEST_CASE("Can create Mixture objects", "[mixture]"){

    SECTION("Mixture object from file"){
        auto pvtFile = std::ifstream("PVT.csv");
        std::vector<std::pair<Component, Precision_t>> components;

        std::string line;

        // Get rid of the header line
        std::getline(pvtFile, line);
        while (std::getline(pvtFile, line)){
            std::istringstream ss(line);
            std::string name;
            Precision_t composition, Pc, Tc, w, MW, Vc;
            ss >> name >> composition >> Pc >> Tc >> w >> MW >> Vc;
            components.push_back({Component(std::move(name), Pc, Tc, Vc, MW, w ), composition});
        }

        Mixture mixture(components);

        REQUIRE(true);
    }

    SECTION("Mixture object from known components"){
        Component CO2 {"CO2",1071,547.91,0.2667,4401,0.0344};
        Component C1 {"CO2",1071,547.91,0.2667,4401,0.0344};

        Mixture((CO2, 0.0031), (C1, 0.6192));
    }

}