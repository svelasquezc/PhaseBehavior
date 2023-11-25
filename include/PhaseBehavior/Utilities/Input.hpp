#ifndef INPUT_HPP
#define INPUT_HPP

#include <string>
#include <fstream>
#include <sstream>
#include <vector>

#include "../Mixture.hpp"


namespace PhaseBehavior::Input {
    using Precision_t = PhaseBehavior::Types::NumericalPrecision;
    using MixComp = std::pair<PhaseBehavior::Component, Precision_t>;

    Mixture createMixtureFromFile(std::string pvtFile, std::string interactionsFile){
        std::vector<MixComp> components;

        {
            auto pvtFile = std::ifstream("PVT.csv");
            std::string line;

            // Get rid of the header line
            std::getline(pvtFile, line);
            while (std::getline(pvtFile, line)){
                std::istringstream ss(line);
                std::string name;
                Precision_t composition, Pc, Tc, w, MW, Vc;
                ss >> name >> composition >> Pc >> Tc >> w >> MW >> Vc;
                components.push_back({PhaseBehavior::Component(std::move(name), Pc, Tc, Vc, MW, w ), composition});
            }
            pvtFile.close();
        }
        PhaseBehavior::Mixture mixture(components);
        {
            auto interactionCoefficientsFile = std::ifstream("InteractionCoefficients.csv");
            std::string line;

            std::getline(interactionCoefficientsFile, line);

            for (const auto& mixComponent1 : mixture){
                std::getline(interactionCoefficientsFile, line);
                std::istringstream ss(line);
                std::string name;
                ss >> name;
                for (const auto& mixComponent2 : mixture){
                    Precision_t interactionCoefficientValue = 0.0;
                    ss >> interactionCoefficientValue;
                    if (std::abs(interactionCoefficientValue) > 0.0){
                        mixture.interactionCoefficient(mixComponent1, mixComponent2, interactionCoefficientValue);
                    }
                }
            }
        }
        return mixture;
    }

}
#endif /* INPUT_HPP */