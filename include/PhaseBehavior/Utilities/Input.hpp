#ifndef INPUT_HPP
#define INPUT_HPP

#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <exception>

#include "../Mixture.hpp"


namespace PhaseBehavior::Input {
    using Precision_t = PhaseBehavior::Types::NumericalPrecision;
    using MixComp = std::pair<PhaseBehavior::Component, Precision_t>;

    namespace __utils{
        char detectSeparator(const std::string &line) {
            std::vector<char> separators = {',', ';', '|', '\t', ' '};

            for (char sep : separators) {
                if (line.find(sep) != std::string::npos) {
                    return sep;
                }
            }
            return ',';
        }
    };

    template<bool hasShift=false>
    Mixture createMixtureFromFile(std::string const& pvtFilename, std::string const& interactionsFilename){
        using __utils::detectSeparator;

        std::vector<MixComp> components;

        {
            auto pvtFile = std::ifstream(pvtFilename);
            std::string line;
            // Get rid of the header line
            std::getline(pvtFile, line);

            char separator = detectSeparator(line);

            while (std::getline(pvtFile, line)){
                std::istringstream ss(line);
                std::string name, value;
                Precision_t composition, Pc, Tc, w, MW, Vc, shift;
                // Parse the line using the detected separator
                std::getline(ss, name, separator); // Component name
                std::getline(ss, value, separator); composition = std::stod(value);  // Mole fraction
                std::getline(ss, value, separator); Pc = std::stod(value);  // Critical temperature
                std::getline(ss, value, separator); Tc = std::stod(value);  // Critical pressure
                std::getline(ss, value, separator); w = std::stod(value);  // Acentric factor
                std::getline(ss, value, separator); MW = std::stod(value);  // Molecular weight
                std::getline(ss, value, separator); Vc = std::stod(value);  // Critical Volume
                if constexpr (hasShift){
                    std::getline(ss, value, separator); shift = std::stod(value);  // Volume shift
                    components.push_back({PhaseBehavior::Component(std::move(name), Pc, Tc, Vc, MW, w, shift ), composition});
                }else{
                    components.push_back({PhaseBehavior::Component(std::move(name), Pc, Tc, Vc, MW, w ), composition});
                }
            }
            pvtFile.close();
        }
        PhaseBehavior::Mixture mixture(components);
        {
            auto interactionCoefficientsFile = std::ifstream(interactionsFilename);
            std::string line;

            std::getline(interactionCoefficientsFile, line);

            char separator = detectSeparator(line);

            for (const auto& mixComponent1 : mixture){
                std::getline(interactionCoefficientsFile, line);
                std::istringstream ss(line);
                std::string name, value;
                std::getline(ss, name, separator);
                for (const auto& mixComponent2 : mixture){
                    std::getline(ss, value, separator);
                    Precision_t interactionCoefficientValue = 0.0;
                    if (!value.empty()) {
                        interactionCoefficientValue = std::stod(value);
                        if (interactionCoefficientValue != 0){
                            mixture.interactionCoefficient(mixComponent1, mixComponent2, interactionCoefficientValue);
                            mixture.interactionCoefficient(mixComponent2, mixComponent1, interactionCoefficientValue);
                        }
                    }
                }
            }
        }
        return mixture;
    };

}
#endif /* INPUT_HPP */