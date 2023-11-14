#ifndef COMPONENT_HPP
#define COMPONENT_HPP

#include <string>
#include "Utilities/Types.hpp"

/**
 * @brief Class representing a chemical component with various properties.
 */
class Component {
private:
    using NP_t = TypesDefinition::NumericalPrecision;

    std::string mutable name_;      /**< Name of the chemical component. */
    NP_t criticalPressure_;         /**< Critical pressure property. */
    NP_t criticalTemperature_;      /**< Critical temperature property. */
    NP_t criticalVolume_;           /**< Critical volume property. */
    NP_t molarWeight_;              /**< Molar weight property. */
    NP_t accentricFactor_;          /**< Accentric factor property. */

public:
    /**
     * @brief Constructor to initialize the component properties.
     * @param name Name of the chemical component.
     * @param criticalPressure Critical pressure value.
     * @param criticalTemperature Critical temperature value.
     * @param criticalVolume Critical volume value.
     * @param molarWeight Molar weight value.
     * @param accentricFactor Accentric factor value.
     */
    Component(std::string name, NP_t criticalPressure, NP_t criticalTemperature,
              NP_t criticalVolume, NP_t molarWeight, NP_t accentricFactor)
        : name_(std::move(name)),
          criticalPressure_(criticalPressure),
          criticalVolume_(criticalVolume),
          criticalTemperature_(criticalTemperature),
          molarWeight_(molarWeight),
          accentricFactor_(accentricFactor)
    {}

    /// Default copy constructor
    Component(Component& rhs) = default;
    /// Default move constructor
    Component(Component&& rhs) = default;
    
    /// Default copy assignment
    Component& operator=(Component& rhs) = default;
    /// Default move assignment
    Component& operator=(Component&& rhs) = default;

    /**
     * @brief Getter for the name of the component.
     * @return Name of the chemical component.
     */
    const std::string& name() const {
        return name_;
    }

    /**
     * @brief Getter for the critical pressure of the component.
     * @return Critical pressure property.
     */
    NP_t criticalPressure() const {
        return criticalPressure_;
    }

    /**
     * @brief Getter for critical temperature of the component.
     * @return Critical temperature property.
     */
    NP_t criticalTemperature() const {
        return criticalTemperature_;
    }

    /**
     * @brief Getter for critical volume of the component.
     * @return Critical volume property.
     */
    NP_t criticalVolume() const {
        return criticalVolume_;
    }

    /**
     * @brief Getter for molar weight of the component.
     * @return Molar weight property.
     */
    NP_t molarWeight() const {
        return molarWeight_;
    }

    /**
     * @brief Getter for accentric factor of the component.
     * @return Accentric factor property.
     */
    NP_t accentricFactor() const {
        return accentricFactor_;
    }

    /**
     * @brief Calculates the reduced pressure of the component.
     * @param absolutePressure The absolute pressure.
     * @return Reduced pressure.
     */
    NP_t reducedPressure(NP_t const& absolutePressure) const {
        return absolutePressure / criticalPressure_;
    }

    /**
     * @brief Calculates the reduced temperature of the component.
     * @param absoluteTemperature The absolute temperature.
     * @return Reduced temperature.
     */
    NP_t reducedTemperature(NP_t const& absoluteTemperature) const {
        return absoluteTemperature / criticalTemperature_;
    }

    /// Default destructor
    ~Component() = default;
};

#endif /* COMPONENT_HPP */
