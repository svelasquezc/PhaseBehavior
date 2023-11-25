#ifndef TYPES_HPP
#define TYPES_HPP

#include "Math.hpp"

// Namespace defining custom types related to numerical precision
namespace PhaseBehavior::Types {
    /**
     * @brief Alias for the numerical precision type.
     *
     * This type represents the level of numerical precision used
     * throughout the codebase.
     */
    using NumericalPrecision = long double;
    using ConstexprConstant = decltype(Math::constants::pi<NumericalPrecision>);   
}

#endif /* TYPES_HPP */
