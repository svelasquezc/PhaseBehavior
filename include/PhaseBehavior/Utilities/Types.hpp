#ifndef TYPES_HPP
#define TYPES_HPP

#include <queue>
#include <deque>

#include "Math.hpp"

// Namespace defining custom types related to numerical precision
namespace PhaseBehavior::Types {
    /**
     * @brief Alias for the numerical precision type.
     *
     * This type represents the level of numerical precision used
     * throughout the codebase.
     */
    using NumericalPrecision = double;
    using ConstexprConstant = decltype(Math::constants::pi<NumericalPrecision>);

    template <typename T, int MaxLen, typename Container=std::deque<T>>
    class FixedQueue : public std::queue<T, Container> {
    public:
        void push(const T& value) {
            if (this->size() == MaxLen) {
            this->c.pop_front();
            }
            std::queue<T, Container>::push(value);
        }
    };
}

#endif /* TYPES_HPP */
