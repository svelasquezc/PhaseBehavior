#ifndef NUMERICS_HPP
#define NUMERICS_HPP

#include <type_traits>
#include <functional>
#include <cassert>

#include <Eigen/Dense>

#include "Types.hpp"

using NP_t = PhaseBehavior::Types::NumericalPrecision;

namespace PhaseBehavior::Numerics{
    
    class NewtonRaphson{
        private:
            Eigen::MatrixXd jacobian_;
            Eigen::VectorXd residualVector_, solutionVector_;
            std::size_t size_ = 0;
        public:
            template<template<typename> typename ContainerType, typename FunctionLike>
            std::optional<ContainerType<NP_t>> operator()(ContainerType<NP_t> const& initialGuess, ContainerType<FunctionLike> const& residualSystem){
                assert(initialGuess.size()==residualSystem.size());
                if (initialGuess.size() != size_) {
                    size_ = initialGuess.size();
                    jacobian_.resize(size_, size_);
                    residualVector_.resize(size_);
                    solutionVector_.resize(size_);
                }
                jacobian_.setZero();
                residualVector_.setZero();
                solutionVector_.setZero();

                Eigen::VectorXd newGuess = initialGuess;
                Eigen::VectorXd oldGuess;
                oldValue.setZero();

                // Relative error convergence criterium
                for (std::size_t i = 0; i<500; ++i){
                    // Centered first-order derivative
                    PT goalDerivativeValue = (goalFunction(newValue + scaledEpsilon) - goalFunction(newValue - scaledEpsilon))/(2*scaledEpsilon);
                    oldValue = newValue;
                    newValue = oldValue - (goalFunction(oldValue)/goalDerivativeValue);

                    if ((newGuess - oldGuess).squaredNorm() < 1e-10) return std::make_optional(newValue); 
                }
                return std::optional<PT>();
            }

            return ContainerType<> operator()
    };
}

#endif /* NUMERICS_HPP */