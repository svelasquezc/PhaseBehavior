#ifndef NUMERICS_HPP
#define NUMERICS_HPP

#include <Eigen/Dense>



namespace PhaseBehavior::Numerics{
    template<typename ContainerType>
    class NewtonRaphson{
        private:
            Eigen::MatrixXd jacobian_;
            Eigen::VectorXd residualVector_, solutionVector_;
        public:
            NewtonRaphson(ContainerType initialGuess ){

            }
    };
}

#endif /* NUMERICS_HPP */