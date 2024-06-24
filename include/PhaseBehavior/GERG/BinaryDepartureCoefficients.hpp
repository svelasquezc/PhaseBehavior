#ifndef BINARYDEPARTURECOEFFICIENTS_HPP
#define BINARYDEPARTURECOEFFICIENTS_HPP

#include <tuple>
#include <vector>
#include <stdexcept>

#include "../Utilities/Types.hpp"

#include "Types.hpp"

using NP_t = PhaseBehavior::Types::NumericalPrecision;

namespace PhaseBehavior::EoS::GERG::Coefficients::Departure {

    static constexpr unsigned int reducedDensityExponent = 0;  
    static constexpr unsigned int reducedTemperatureExponent = 1;
    static constexpr unsigned int exponentialReducedDensityExponent = 2;
    static constexpr unsigned int eta = 2;
    static constexpr unsigned int varepsilon = 3;
    static constexpr unsigned int beta = 4;
    static constexpr unsigned int gamma = 5;

    using BinaryFParameter = GERG::Types::BinaryParameter<NP_t>;

    BinaryFParameter const f = {
        {{"CH4", "N2"}, 1.0}, {{"CH4", "CO2"}, 1.0}, {{"CH4", "C2H6"}, 1.0}, {{"CH4", "C3H8"}, 1.0}, {{"CH4", "n-C4H10"}, 1.0}, {{"CH4", "i-C4H10"}, 0.771035405688}, {{"CH4", "H2"}, 1.0},
        {{"N2", "CO2"}, 1.0}, {{"N2", "C2H6"}, 1.0},
        {{"C2H6", "C3H8"}, 0.130424765150}, {{"C2H6", "n-C4H10"}, 0.281570073085}, {{"C2H6", "i-C4H10"}, 0.260632376098},
        {{"C3H8", "n-C4H10"}, 0.312572600489e-1}, {{"C3H8", "i-C4H10"}, -0.551609771024e-1}, 
        {{"n-C4H10", "i-C4H10"}, -0.551240293009e-1}
        };

    const inline static NP_t fParameter(std::string_view component1, std::string_view component2) {
        try {
            return f.at({component1, component2});
        }catch(std::out_of_range e){
            return 0;
        }
    }

    struct DepartureFunction {
        using DepartureCoefficient = std::vector<NP_t>;
        using PolynomialExponent = std::vector<std::pair<int, NP_t>>;
        using ExponentialExponent = std::vector<std::tuple<int, NP_t, NP_t, NP_t, NP_t, NP_t>>;

        DepartureCoefficient const coefficient;
        PolynomialExponent const polynomialExponent;
        ExponentialExponent const exponentialExponent;
    };

    using BinaryDeparture = GERG::Types::BinaryParameter<DepartureFunction const>;

    static const BinaryDeparture binaryDeparture = {
        {
            {"CH4", "N2"}, {
               /* Departure coefficients n_ij */ 
               {-0.98038985517335e-2, 0.42487270143005e-3, -0.34800214576142e-1, -0.13333813013896, -0.11993694974627e-1, 0.69243379775168e-1, -0.31022508148249, 0.24495491753226, 0.22369816716981},
               /* Polynomial exponents d_ij, t_ij */
               {{1, 0.000}, {4, 1.850}},
               /* Exponential exponents d_ij, t_ij, eta_ij, epsilon_ij, beta_ij, gamma_ij*/
               {{1, 7.850, 1.000, 0.5, 1.000, 0.5}, {2, 5.400, 1.000, 0.5, 1.000, 0.5}, {2, 0.000, 0.250, 0.5, 0.250, 0.5}, {2, 0.750, 0.000, 0.5, 3.000, 0.5}, {2, 2.800, 0.000, 0.5, 3.000, 0.5}, {2, 4.450, 0.000, 0.5, 3.000, 0.5}, {3, 4.250, 0.000, 0.5, 3.000, 0.5} }
            },
        },
        {
            {"CH4", "CO2"}, {
               /* Departure coefficients n_ij */ 
               {-0.10859387354942, 0.80228576727389E-1, -0.93303985115717E-2, 0.40989274005848E-1, -0.24338019772494, 0.23855347281124},
               /* Polynomial exponents d_ij, t_ij */
               {{1, 2.600}, {2, 1.950}, {3, 0.000}},
               /* Exponential exponents d_ij, t_ij, eta_ij, epsilon_ij, beta_ij, gamma_ij*/
               {{1, 3.950, 1.000, 0.5, 1.000, 0.5}, {2, 7.950, 0.500, 0.5, 2.000, 0.5}, {3, 8.000, 0.000, 0.5, 3.000, 0.5}}
            },
        },
        {
            {"CH4", "C2H6"}, {
               /* Departure coefficients n_ij */ 
               {-0.80926050298746E-3, -0.75381925080059E-3, -0.41618768891219E-1, -0.23452173681569, 0.14003840584586, 0.63281744807738E-1, -0.34660425848809E-1, -0.23918747334251, 0.19855255066891E-2, 0.61777746171555E1, -0.69575358271105E1, 0.10630185306388E1},
               /* Polynomial exponents d_ij, t_ij */
               {{3, 0.650}, {4, 1.550}},
               /* Exponential exponents d_ij, t_ij, eta_ij, epsilon_ij, beta_ij, gamma_ij*/
               {{1, 3.100, 1.000, 0.5, 1.000, 0.5}, {2, 5.900, 1.000, 0.5, 1.000, 0.5}, {2, 7.050, 1.000, 0.5, 1.000, 0.5}, {2, 3.350, 0.875, 0.5, 1.250, 0.5}, {2, 1.200, 0.750, 0.5, 1.500, 0.5}, {2, 5.800, 0.500, 0.5, 2.000, 0.5}, {2, 2.700, 0.000, 0.5, 3.000, 0.5}, {3, 0.450, 0.000, 0.5, 3.000, 0.5}, {3, 0.550, 0.000, 0.5, 3.000, 0.5}, {3, 1.950, 0.000, 0.5, 3.000, 0.5}}
            },
        },
        {
            {"CH4", "C3H8"}, {
               /* Departure coefficients n_ij */ 
               {0.13746429958576E-1, -0.74425012129552E-2, -0.45516600213685E-2, -0.54546603350237E-2, 0.23682016824471E-2, 0.18007763721438, -0.44773942932486, 0.19327374888200E-1, -0.30632197804624},
               /* Polynomial exponents d_ij, t_ij */
               {{3, 1.850}, {3, 3.950}, {4, 0.000}, {4, 1.850}, {4, 3.850}},
               /* Exponential exponents d_ij, t_ij, eta_ij, epsilon_ij, beta_ij, gamma_ij*/
               {{1, 5.250, 0.250, 0.5, 0.750, 0.5}, {1, 3.850, 0.250, 0.5, 1.000, 0.5}, {1, 0.200, 0.000, 0.5, 2.000, 0.5}, {2, 6.500, 0.000, 0.5, 3.000, 0.5}}
            },
        },
        {
            {"N2", "CO2"}, {
               /* Departure coefficients n_ij */ 
               {0.28661625028399, -0.10919833861247, -0.11374032082270E1, 0.76580544237358, 0.42638000926819E-2, 0.17673538204534},
               /* Polynomial exponents d_ij, t_ij */
               {{2, 1.850}, {3, 1.400}},
               /* Exponential exponents d_ij, t_ij, eta_ij, epsilon_ij, beta_ij, gamma_ij*/
               {{1, 3.200, 0.250, 0.5, 0.750, 0.5}, {1, 2.500, 0.250, 0.5, 1.000, 0.5}, {1, 8.000, 0.000, 0.5, 2.000, 0.5}, {2, 3.750, 0.000, 0.5, 3.000, 0.5}}
            },
        },
        {
            {"N2", "C2H6"}, {
               /* Departure coefficients n_ij */ 
               {-0.47376518126608, 0.48961193461001, -0.57011062090535E-2, -0.19966820041320, -0.69411103101723, 0.69226192739021},
               /* Polynomial exponents d_ij, t_ij */
               {{2, 0.000}, {2, 0.050}, {3, 0.000}},
               /* Exponential exponents d_ij, t_ij, eta_ij, epsilon_ij, beta_ij, gamma_ij*/
               {{1, 3.650, 1.000, 0.5, 1.000, 0.5}, {2, 4.900, 1.000, 0.5, 1.000, 0.5}, {2, 4.450, 0.875, 0.5, 1.250, 0.5}}
            },
        },
        {
            {"CH4", "H2"}, {
               /* Departure coefficients n_ij */ 
               {-0.25157134971934, -0.62203841111983E-2, 0.88850315184396E-1, -0.35592212573239E-1},
               /* Polynomial exponents d_ij, t_ij */
               {{1, 2.000}, {3, -1.000}, {3, 1.750}, {4, 1.400}},
               /* Exponential exponents d_ij, t_ij, eta_ij, epsilon_ij, beta_ij, gamma_ij*/
               {}
            },
        },
        {
            {"CH4", "n-C4H10"}, {
               /* Departure coefficients n_ij */ 
               {0.25574776844118E1, -0.79846357136353E1, 0.47859131465806E1, -0.73265392369587, 0.13805471345312E1, 0.28349603476365, -0.49087385940425, -0.10291888921447, 0.11836314681968, 0.55527385721943e-4},
               /* Polynomial exponents d_ij, t_ij */
               { {1,1.000}, {1,1.550}, {1, 1.700}, {2, 0.250}, {2, 1.350}, {3, 0.000}, {3, 1.250}, {4, 0.000}, {4, 0.700}, {4, 5.400} },
               /* Exponential exponents d_ij, t_ij, eta_ij, epsilon_ij, beta_ij, gamma_ij*/
               {}
            },
        },
        {
            {"CH4", "i-C4H10"}, {
               /* Departure coefficients n_ij */ 
               {0.25574776844118E1, -0.79846357136353E1, 0.47859131465806E1, -0.73265392369587, 0.13805471345312E1, 0.28349603476365, -0.49087385940425, -0.10291888921447, 0.11836314681968, 0.55527385721943e-4},
               /* Polynomial exponents d_ij, t_ij */
               { {1,1.000}, {1,1.550}, {1, 1.700}, {2, 0.250}, {2, 1.350}, {3, 0.000}, {3, 1.250}, {4, 0.000}, {4, 0.700}, {4, 5.400} },
               /* Exponential exponents d_ij, t_ij, eta_ij, epsilon_ij, beta_ij, gamma_ij*/
               {}
            },
        },
        {
            {"C2H6", "C3H8"}, {
               /* Departure coefficients n_ij */ 
               {0.25574776844118E1, -0.79846357136353E1, 0.47859131465806E1, -0.73265392369587, 0.13805471345312E1, 0.28349603476365, -0.49087385940425, -0.10291888921447, 0.11836314681968, 0.55527385721943e-4},
               /* Polynomial exponents d_ij, t_ij */
               { {1,1.000}, {1,1.550}, {1, 1.700}, {2, 0.250}, {2, 1.350}, {3, 0.000}, {3, 1.250}, {4, 0.000}, {4, 0.700}, {4, 5.400} },
               /* Exponential exponents d_ij, t_ij, eta_ij, epsilon_ij, beta_ij, gamma_ij*/
               {}
            },
        },
        {
            {"C2H6", "n-C4H10"}, {
               /* Departure coefficients n_ij */ 
               {0.25574776844118E1, -0.79846357136353E1, 0.47859131465806E1, -0.73265392369587, 0.13805471345312E1, 0.28349603476365, -0.49087385940425, -0.10291888921447, 0.11836314681968, 0.55527385721943e-4},
               /* Polynomial exponents d_ij, t_ij */
               { {1,1.000}, {1,1.550}, {1, 1.700}, {2, 0.250}, {2, 1.350}, {3, 0.000}, {3, 1.250}, {4, 0.000}, {4, 0.700}, {4, 5.400} },
               /* Exponential exponents d_ij, t_ij, eta_ij, epsilon_ij, beta_ij, gamma_ij*/
               {}
            },
        },
        {
            {"C2H6", "i-C4H10"}, {
               /* Departure coefficients n_ij */ 
               {0.25574776844118E1, -0.79846357136353E1, 0.47859131465806E1, -0.73265392369587, 0.13805471345312E1, 0.28349603476365, -0.49087385940425, -0.10291888921447, 0.11836314681968, 0.55527385721943e-4},
               /* Polynomial exponents d_ij, t_ij */
               { {1,1.000}, {1,1.550}, {1, 1.700}, {2, 0.250}, {2, 1.350}, {3, 0.000}, {3, 1.250}, {4, 0.000}, {4, 0.700}, {4, 5.400} },
               /* Exponential exponents d_ij, t_ij, eta_ij, epsilon_ij, beta_ij, gamma_ij*/
               {}
            },
        },
        {
            {"C3H8", "n-C4H10"}, {
               /* Departure coefficients n_ij */ 
               {0.25574776844118E1, -0.79846357136353E1, 0.47859131465806E1, -0.73265392369587, 0.13805471345312E1, 0.28349603476365, -0.49087385940425, -0.10291888921447, 0.11836314681968, 0.55527385721943e-4},
               /* Polynomial exponents d_ij, t_ij */
               { {1,1.000}, {1,1.550}, {1, 1.700}, {2, 0.250}, {2, 1.350}, {3, 0.000}, {3, 1.250}, {4, 0.000}, {4, 0.700}, {4, 5.400} },
               /* Exponential exponents d_ij, t_ij, eta_ij, epsilon_ij, beta_ij, gamma_ij*/
               {}
            },
        },
        {
            {"C3H8", "i-C4H10"}, {
               /* Departure coefficients n_ij */ 
               {0.25574776844118E1, -0.79846357136353E1, 0.47859131465806E1, -0.73265392369587, 0.13805471345312E1, 0.28349603476365, -0.49087385940425, -0.10291888921447, 0.11836314681968, 0.55527385721943e-4},
               /* Polynomial exponents d_ij, t_ij */
               { {1,1.000}, {1,1.550}, {1, 1.700}, {2, 0.250}, {2, 1.350}, {3, 0.000}, {3, 1.250}, {4, 0.000}, {4, 0.700}, {4, 5.400} },
               /* Exponential exponents d_ij, t_ij, eta_ij, epsilon_ij, beta_ij, gamma_ij*/
               {}
            },
        },
        {
            {"n-C4H10", "i-C4H10"}, {
               /* Departure coefficients n_ij */ 
               {0.25574776844118E1, -0.79846357136353E1, 0.47859131465806E1, -0.73265392369587, 0.13805471345312E1, 0.28349603476365, -0.49087385940425, -0.10291888921447, 0.11836314681968, 0.55527385721943e-4},
               /* Polynomial exponents d_ij, t_ij */
               { {1,1.000}, {1,1.550}, {1, 1.700}, {2, 0.250}, {2, 1.350}, {3, 0.000}, {3, 1.250}, {4, 0.000}, {4, 0.700}, {4, 5.400} },
               /* Exponential exponents d_ij, t_ij, eta_ij, epsilon_ij, beta_ij, gamma_ij*/
               {}
            },
        },
    };
}

#endif /* BINARYDEPARTURECOEFFICIENTS_HPP */