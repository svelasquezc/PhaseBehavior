#ifndef PUREFLUIDHELMHOLTZCOEFFICIENTS_HPP
#define PUREFLUIDHELMHOLTZCOEFFICIENTS_HPP

#include <map>
#include <tuple>
#include <string_view>
#include "../Utilities/Types.hpp"

using NP_t = PhaseBehavior::Types::NumericalPrecision;

namespace PhaseBehavior::EoS::GERG::Coefficients::PureFluid {

    namespace Ideal {

        struct PureFluid {
            std::array<NP_t, 7> const polynomialCoefficient;
            std::array<NP_t, 4> const hyperbolicCoefficient;
        };

        const static std::map<std::string_view, PureFluid const> pureFluid = {
            {
                "CH4", { // Methane
                        {19.597508817, -83.959667892, 3.00088, 0.76315, 0.00460, 8.74432, -4.46921},
                        {4.306474465, 0.936220902, 5.577233895, 5.722644361}
                }
            }, 
            {
                "N2", { // Nitrogen
                        {11.083407489, -22.202102428, 2.50031, 0.13732, -0.14660, 0.90066, 0},
                        {5.251822620, -5.393067706, 13.788988208, 0}
                }
            }, 
            {
                "CO2", {
                    {11.925152758, -16.118762264, 2.50002, 2.04452, -1.06044, 2.03366, 0.01393},
                    {3.022758166, -2.844425476, 1.589964364, 1.121596090}
                }
            }, // Carbon Dioxide
            {
                "C2H6", { // Ethane
                    {24.675437527, -77.425313760, 3.00263, 4.33939, 1.23722, 13.19740, -6.01989},
                    {1.831882406, 0.731306621, 3.378007481, 3.508721939}
                }
            }, 
            {"C3H8", {{31.602908195, -84.463284382, 3.02939, 6.60569, 3.19700, 19.19210, -8.37267}, {1.297521801, 0.543210978, 2.583146083, 2.777773271}}}, // Propane
            {"n-C4H10", {{20.884143364, -91.638478026, 3.33944, 9.44893, 6.89406, 24.46180, 14.78240}, {1.101487798, 0.431957660, 4.502440459, 2.124516319}}}, // n-Butane,
            {"i-C4H10", {{20.413726078, -94.467620036, 3.06714, 8.97575, 5.25156, 25.14230, 16.13880}, {1.074673199, 0.485556021, 4.671261865, 2.191583480}}}, // isoButane
            {"n-C5H12", {{14.536611217, -89.919548319, 3.00000, 8.95043, 21.83600, 33.40320, 0}, {0.380391739, 1.789520971, 3.777411113, 0}}}, // n-Pentane
            {"i-C5H12", {{15.449907693, -101.298172792, 3.00000, 11.76180, 20.11010, 33.16880, 0}, {0.635392636, 1.977271641, 4.169371131, 0}}}, // isoPentane
            {"n-C6H14", {{14.345969349, -96.165722367, 3.00000, 11.69770, 26.81420, 38.61640, 0}, {0.359036667, 1.691951873, 3.596924107, 0}}}, // n-Hexane
            {"n-C7H16", {{15.063786601, -97.345252349, 3.00000, 13.72660, 30.47070, 43.55610, 0}, {0.314348398, 1.548136560, 3.259326458, 0}}}, // n-Heptane
            {"n-C8H18", {{15.864687161, -97.370667555, 3.00000, 15.68650, 33.80290, 48.17310, 0}, {0.279143540, 1.431644769, 2.973845992, 0}}}, // n-Octane
            {"n-C9H20", {{16.313913248, -102.160247463, 3.00000, 18.02410, 38.12350, 53.34150, 0}, {0.263819696, 1.370586158, 2.848860483, 0}}}, //n-Nonane
            {"n-C10H22", {{15.870791919, -108.858547525, 3.00000, 21.00690, 43.49310, 58.36570, 0}, {0.267034159, 1.353835195, 2.833479035, 0}}}, // n-Decane
            {"H2", {{13.796443393, -175.864487294, 1.47906, 0.95806, 0.45444, 1.56039, -1.37560}, {6.891654113, 9.847634830, 49.765290750, 50.367279301}}}, // Hydrogen
            {"O2", {{10.001843586, -14.996095135, 2.50146, 1.07558, 1.01334, 0, 0}, {14.461722565, 7.223325463, 0, 0}}}, // Oxygen
            {"CO", {{10.813340744, -19.834733959, 2.50055, 1.02865, 0.00493, 0, 0}, {11.669802800, 5.302762306, 0, 0}}}, // Carbon Monoxide
            {"H2O", {{8.203520690, -11.996306443, 3.00392, 0.01059, 0.98763, 3.06904, 0}, {0.415386589, 1.763895929, 3.874803739, 0}}}, // Water
            {"H2S", {{9.336197742, -16.266508995, 3.00000, 3.11942, 1.00243, 0, 0}, {4.914580541, 2.270653980, 0, 0}}}, // Hydrogen Sulfide
            {"He", {{13.628409737, -143.470759602, 1.50000, 0, 0, 0, 0}, {0 ,0, 0, 0}}}, // Helium
            {"Ar", {{8.316631500, -4.946502600, 1.50000, 0, 0, 0, 0}, {0 ,0, 0, 0}}} // Argon
        };
    }

    namespace Residual {

        struct PureFluid {
            std::vector<NP_t> const residualCoefficient;
            std::vector<std::pair<int, NP_t>> const polynomialExponent;
            std::vector<std::tuple<int, NP_t, int>> const exponentialExponent;
        };

         const static std::map<std::string_view, PureFluid const> pureFluid = {
            {
                "CH4", { // Methane
                    {0.57335704239162, -0.16760687523730e1, 0.23405291834916, -0.21947376343441, 0.16369201404128e-1, 0.15004406389280e-1, 0.98990489492918e-1, 0.58382770929055, -0.74786867560390, 0.30033302857974, 0.20985543806568, -0.18590151133061e-1, -0.15782558339049, 0.12716735220791, -0.32019743894346e-1, -0.68049729364536e-1, 0.24291412853736e-1, 0.51440451639444e-2, -0.19084949733532e-1, 0.55229677241291e-2, -0.44197392976085e-2, 0.40061416708429e-1, -0.33752085907575e-1, -0.25127658213357e-2},
                    {{1, 0.125}, {1, 1.125}, {2, 0.375}, {2, 1.125}, {4, 0.625}, {4, 1.500}},
                    {{1, 0.625, 1}, {1, 2.625, 1}, {1, 2.750, 1}, {2, 2.125, 1}, {3, 2.000, 1}, {6, 1.750, 1}, {2, 4.500, 2}, {3, 4.750, 2}, {3, 5.000, 2}, {4, 4.000, 2}, {4, 4.500, 2}, {2, 7.500, 3}, {3, 14.000, 3}, {4, 11.500, 3}, {5, 26.000, 6}, {6, 28.000, 6}, {6, 30.000, 6}, {7, 16.000, 6}}
                }
            }, 
            {
                "N2", { // Nitrogen
                    {0.59889711801201, -0.16941557480731e1, 0.24579736191718, -0.23722456755175, 0.17954918715141e-1, 0.14592875720215e-1, 0.10008065936206, 0.73157115385532, -0.88372272336366, 0.31887660246708, 0.20766491728799, -0.19379315454158e-1, -0.16936641554983, 0.13546846041701, -0.33066712095307e-1, -0.60690817018557e-1, 0.12797548292871e-1, 0.58743664107299e-2, -0.18451951971969e-1, 0.47226622042472e-2, -0.52024079680599e-2, 0.43563505956635e-1, -0.36251690750939e-1, -0.28974026866543e-2},
                    {{1, 0.125}, {1, 1.125}, {2, 0.375}, {2, 1.125}, {4, 0.625}, {4, 1.500}},
                    {{1, 0.625, 1}, {1, 2.625, 1}, {1, 2.750, 1}, {2, 2.125, 1}, {3, 2.000, 1}, {6, 1.750, 1}, {2, 4.500, 2}, {3, 4.750, 2}, {3, 5.000, 2}, {4, 4.000, 2}, {4, 4.500, 2}, {2, 7.500, 3}, {3, 14.000, 3}, {4, 11.500, 3}, {5, 26.000, 6}, {6, 28.000, 6}, {6, 30.000, 6}, {7, 16.000, 6}}
                }
            }, 
            {
                "CO2", { // Carbon Dioxide
                    {0.52646564804653, -0.14995725042592e1, 0.27329786733782, 0.12949500022786, 0.15404088341841, -0.58186950946814, -0.18022494838296, -0.95389904072812e-1, -0.80486819317679e-2, -0.35547751273090e-1, -0.28079014882405, -0.82435890081677e-1, 0.10832427979006e-1, -0.67073993161097e-2, -0.46827907600524e-2, -0.28359911832177e-1, 0.19500174744098e-1, -0.21609137507166, 0.43772794926972, -0.22130790113593, 0.15190189957331e-1, -0.15380948953300e-1},
                    {{1, 0.000}, {1, 1.250}, {2, 1.625}, {3, 0.375}},
                    {{3, 0.375, 1}, {3, 1.375, 1}, {4, 1.125, 1}, {5, 1.375, 1}, {6, 0.125, 1}, {6, 1.625, 1}, {1, 3.750, 2}, {4, 3.500, 2}, {1, 7.500, 3}, {1, 8.000, 3}, {3, 6.000 , 3}, {3, 16.000 , 6}, {4, 11.000, 3}, {5, 24.000, 5}, {5, 26.000, 5}, {5, 28.000, 5}, {5, 24.000, 6}, {5, 26.000, 6}}
                }
            },
            {
                "C2H6", { // Ethane
                    {0.63596780450714, -0.17377981785459e1, 0.28914060926272, -0.33714276845694, 0.22405964699561e-1, 0.15715424886913e-1, 0.11450634253745, 0.10612049379745e1, -0.12855224439423e1, 0.39414630777652, 0.31390924682041, -0.21592277117247e-1, -0.21723666564905, -0.28999574439489, 0.42321173025732, 0.46434100259260e-1, -0.13138398329741, 0.11492850364368e-1, -0.33387688429909e-1, 0.15183171583644e-1, -0.47610805647657e-2, 0.46917166277885e-1, -0.39401755804649e-1, -0.32569956247611e-2},
                    {{1, 0.125}, {1, 1.125}, {2, 0.375}, {2, 1.125}, {4, 0.625}, {4, 1.500}},
                    {{1, 0.625, 1}, {1, 2.625, 1}, {1, 2.750, 1}, {2, 2.125, 1}, {3, 2.000, 1}, {6, 1.750, 1}, {2, 4.500, 2}, {3, 4.750, 2}, {3, 5.000, 2}, {4, 4.000, 2}, {4, 4.500, 2}, {2, 7.500, 3}, {3, 14.000, 3}, {4, 11.500, 3}, {5, 26.000, 6}, {6, 28.000, 6}, {6, 30.000, 6}, {7, 16.000, 6}}
                }
            }, 
            {
                "C3H8", {  // Propane
                    {0.10403973107358e1, -0.28318404081403e1, 0.84393809606294, -0.76559591850023e1, 0.94697373057280e1, 0.24796475497006e-3, 0.27743760422870, -0.43846000648377e-1, -0.26991064784350, -0.69313413089860e-1, -0.29632145981653e-1, 0.14040126751380e-1},
                    {{1, 0.250}, {1, 1.125}, {1, 1.500}, {2, 1.375}, {3, 0.250}, {7, 0.875}},
                    {{2, 0.625, 1}, {5, 1.750, 1}, {1, 3.625, 2}, {4, 3.625, 2}, {3, 14.500, 3}, {4, 12.000, 3}}
                }
            },
            {
                "n-C4H10", { // n-Butane
                    {0.10626277411455e-1, -0.28620951828350e1, 0.88738233403777, -0.12570581155345, 0.10286308708106, 0.25358040602654e-3, 0.32325200233982, -0.37950761057432e-1, -0.32534802014452, -0.79050969051011e-1, -0.20636720547775e-1, 0.57053809334750e-2},
                    {{1, 0.250}, {1, 1.125}, {1, 1.500}, {2, 1.375}, {3, 0.250}, {7, 0.875}},
                    {{2, 0.625, 1}, {5, 1.750, 1}, {1, 3.625, 2}, {4, 3.625, 2}, {3, 14.500, 3}, {4, 12.000, 3}}
                }
            },
            {
                "i-C4H10", { // isoButane
                    {0.10429331589100e1, -0.28184272548892e1, 0.86176232397850, -0.10613619452487, 0.98615749302134e-1, 0.23948208682322e-3, 0.30330004856950, -0.41598156135099e-1, -0.29991937470058, -0.80369342764109e-1, -0.29761373251151e-1, 0.13059630303140e-1},
                    {{1, 0.250}, {1, 1.125}, {1, 1.500}, {2, 1.375}, {3, 0.250}, {7, 0.875}},
                    {{2, 0.625, 1}, {5, 1.750, 1}, {1, 3.625, 2}, {4, 3.625, 2}, {3, 14.500, 3}, {4, 12.000, 3}}
                }
            }, 
            {
                "n-C5H12", { // n-Pentane
                    {0.10968643098001e1, -0.29988888298061e1, 0.99516886799212, -0.16170708558539, 0.11334460072775, 0.26760595150748e-3, 0.40979881986931, -0.40876423083075e-1, -0.38169482469447, -0.10931956843993, -0.32073223327990e-1, 0.16877016216975e-1},
                    {{1, 0.250}, {1, 1.125}, {1, 1.500}, {2, 1.375}, {3, 0.250}, {7, 0.875}},
                    {{2, 0.625, 1}, {5, 1.750, 1}, {1, 3.625, 2}, {4, 3.625, 2}, {3, 14.500, 3}, {4, 12.000, 3}}
                }
            }, 
            {
                "i-C5H12", {  // isoPentane
                    {0.10963e1, -0.30402e1, 0.10317e1, -0.15410, 0.11535, 0.29809e-3, 0.39571, -0.45881e-1, -0.35804, -0.10107, -0.35484e-1, 0.18156e-1},
                    {{1, 0.250}, {1, 1.125}, {1, 1.500}, {2, 1.375}, {3, 0.250}, {7, 0.875}},
                    {{2, 0.625, 1}, {5, 1.750, 1}, {1, 3.625, 2}, {4, 3.625, 2}, {3, 14.500, 3}, {4, 12.000, 3}}
                }
            },
            {
                "n-C6H14", { // n-Hexane
                    {0.10553238013661e1, -0.26120615890629e1, 0.76613882967260, -0.29770320622459, 0.11879907733358, 0.27922861062617e-3, 0.46347589844105, 0.11433196980297e-1, -0.48256968738131, -0.93750558924659e-1, -0.67273247155994e-2, -0.51141583585428e-2},
                    {{1, 0.250}, {1, 1.125}, {1, 1.500}, {2, 1.375}, {3, 0.250}, {7, 0.875}},
                    {{2, 0.625, 1}, {5, 1.750, 1}, {1, 3.625, 2}, {4, 3.625, 2}, {3, 14.500, 3}, {4, 12.000, 3}}
                }
            }, 
            {
                "n-C7H16", { // n-Heptane
                    {0.10543747645262e1, -0.26500681506144, 0.81730047827543, -0.30451391253428, 0.12253868710800, 0.27266472743928e-3, 0.49865825681670, -0.71432815084176e-3, -0.54236895525450, -0.13801821610756, -0.61595287380011e-2, 0.48602510393022e-3},
                    {{1, 0.250}, {1, 1.125}, {1, 1.500}, {2, 1.375}, {3, 0.250}, {7, 0.875}},
                    {{2, 0.625, 1}, {5, 1.750, 1}, {1, 3.625, 2}, {4, 3.625, 2}, {3, 14.500, 3}, {4, 12.000, 3}}
                }
            }, 
            {
                "n-C8H18", { // n-Octane
                    {0.10722544875633e1, -0.24632951172003e1, 0.65386674054928, -0.36324974085628, 0.12713269626764, 0.30713572777930e-3, 0.52656856987540, 0.19362862857653e-1, -0.58939426849155, -0.14069963991934, -0.78966330500036e-2, 0.33036597968109e-2},
                    {{1, 0.250}, {1, 1.125}, {1, 1.500}, {2, 1.375}, {3, 0.250}, {7, 0.875}},
                    {{2, 0.625, 1}, {5, 1.750, 1}, {1, 3.625, 2}, {4, 3.625, 2}, {3, 14.500, 3}, {4, 12.000, 3}}
                }
            }, 
            {
                "n-C9H20", { //n-Nonane
                    {0.11151e1, -0.27020e1, 0.83416, -0.38828, 0.13760, 0.28185e-3, 0.62037, 0.15847e-1, -0.61726, -0.15043, -0.12982e-1, 0.44325e-2},
                    {{1, 0.250}, {1, 1.125}, {1, 1.500}, {2, 1.375}, {3, 0.250}, {7, 0.875}},
                    {{2, 0.625, 1}, {5, 1.750, 1}, {1, 3.625, 2}, {4, 3.625, 2}, {3, 14.500, 3}, {4, 12.000, 3}}
                }
            }, 
            {
                "n-C10H22", { // n-Decane
                    {0.10461e1, -0.24807e1, 0.74372, -0.52579, 0.15315, 0.32865e-3, 0.84178, 0.55424e-1, -0.73555, -0.18507, -0.20775e-1, 0.12335e-1},
                    {{1, 0.250}, {1, 1.125}, {1, 1.500}, {2, 1.375}, {3, 0.250}, {7, 0.875}},
                    {{2, 0.625, 1}, {5, 1.750, 1}, {1, 3.625, 2}, {4, 3.625, 2}, {3, 14.500, 3}, {4, 12.000, 3}}
                }
            }, 
            {
                "H2", { // Hydrogen
                    {0.53579928451252e1, -0.62050252530595e1, 0.13830241327086, -0.71397954896129e-1, 0.15474053959733e-1, -0.14976806405771, -0.26368723988451e-1, 0.56681303156066e-1, -0.60063958030436e-1, -0.45043942027132, 0.42478840244500, -0.21997640827139e-1, -0.10499521374530e-1, -0.28955902866816e-2},
                    {{1, 0.500}, {1, 0.625}, {2, 0.375}, {2, 0.625}, {4, 1.125}},
                    {{1, 2.625, 1}, {5, 0.000, 1}, {5, 0.250, 1}, {5, 1.375, 1}, {1, 4.000, 2}, {1, 4.250, 2}, {2, 5.000, 3}, {5, 8.000, 3}, {1, 8.000, 5}}
                }
            }, 
            {
                "O2", { // Oxygen
                    {0.88878286369701, -0.24879433312148e1, 0.59750190775886, 0.96501817061881e-2, 0.71970428712770e-1, 0.22337443000195e-3, 0.18558686391474, -0.38129368035760e-1, -0.15352245383006, -0.26726814910919e-1, -0.25675298677127e-1, 0.95714302123668e-2},
                    {{1, 0.250}, {1, 1.125}, {1, 1.500}, {2, 1.375}, {3, 0.250}, {7, 0.875}},
                    {{2, 0.625, 1}, {5, 1.750, 1}, {1, 3.625, 2}, {4, 3.625, 2}, {3, 14.500, 3}, {4, 12.000, 3}}
                }
            }, 
            {
                "CO", { // Carbon Monoxide
                    {0.90554, -0.24515e1, 0.53149, 0.24173e-1, 0.72156e-1, 0.18818e-3, 0.19405, -0.43268e-1, -0.12778, -0.27896e-1, -0.34154e-1, 0.16329e-1},
                    {{1, 0.250}, {1, 1.125}, {1, 1.500}, {2, 1.375}, {3, 0.250}, {7, 0.875}},
                    {{2, 0.625, 1}, {5, 1.750, 1}, {1, 3.625, 2}, {4, 3.625, 2}, {3, 14.500, 3}, {4, 12.000, 3}}
                }
            }, 
            {
                "H2O", { // Water
                    {0.82728408749586, -0.18602220416584e1, -0.11199009613744e1, 0.15635753976056, 0.87375844859025, -0.36674403715731, 0.53987893432436e-1, 0.10957690214499e1, 0.53213037828563e-1, 0.13050533930825e-1, -0.41079520434476, 0.14637443344120, -0.55726838623719e-1, -0.11201774143800e-1, -0.66062758068099e-2, 0.46918522004538e-2},
                    {{1, 0.500}, {1, 1.250}, {1, 1.875}, {2, 0.125}, {2, 1.500}, {3, 1.000}, {4, 0.750}},
                    {{1, 1.500, 1}, {5, 0.625, 1}, {5, 2.625, 1}, {1, 5.000, 2}, {2, 4.000, 2}, {4, 4.500, 2}, {4, 3.000, 3}, {1, 4.000, 5}, {1, 6.000, 5}}
                }
            }, 
            {
                "H2S", { // Hydrogen Sulfide
                    {0.87641, -0.20367e1, 0.21634, -0.50199e-1, 0.66994e-1, 0.19076e-3, 0.20227, -0.45348e-2, -0.22230, -0.34714e-1, -0.14885e-1, 0.74154e-2},
                    {{1, 0.250}, {1, 1.125}, {1, 1.500}, {2, 1.375}, {3, 0.250}, {7, 0.875}},
                    {{2, 0.625, 1}, {5, 1.750, 1}, {1, 3.625, 2}, {4, 3.625, 2}, {3, 14.500, 3}, {4, 12.000, 3}}
                }
            }, 
            {
                "He", { // Helium
                    {-0.45579024006737, 0.12516390754925e1, -0.15438231650621e1, 0.20467489707221e-1, -0.34476212380781, -0.20858459512787e-1, 0.16227414711778e-1, -0.57471818200892e-1, 0.19462416430715e-1, -0.33295680123020e-1, -0.10863577372367e-1, -0.22173365245954e-1},
                    {{1, 0.000}, {1, 0.125}, {1, 0.750}, {4, 1.000}},
                    {{1, 0.750, 1}, {3, 2.625, 1}, {5, 0.125, 1}, {5, 1.250, 1}, {5, 2.000, 1}, {2, 1.000, 2}, {1, 4.500, 3}, {2, 5.000, 2}}
                }
            }, 
            {
                "Ar", {// Argon
                    {0.85095714803969, -0.24003222943480e1, 0.54127841476466, 0.16919770692538e-1, 0.68825965019035e-1, 0.21428032815338e-3, 0.17429895321992, -0.33654495604194e-1, -0.13526799857691, -0.16387350791552e-1, -0.24987666851475e-1, 0.88769204815709e-2},
                    {{1, 0.250}, {1, 1.125}, {1, 1.500}, {2, 1.375}, {3, 0.250}, {7, 0.875}},
                    {{2, 0.625, 1}, {5, 1.750, 1}, {1, 3.625, 2}, {4, 3.625, 2}, {3, 14.500, 3}, {4, 12.000, 3}}
                } 
            }
         };
    }
}

#endif /* PUREFLUIDHELMHOLTZCOEFFICIENTS_HPP */