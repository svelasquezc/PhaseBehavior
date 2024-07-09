#ifndef GERG_COMPONENTS_HPP
#define GERG_COMPONENTS_HPP

#include "../Component.hpp"

namespace PhaseBehavior::EoS::GERG::Components {
                                                 /*[kPa]*/  /*[K]*/    /*  [m3/kg]  */ /*[g/mol]*/ /*[-]*/
    static Component const CH4      = {"CH4",      4599,    190.564,   6.1477929425e-3, 16.04246,  0.011};
    static Component const N2       = { "N2",      3390,    126.192,   3.1918384511e-3, 28.0134,   0.039};
    static Component const CO2      = {"CO2",      7380,    304.1282,  2.1385799828e-3, 44.0095,   0.239};
    static Component const C2H6     = {"C2H6",     4900,    305.322,   4.8402710552e-3, 30.06904,  0.099};
    static Component const C3H8     = {"C3H8",     4250,    369.825,   4.5355587810e-3, 44.09562,  0.153};
    static Component const n_C4H10  = {"n-C4H10",  3796,    425.125,   4.3890449440e-3, 58.1222,   0.199};
    static Component const i_C4H10  = {"i-C4H10",  3647,    407.817,   4.4571224820e-3, 58.1222,   0.183};
    static Component const n_C5H12  = {"n-C5H12",  3369,    469.7,     4.3103448281e-3, 72.14878,  0.251};
    static Component const i_C5H12  = {"i-C5H12",  3380,    460.35,    4.2373122060e-3, 72.14878,  0.227};
    static Component const n_C6H14  = {"n-C6H14",  3020,    507.82,    4.2885324643e-3, 86.17536,  0.299};
    static Component const n_C7H16  = {"n-C7H16",  2736,    540.13,    4.3103448273e-3, 100.20194, 0.349};
    static Component const n_C8H18  = {"n-C8H18",  2486,    569.32,    4.2571306948e-3, 114.22852, 0.398};
    static Component const n_C9H20  = {"n-C9H20",  2300,    594.55,    4.3077132047e-3, 128.2551,  0.4433};
    static Component const n_C10H22 = {"n-C10H22", 2110,    617.7,     4.2855559307e-3, 142.28168, 0.4884};
    static Component const H2       = {"H2",       1300,    33.19,     3.3203565829e-2, 2.01588,  -0.216};
    static Component const O2       = {"O2",       5043,    154.595,   2.2928225913e-3, 31.9988,   0.025};
    static Component const CO       = {"CO",       3498.75, 132.86,    3.2904523232e-3, 28.0101,   0.066};
    static Component const H2O      = {"H2O",      22064,   647.096,   3.1055900621e-3, 18.01528,  0.344};
    static Component const H2S      = {"H2S",      8970,    373.1,     2.8794862952e-3, 34.08088,  0.081};
    static Component const He       = {"He",       227.4,   5.1953,    1.4359301151e-2, 4.002602, -0.365};
    static Component const Ar       = {"Ar",       4897.9,  150.687,   1.8670649738e-3, 39.948,    0.001};
}

#endif /* GERG_COMPONENTS_HPP */