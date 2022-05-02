// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.
//
// Copyright (c) Roman Ishchenko 2022.

#include "../include/primitive_integrals.h"
#include <string.h>
#include "../include/integral_tools.h"
#include <gsl/gsl_math.h>


double s00(double exp_a, double exp_b, double deltaR2) {
    double Na = N(0, 0, 0, exp_a);
    double Nb = N(0, 0, 0, exp_b);

    double q = exp_a * exp_b / (exp_a + exp_b);

    double t1 = sqrt(gsl_pow_3(M_PI / (exp_a + exp_b)));

    double e = exp(-q * deltaR2);

    return Na * Nb * t1 * e;
}

//double E(uint8_t l1, uint8_t l2, int8_t t, double exp_a, double exp_b, double deltaR2) {
//    if ((t < 0) || (t > l1 + l2)) return 0.0;
//    else if ((l1 == 0) && (l2 == 0) && (t == 0)) {
//        return s00(exp_a, exp_b, deltaR2);
//    }
//    return -1;
//}



double k00(double exp_a, double exp_b, double deltaR2) {
    double Na = N(0, 0, 0, exp_a);
    double Nb = N(0, 0, 0, exp_b);

    double q = exp_a * exp_b / (exp_a + exp_b);
    return Na * Nb * (3 * q - 2 * gsl_pow_2(q) * deltaR2);
}


