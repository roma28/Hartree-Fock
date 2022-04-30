// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.
//
// Copyright (c) Roman Ishchenko 2022.

#include "../include/primitive_integrals.h"

double N(double alpha) {
    return sqrt(sqrt(gsl_pow_3(2 * alpha * M_1_PI)));
}

double Q(const gsl_vector *Ra, const gsl_vector *Rb) {
    gsl_vector *res = gsl_vector_alloc(3);
    gsl_vector_memcpy(res, Ra);
    gsl_vector_sub(res, Rb);
    double Q;
    gsl_blas_ddot(res, res, &Q);
    gsl_vector_free(res);

    return Q;
}

double s00(const double exp_a, const double exp_b, const gsl_vector *Ra, const gsl_vector *Rb) {
//    double Na = N(exp_a);
//    double Nb = N(exp_b);

    double Na = 1l;
    double Nb = 1l;


    double t1 = gsl_pow_3(M_SQRTPI / sqrt(exp_a + exp_b));

    double q = -exp_a * exp_b / (exp_a + exp_b);

    double e = exp(q * Q(Ra, Rb));

    return Na * Nb * t1 * e;

}
