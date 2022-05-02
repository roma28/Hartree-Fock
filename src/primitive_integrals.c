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
    return pow(gsl_pow_3(M_2_PI * alpha), 0.25);
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

double q(const double exp_a, const double exp_b) {
    return -exp_a * exp_b / (exp_a + exp_b);
}


double s00(double exp_a, double exp_b, double Q) {
    double Na = N(exp_a);
    double Nb = N(exp_b);

    double t1 = sqrt(gsl_pow_3(M_PI / (exp_a + exp_b)));

    double e = exp(q(exp_a, exp_b) * Q);

    return Na * Nb * t1 * e;

}

double E(uint8_t i, uint8_t j, uint8_t t, double exp_a, double exp_b, double Q) {
    if ((t < 0) || (t > i + j)) return 0.0;
    else if ((i == 0) && (j == 0) && (t == 0)) {
        return s00(exp_a, exp_b, Q);
    } else if (j == 0) {
    }
    return -1;
}


double S(const basis_function *a, const basis_function *b) {
    double S = 0.0;
    double q = Q(a->origin, b->origin);

    for (size_t i = 0; i < a->n_primitives; ++i) {
        for (size_t j = 0; j < b->n_primitives; ++j) {
            S += a->contractions->data[i] * b->contractions->data[j] *
                 s00(a->exponents->data[i], b->exponents->data[j], q);
        }
    }

    return S;
}
