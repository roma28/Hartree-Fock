// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.
//
// Copyright (c) Roman Ishchenko 2022.

#include "../include/integral_tools.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>


double normalization_constant(uint8_t i, uint8_t j, uint8_t k, double alpha) {
    double eN = pow(gsl_pow_3(M_2_PI * alpha), 0.25);
    if ((i == 0) && (j == 0) && (k == 0)) { return eN; }
    else {
        double fac = sqrt(pow(4 * alpha, i + j + k) / gsl_sf_doublefact(2 * i - 1) /
                          gsl_sf_doublefact(2 * j - 1) / gsl_sf_doublefact(2 * k - 1));

        return eN * fac;
    }
}

double deltaR2(const gsl_vector *Ra, const gsl_vector *Rb) {
    double dR2 = 0;
    for (size_t i = 0; i < 3; ++i) {
        dR2 += gsl_pow_2(Ra->data[i] - Rb->data[i]);
    }
    return dR2;
}

double boys_function(double n, double T) {
    return gsl_sf_hyperg_1F1(n + 0.5, n + 1.5, T);
}


gsl_vector *gaussian_center(const gsl_vector *Ra, const gsl_vector *Rb, double exp_a, double exp_b) {
    gsl_vector *res = gsl_vector_alloc(3);

    for (size_t i = 0; i < 3; ++i) {
        res->data[i] = (exp_a * Ra->data[i] + exp_b * Rb->data[i]) / (exp_a + exp_b);
    }

    return res;
}
