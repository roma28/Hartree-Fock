// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.
//
// Copyright (c) Roman Ishchenko 2022.

//
// Created by Roman Ishchenko on 02.05.2022.
//

#include "../include/integral_tools.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>


double N(uint8_t i, uint8_t j, uint8_t k, double alpha) {
    double eN = pow(gsl_pow_3(M_2_PI * alpha), 0.25);
    if ((i == 0) && (j == 0) && (k == 0)) { return eN; }
    else {
        double fac = sqrt(pow(4 * alpha, i + j + k) / gsl_sf_doublefact(2 * i - 1) /
                          gsl_sf_doublefact(2 * j - 1) / gsl_sf_doublefact(2 * k - 1));

        return eN * fac;
    }
}

double deltaR2(const gsl_vector *Ra, const gsl_vector *Rb) {
    gsl_vector *res = gsl_vector_alloc(3);
    gsl_vector_memcpy(res, Ra);
    gsl_vector_sub(res, Rb);

    double dR2;
    gsl_blas_ddot(res, res, &dR2);
    gsl_vector_free(res);

    return dR2;
}

double boys_function_hyperg(double n, double T) {
    return gsl_sf_hyperg_1F1(n + 0.5, n + 1.5, T);
}

double boys_function_incomplete_gamma(double n, double T) {
    return pow(2 * T, -(n + 0.5)) * gsl_sf_gamma_inc(n + 0.5, T);
}
