// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.
//
// Copyright (c) Roman Ishchenko 2022.

//
// Created by Roman Ishchenko on 30.04.2022.
//

#ifndef SCF_PRIMITIVE_INTEGRALS_H
#define SCF_PRIMITIVE_INTEGRALS_H
#include "atom_basis.h"
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>

/**
 * Computes normalization constant of 3D gaussian primitive
 * @param alpha exponent of gaussian primitive
 * @return normalization constant
 */
double N(double alpha);

double Q(const gsl_vector *Ra, const gsl_vector *Rb);

/**
 * Computes S00 type integral of two 3D gaussian primitives of L=0
 * @param exp_a exponent of first gaussian primitive
 * @param exp_b exponent of second gaussian primitive
 * @param Q squared distance between two gaussians' origins
 * @return S00 integral value
 */
double s00( double exp_a, double exp_b, double Q);

double E(uint8_t i, uint8_t j, uint8_t t, double exp_a, double exp_b, double Q);

/**
 * Computes overlap integral of two contracted basis functions
 * @param a pointer to first function
 * @param b pointer to second function
 * @return overlap integral values
 */
double S(const basis_function *a, const basis_function *b);

#endif //SCF_PRIMITIVE_INTEGRALS_H
