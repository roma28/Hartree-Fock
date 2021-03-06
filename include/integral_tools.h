// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.
//
// Copyright (c) Roman Ishchenko 2022.

#ifndef SCF_INTEGRAL_TOOLS_H
#define SCF_INTEGRAL_TOOLS_H

#include <stdint.h>
#include <gsl/gsl_vector.h>


/**
 * Computes normalization constant of 3D gaussian primitive @n
 * of form  y = x^i * y^j * z^k * exp(-alpha * R^2)
 * @param alpha exponent of gaussian primitive
 * @param i x cartesian quantum number
 * @param j y cartesian quantum number
 * @param k z cartesian quantum number
 * @return normalization constant
 */
double normalization_constant(uint8_t i, uint8_t j, uint8_t k, double alpha);

/**
 * Computes normalization constant of 3D gaussian primitive @n
 * of form  y = x^i * y^j * z^k * exp(-alpha * R^2)
 * @param Ra first gaussian origin (gsl_vector with 3 elements)
 * @param Rb second gaussian origin (gsl_vector with 3 elements)
 * @return squared distance between two origins
 */
double deltaR2(const gsl_vector *Ra, const gsl_vector *Rb);

double boys_function(double n, double T);

gsl_vector *gaussian_center(const gsl_vector *Ra, const gsl_vector *Rb, double exp_a, double exp_b);

#endif //SCF_INTEGRAL_TOOLS_H
