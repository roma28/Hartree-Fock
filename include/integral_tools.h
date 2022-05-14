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
double N(uint8_t i, uint8_t j, uint8_t k, double alpha);

double deltaR2(const gsl_vector *Ra, const gsl_vector *Rb);


#endif //SCF_INTEGRAL_TOOLS_H
