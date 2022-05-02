// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.
//
// Copyright (c) Roman Ishchenko 2022.

#ifndef SCF_PRIMITIVE_INTEGRALS_H
#define SCF_PRIMITIVE_INTEGRALS_H
#include "atom_basis.h"
#include <string.h>
#include "integral_tools.h"

/**
 * Computes overlap integral of two 3D gaussian primitives of L=0
 * @param exp_a exponent of first gaussian primitive
 * @param exp_b exponent of second gaussian primitive
 * @param deltaR2 squared distance between two gaussians' origins
 * @return S00 integral value
 */
double s00(double exp_a, double exp_b, double deltaR2);

double k00(double exp_a, double exp_b, double deltaR2);

///**
// * Computes S00 type integral of two 3D gaussian primitives of L=0
// * @param l1 exponent of first gaussian primitive
// * @param l2 exponent of second gaussian primitive
// * @param t
// * @param deltaR2 squared distance between two gaussians' origins
// * @return S00 integral value
// *
// * @todo recurrent expansion
// */
//double E(uint8_t l1, uint8_t l2, int8_t t, double exp_a, double exp_b, double deltaR2);

#endif //SCF_PRIMITIVE_INTEGRALS_H
