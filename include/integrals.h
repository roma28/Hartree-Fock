// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.
//
// Copyright (c) Roman Ishchenko 2022.

#ifndef SCF_INTEGRALS_H
#define SCF_INTEGRALS_H

#include "primitive_integrals.h"
#include "integral_tools.h"

/**
 * Computes overlap integral of two contracted basis functions
 * @param a pointer to first function
 * @param b pointer to second function
 * @return overlap integral values
 */
double S(const basis_function *a, const basis_function *b);

/**
 * Computes kinetic integral of two contracted basis functions
 * @param a pointer to first function
 * @param b pointer to second function
 * @return kinetic integral values
 */
double T(const basis_function *a, const basis_function *b);


#endif //SCF_INTEGRALS_H
