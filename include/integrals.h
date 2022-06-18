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
 * @return overlap integral value
 */
double overlap_integral(const basis_function *a, const basis_function *b);

/**
 * Computes kinetic energy integral of two contracted basis functions
 * @param a pointer to first function
 * @param b pointer to second function
 * @return kinetic integral value
 */
double kinetic_energy_integral(const basis_function *a, const basis_function *b);

/**
 * Computes nuclear repulsion integral of two contracted basis functions with all the nuclei
 * @param a pointer to first function
 * @param b pointer to second function
 * @param Z nuclei charges array
 * @param atom_origins nuclei coordinates array
 * @param natoms number of atoms
 * @return nuclear repulsion integral value
 */
double
nuclear_repulsion_integral(const basis_function *a, const basis_function *b, const uint8_t *Z,
                           gsl_vector **atom_origins, size_t natoms);

/**
 * Computes nuclear repulsion integral of two contracted basis functions with all the nuclei
 * @param a pointer to first function
 * @param b pointer to second function
 * @param x pointer to first function
 * @param y pointer to second function
 * @return kinetic integral value
 */
double electron_repulsion_integral(const basis_function *a, const basis_function *b, const basis_function *x,
                                   const basis_function *y);

#endif //SCF_INTEGRALS_H
