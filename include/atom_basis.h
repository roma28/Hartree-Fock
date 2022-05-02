// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.
//
// Copyright (c) Roman Ishchenko 2022.


#ifndef SCF_ATOM_BASIS_H
#define SCF_ATOM_BASIS_H

#include <gsl/gsl_vector.h>

#define L_S 0u;
#define L_P 1u;

typedef struct {
    size_t L;
    double exponent;
    gsl_vector *origin;
    double N;
} primitive;

/**
 * \struct Stores a contracted basis_functions function consisting of gaussian primitives
 */
typedef struct {
    size_t L; ///< angular momentum
    size_t n_primitives; ///< number of gaussian primitives
    gsl_vector *exponents; ///< exponents of gaussian primitives
    gsl_vector *contractions; ///< contraction coefficients of gaussian primitives
    gsl_vector *origin;
} basis_function;

/**
 * \struct Stores a basis_functions set for an atom
 */
typedef struct {
    size_t n_contracted; ///< number of contracted functions
    basis_function **basis_functions; ///< contracted basis functions of the atom
} atom_basis;

/**
 * Allocates basis_functions structure consisting of n gaussian primitives
 * @param n number of gaussian primitives
 * @return pointer to basis function
 */
basis_function *basis_function_alloc(size_t n);

void basis_function_free(basis_function *p);

/**
 * Allocates atom_basis structure consisting of n contracted basis functions
 * @param n number of gaussian primitives
 * @return pointer to atom_basis
 */
atom_basis *atom_basis_alloc(size_t n);

void atom_basis_free(atom_basis *p);

#endif //SCF_ATOM_BASIS_H
