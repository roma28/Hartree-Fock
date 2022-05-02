// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.
//
// Copyright (c) Roman Ishchenko 2022.

#include "../include/atom_basis.h"

basis_function *basis_function_alloc(size_t n) {
    basis_function *f = malloc(sizeof(basis_function));
    if (f) {
        f->n_primitives = n;
        f->exponents = gsl_vector_alloc(n);
        f->contractions = gsl_vector_alloc(n);
        f->origin = gsl_vector_alloc(3);
    }
    return f;
}

void basis_function_free(basis_function *p) {
    gsl_vector_free(p->exponents);
    gsl_vector_free(p->contractions);
    gsl_vector_free(p->origin);
    free(p);
    p = NULL;
}

atom_basis *atom_basis_alloc(size_t nfunc) {
    atom_basis *b = malloc(sizeof(atom_basis));
    if (b) {
        b->n_contracted = nfunc;
        b->basis_functions = (basis_function**) malloc(nfunc * sizeof(basis_function *));
    }
    return b;
}

void atom_basis_free(atom_basis *p) {
    for (size_t i = 0; i < p->n_contracted; ++i) {
        basis_function_free(p->basis_functions[i]);
    }
    free(p->basis_functions);
    free(p);
    p = NULL;
}