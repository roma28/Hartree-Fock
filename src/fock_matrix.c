// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.
//
// Copyright (c) Roman Ishchenko 2022.

#include "../include/fock_matrix.h"
#include "../include/integrals.h"
//#include "../include/atom.h"

gsl_matrix *fock_matrix(const atom **atoms, const uint32_t natoms) {
    int nfunc = 0;
    gsl_vector **origins = malloc(natoms * sizeof(gsl_vector *));
    int *Z = malloc(natoms * sizeof(int));

    for (size_t i = 0; i < natoms; ++i) {
        nfunc += atoms[i]->basis->n_contracted;
        origins[i] = atoms[i]->coords;
        Z[i] = atoms[i]->Z;
    }


    free(origins);
    free(Z);
}

gsl_matrix *overlap_matrix(const atom **atoms, const uint32_t nfunc) {
    gsl_matrix *s = gsl_matrix_alloc(nfunc, nfunc);
    for (size_t i = 0; i < nfunc; ++i) {
        for (size_t j = 0; j < nfunc; ++j) {
            for (size_t a = 0; a < atoms[i]->basis->n_contracted; ++a) {
                for (size_t b = 0; b < atoms[j]->basis->n_contracted; ++b) {
                    gsl_matrix_set(s, i * atoms[i]->basis->n_contracted + a, j * atoms[j]->basis->n_contracted + b,
                                   overlap_integral(atoms[i]->basis->basis_functions[a],
                                                    atoms[j]->basis->basis_functions[b]));
                }
            }
        }
    }

    return s;
}


void nuclear_repulsion_matrix(gsl_matrix *t, const atom **atoms, const uint32_t natoms) {
    for (size_t i = 0; i < natoms; ++i) {
        for (size_t j = 0; j < natoms; ++j) {
            for (size_t a = 0; a < atoms[i]->basis->n_contracted; ++a) {
                for (size_t b = 0; b < atoms[j]->basis->n_contracted; ++b) {
                    gsl_matrix_set(t, i * atoms[i - 1]->basis->n_contracted + a,
                                   j * atoms[j - 1]->basis->n_contracted + b,
                                   nuclear_repulsion_integral(atoms[i]->basis->basis_functions[a],
                                                              atoms[j]->basis->basis_functions[b]));
                }
            }
        }
    }
}


