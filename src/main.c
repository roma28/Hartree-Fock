// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.
//
// Copyright (c) Roman Ishchenko 2022.

#include <stdio.h>
#include <string.h>
#include <gsl/gsl_vector.h>

#include "../include/atom.h"
#include "../include/integrals.h"


int main() {


    atom *a = atom_alloc();
    atom *b = atom_alloc();
    a->Z = 1;
    b->Z = 1;

    // STO-3G for hydrogen
    double e[3] = {0.3425250914E+01, 0.6239137298E+00, 0.1688554040E+00};
    double c[3] = {0.1543289673E+00, 0.5353281423E+00, 0.4446345422E+00};


    a->basis = atom_basis_alloc(1);
    b->basis = atom_basis_alloc(1);
    for (size_t i = 0; i < a->basis->n_contracted; ++i) {
        a->basis->basis_functions[i] = basis_function_alloc(3);
        b->basis->basis_functions[i] = basis_function_alloc(3);
    }

    a->basis->basis_functions[0]->L = L_S;
    b->basis->basis_functions[0]->L = L_S;


    memcpy(a->basis->basis_functions[0]->exponents->data, e, 3 * sizeof(double));
    memcpy(a->basis->basis_functions[0]->contractions->data, c, 3 * sizeof(double));
    memcpy(b->basis->basis_functions[0]->exponents->data, e, 3 * sizeof(double));
    memcpy(b->basis->basis_functions[0]->contractions->data, c, 3 * sizeof(double));


    gsl_vector *ca = gsl_vector_alloc(3);

    for (size_t i = 0; i < 3; ++i) {
        ca->data[i] = 0.1 * i;
    }

    gsl_vector_memcpy(a->basis->basis_functions[0]->origin, ca);

    double s;
    double t;
    double v;
    double g;

    for (size_t i = 0; i < 1000; ++i) {
        s = overlap_integral(a->basis->basis_functions[0], a->basis->basis_functions[0]);
        t = kinetic_energy_integral(a->basis->basis_functions[0], a->basis->basis_functions[0]);
        v = nuclear_repulsion_integral(a->basis->basis_functions[0], a->basis->basis_functions[0], &a->Z, &ca, 1);
        g = electron_repulsion_integral(a->basis->basis_functions[0], a->basis->basis_functions[0],
                                        a->basis->basis_functions[0],
                                        a->basis->basis_functions[0]);

    }

    printf("S = %e\nT = %e\nV = %e\nG = %e", s, t, v, g);
//    char str[255];

    return 0;
}
