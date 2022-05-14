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
    a->Z = 1;

    // STO-3G for hydrogen
    double e[3] = {0.3425250914E+01, 0.6239137298E+00, 0.1688554040E+00};
    double c[3] = {0.1543289673E+00, 0.5353281423E+00, 0.4446345422E+00};


    a->basis = atom_basis_alloc(1);
    for (size_t i = 0; i < a->basis->n_contracted; ++i) {
        a->basis->basis_functions[i] = basis_function_alloc(3);
    }

    a->basis->basis_functions[0]->L = L_S;


    memcpy(a->basis->basis_functions[0]->exponents->data, e, 3 * sizeof(double));
    memcpy(a->basis->basis_functions[0]->contractions->data, c, 3 * sizeof(double));

    gsl_vector *ca = gsl_vector_alloc(3);

    for (size_t i = 0; i < 3; ++i) {
        ca->data[i] = 0.1 * i;
    }

    gsl_vector_memcpy(a->basis->basis_functions[0]->origin, ca);

    for (size_t i = 0; i < 100000; ++i) {
        double s = S(a->basis->basis_functions[0], a->basis->basis_functions[0]);
        double t = T(a->basis->basis_functions[0], a->basis->basis_functions[0]);
        double v = V(a->basis->basis_functions[0], a->basis->basis_functions[0], &a->Z, &ca, 1);
        double g = G(a->basis->basis_functions[0], a->basis->basis_functions[0], a->basis->basis_functions[0],
                     a->basis->basis_functions[0]);

    }

//    printf("S = %e\nT = %e\nV = %e\nG = %e", s, t, v, g);
    char str[255];

    return 0;
}
