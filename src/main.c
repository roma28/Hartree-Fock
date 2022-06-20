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
//#include <gsl/gsl_vector.h>
#include "../include/atom.h"
#include "../include/integrals.h"
#include "../include/xyz_reader.h"
#include "c-logger/logger.h"
#include "../include/system.h"

print_system(const System *s) {
    printf("%d atoms in total\n", s->natoms);
    for (size_t i = 0; i < s->natoms; ++i) {
        printf("%d\t%d\t%f\t%f\t%f\n", i, s->atoms[i]->Z, s->atoms[i]->coords->data[0], s->atoms[i]->coords->data[1],
               s->atoms[i]->coords->data[2]);
    }
}


int main() {
    logger_initConsoleLogger(stdout);
    logger_setLevel(LogLevel_DEBUG);

    LOG_INFO("Started");

    System *s = malloc(sizeof(System));
    if (!s) {
        LOG_FATAL("Allocation error");
        return -1;
    }

    // STO-3G for hydrogen
    double e[3] = {0.3425250914E+01, 0.6239137298E+00, 0.1688554040E+00};
    double c[3] = {0.1543289673E+00, 0.5353281423E+00, 0.4446345422E+00};

    atom_basis *h_bas = atom_basis_alloc(1);
    h_bas->basis_functions[0] = basis_function_alloc(3);
    memcpy(h_bas->basis_functions[0]->exponents, e, 3 * sizeof(double));
    memcpy(h_bas->basis_functions[0]->contractions, c, 3 * sizeof(double));


    FILE *f = fopen("../struct.xyz", "r");
    if (f) {
        read_xyz(f, s);
    } else {
        LOG_ERROR("File opening error");
    }

    print_system(s);

    return 0;
}
