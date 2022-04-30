// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.
//
// Copyright (c) Roman Ishchenko 2022.

#include "../include/basis_reader.h"

int a[255];

void parse_block(FILE *f) {
    char s[255];
    fgets(s, 255, f);
    char moment;
    int nf;
    float norm;
    sscanf(s, "%c\t%d\t%f", &moment, &nf, &norm);

    basis_function *bf = basis_function_alloc(nf);

    int n;

    for (size_t i = 0; i < nf; ++i) {
        fgets(s, 255, f);
        sscanf(s, "\t%le\t%le\n", &(bf->exponents->data[i]), &(bf->contractions->data[i]));
    }

    return;

}

void parse_basis_file(FILE *f) {
    char s[255];
    do {
        fgets(s, 255, f);
        if (s[0] == '!') break;
        if (s == "****") break;
        char atom_symbol;
        sscanf(s, "%c\t", &atom_symbol);
        parse_block(f);
    } while (strnlen(s, 255) > 0);

}


