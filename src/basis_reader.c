// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.
//
// Copyright (c) Roman Ishchenko 2022.

#include "../include/basis_reader.h"


void parse_block(FILE *f, basis_function *bf) {
    char s[255];

    fgets(s, 255, f);
    if (s[0] == '\0') return;
    char atom_symbol;
    sscanf(s, "%c\t", &atom_symbol);

    fgets(s, 255, f);
    char moment;
    int n_p;

    sscanf(s, "%c\t%d\t", &moment, &n_p);
    bf = basis_function_alloc(n_p);

    for (size_t i = 0; i < n_p; ++i) {
        fgets(s, 255, f);
        sscanf(s, GAUSSIAN_PRIMITIVE_FORMAT_STRING, &(bf->exponents->data[i]), &(bf->contractions->data[i]));
    }

    return;

}


basis_function **parse_basis_file(FILE *f) {
    char s[255];
    do {
        fgets(s, 255, f);
        if (s[0] == COMMENT_START) break;
        if (strcmp(s, BLOCK_START) != 0) break;
        basis_function *bf;
        parse_block(f, bf);
    } while (strnlen(s, 255) > 0);

}


