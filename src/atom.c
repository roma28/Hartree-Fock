// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.
//
// Copyright (c) Roman Ishchenko 2022.


#include "../include/atom.h"

atom *atom_alloc() {
    atom *a = malloc(sizeof(atom));
    if (a) a->coords = gsl_vector_alloc(3);
    return a;
}

void atom_dealloc(atom *p) {
    atom_basis_dealloc(p->basis);
    gsl_vector_free(p->coords);
    free(p);
    p = NULL;
}