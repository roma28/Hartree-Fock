// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.
//
// Copyright (c) Roman Ishchenko 2022.

#ifndef SCF_ATOM_H
#define SCF_ATOM_H

#include "atom_basis.h"

typedef struct {
    uint8_t atomic_number;
    gsl_vector *coords;
    atom_basis *basis;
    uint8_t Z;
} atom;

atom *atom_alloc();

void atom_free(atom *p);

#endif //SCF_ATOM_H
