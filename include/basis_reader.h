// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.
//
// Copyright (c) Roman Ishchenko 2022.


#ifndef SCF_BASIS_READER_H
#define SCF_BASIS_READER_H
#include <stdio.h>
#include "atom_basis.h"
#include <string.h>

void parse_basis_file(FILE *f);

#endif //SCF_BASIS_READER_H