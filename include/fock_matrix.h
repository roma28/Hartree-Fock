// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.
//
// Copyright (c) Roman Ishchenko 2022.

#ifndef SCF_FOCK_MATRIX_H
#define SCF_FOCK_MATRIX_H

#include <gsl/gsl_matrix.h>
#include "atom.h"

gsl_matrix *make_fock_matrix(const atom **atoms, const uint32_t natoms);

#endif //SCF_FOCK_MATRIX_H
