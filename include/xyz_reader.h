// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.
//
// Copyright (c) Roman Ishchenko 2022.
#ifndef SCF_XYZ_READER_H
#define SCF_XYZ_READER_H

#include <stdio.h>
#include "atom.h"

#include "system.h"

/**
 * Reads atom types and coordinates from XYZ file and puts them to System structure
 * @param f pointer to a XYZ file
 * @param s pointer to a System structure
 */
int read_xyz(FILE *f, System *s);

#endif //SCF_XYZ_READER_H

