// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.
//
// Copyright (c) Roman Ishchenko 2022.

#include "../include/xyz_reader.h"
#include "c-logger/logger.h"
#include "c-logger/loggerconf.h"


int read_xyz(FILE *f, System *s) {
    char buff[255];

    LOG_DEBUG("Reading xyz file");
    fgets(buff, 255, f);
    s->natoms = atoi(buff);
    if (s->natoms != 0) {
        s->atoms = malloc(s->natoms * sizeof(atom *));
        if (!s->atoms) {
            LOG_FATAL("Allocation error");
            return -1;
        }
        for (size_t i = 0; i < s->natoms; ++i) {
            fgets(buff, 255, f);
            uint32_t atom_number;
            double x, y, z;
            sscanf(buff, "%d\t%lf\t%lf\t%lf\n", &atom_number, &x, &y, &z);
            s->atoms[i] = atom_alloc();
            s->atoms[i]->Z = atom_number;
            s->atoms[i]->coords->data[0] = x;
            s->atoms[i]->coords->data[1] = y;
            s->atoms[i]->coords->data[2] = z;
        }
    }
    return 0;
}