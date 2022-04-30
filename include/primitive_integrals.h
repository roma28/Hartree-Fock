// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.
//
// Copyright (c) Roman Ishchenko 2022.

//
// Created by Roman Ishchenko on 30.04.2022.
//

#ifndef SCF_PRIMITIVE_INTEGRALS_H
#define SCF_PRIMITIVE_INTEGRALS_H
#include "atom_basis.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>

double N(double alpha);

double Q(const gsl_vector *Ra, const gsl_vector *Rb);

double s00(const double exp_a, const double exp_b, const gsl_vector *Ra, const gsl_vector *Rb);

#endif //SCF_PRIMITIVE_INTEGRALS_H
