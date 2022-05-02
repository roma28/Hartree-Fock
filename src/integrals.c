// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.
//
// Copyright (c) Roman Ishchenko 2022.

#include "../include/integrals.h"
#include <assert.h>


double S(const basis_function *a, const basis_function *b) {
    assert((a->L == 0) && (b->L == 0));

    double S = 0.0;
    double dR2 = deltaR2(a->origin, b->origin);

    for (size_t i = 0; i < a->n_primitives; ++i) {
        for (size_t j = 0; j < b->n_primitives; ++j) {
            S += a->contractions->data[i] * b->contractions->data[j] *
                 s00(a->exponents->data[i], b->exponents->data[j], dR2);
        }
    }

    return S;
}

double T(const basis_function *a, const basis_function *b) {
    assert((a->L == 0) && (b->L == 0));

    double T = 0.0;
    double dR2 = deltaR2(a->origin, b->origin);

    for (size_t i = 0; i < a->n_primitives; ++i) {
        for (size_t j = 0; j < b->n_primitives; ++j) {
            T += a->contractions->data[i] * b->contractions->data[j] *
                 s00(a->exponents->data[i], b->exponents->data[j], dR2) *
                 k00(a->exponents->data[i], b->exponents->data[j], dR2);
        }
    }

    return T;

};

