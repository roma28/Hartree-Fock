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
#include <gsl/gsl_math.h>


double overlap_integral(const basis_function *a, const basis_function *b) {
    assert((a->L == 0) && (b->L == 0)); //implemented only for S functions

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

double kinetic_energy_integral(const basis_function *a, const basis_function *b) {
    assert((a->L == 0) && (b->L == 0)); //implemented only for S functions

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
}

double
nuclear_repulsion_integral(const basis_function *a, const basis_function *b, const uint8_t *Z,
                           gsl_vector **atom_origins, size_t natoms) {
    assert((a->L == 0) && (b->L == 0)); //implemented only for S functions

    double V = 0.0;
    double dR2 = deltaR2(a->origin, b->origin);

    for (size_t i = 0; i < a->n_primitives; ++i) {
        for (size_t j = 0; j < b->n_primitives; ++j) {

            gsl_vector *p = gaussian_center(a->origin, b->origin, a->exponents->data[i], b->exponents->data[j]);
            double theta = M_2_SQRTPI * sqrt(a->exponents->data[i] + b->exponents->data[j]);

            for (size_t n = 0; n < natoms; n++) {

                uint8_t z = Z[n];
                const gsl_vector *c = atom_origins[n];
                double t = (a->exponents->data[i] + b->exponents->data[j]) * deltaR2(p, c);

                V -= theta * a->contractions->data[i] * b->contractions->data[j] *
                     v00(a->exponents->data[i], b->exponents->data[j], dR2, z, t);
            }
            gsl_vector_free(p);
        }
    }

    return V;
}

double electron_repulsion_integral(const basis_function *a, const basis_function *b, const basis_function *x,
                                   const basis_function *y) {
    assert((a->L == 0) && (b->L == 0) && ((x->L == 0)) && ((y->L == 0))); //implemented only for S functions
    double G = 0.0;
    double dab2 = deltaR2(a->origin, b->origin);
    double dxy2 = deltaR2(x->origin, y->origin);

    for (size_t i = 0; i < a->n_primitives; ++i) {
        for (size_t j = 0; j < b->n_primitives; ++j) {

            gsl_vector *p = gaussian_center(a->origin, b->origin, a->exponents->data[i], b->exponents->data[j]);
            double s1 = a->exponents->data[i] * b->exponents->data[j];

            for (size_t k = 0; k < x->n_primitives; ++k) {
                for (size_t l = 0; l < y->n_primitives; ++l) {

                    double s2 = x->exponents->data[k] * y->exponents->data[l];
                    gsl_vector *q = gaussian_center(x->origin, y->origin, a->exponents->data[k], b->exponents->data[l]);
                    double dpq2 = deltaR2(p, q);
                    double t = s1 * s2 * dpq2 / (s1 + s2);
                    double lambda = M_2_SQRTPI * s1 * s2 / (s1 + s2);

                    G += a->contractions->data[i] * b->contractions->data[j] *
                         x->contractions->data[k] * y->contractions->data[l] *
                         lambda * g0000(a->exponents->data[i], b->exponents->data[j], x->exponents->data[k],
                                        y->exponents->data[l],
                                        dab2, dxy2, t);

                    gsl_vector_free(q);
                }
            }
            gsl_vector_free(p);
        }
    }

    return G;
}

