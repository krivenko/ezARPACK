/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2018 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/

#include "common.hpp"

////////////////////////////////////////////////
// Eigenproblems with real symmetric matrices //
////////////////////////////////////////////////

using params_t = arpack_worker<Symmetric>::params_t;

const int N = 100;
const double diag_coeff = 0.75;
const int offdiag_offset = 3;
const double offdiag_coeff = 1.0;
const int nev = 10;

// Symmetric matrix A
auto A = make_sparse_matrix<Symmetric>(N, diag_coeff, offdiag_offset, offdiag_coeff);
// Inner product matrix
auto M = make_inner_prod_matrix<Symmetric>(N);

auto spectrum_parts = {params_t::Smallest, params_t::Largest,
                       params_t::SmallestMagnitude, params_t::LargestMagnitude,
                       params_t::BothEnds};

TEST(arpack_worker_symmetric, InnerProduct) {
 ASSERT_GT(eigenvalues(M)(0),.0);
}

// Standard eigenproblem
TEST(arpack_worker_symmetric, Standard) {
 auto Aop = [](vector_const_view<double> from, int, vector_view<double> to, int) {
  to = A*from;
 };

 arpack_worker<Symmetric> ar(first_dim(A));

 for(auto e : spectrum_parts) {
  params_t params(nev, e, true);
  ar(Aop, params);
  check_eigenvectors(ar,A);
 }
}

// Generalized eigenproblem: invert mode
TEST(arpack_worker_symmetric, Invert) {
 decltype(A) invM = inverse(M);

 auto op = [&invM](vector_view<double> from, int, vector_view<double> to, int, bool) {
  from = A * from;
  to = invM * from;
 };
 auto Bop = [](vector_const_view<double> from, int, vector_view<double> to, int) {
  to = M * from;
 };

 arpack_worker<Symmetric> ar(first_dim(A));

 for(auto e : spectrum_parts) {
  params_t params(nev, e, true);
  ar(op, Bop, arpack_worker<Symmetric>::Invert, params);
  check_eigenvectors(ar,A,M);
 }
}

// Generalized eigenproblem: Shift-and-Invert mode
TEST(arpack_worker_symmetric, ShiftAndInvert) {
 double sigma = 1.0;
 decltype(A) inv = inverse(A - sigma*M) * M;

 auto op = [&inv](vector_view<double> from, int, vector_view<double> to, int, bool) {
  to = inv * from;
 };
 auto Bop = [](vector_const_view<double> from, int, vector_view<double> to, int) {
  to = M * from;
 };

 arpack_worker<Symmetric> ar(first_dim(A));

 for(auto e : spectrum_parts) {
  params_t params(nev, e, true);
  params.sigma = sigma;
  ar(op, Bop, arpack_worker<Symmetric>::ShiftAndInvert, params);
  check_eigenvectors(ar,A,M);
 }
}

// Generalized eigenproblem: Buckling mode
TEST(arpack_worker_symmetric, Buckling) {
 double sigma = 1.0;
 decltype(A) inv = inverse(A - sigma*M) * A;

 auto op = [&inv](vector_view<double> from, int, vector_view<double> to, int, bool) {
  to = inv * from;
 };
 auto Bop = [](vector_const_view<double> from, int, vector_view<double> to, int) {
  to = M * from;
 };

 arpack_worker<Symmetric> ar(first_dim(A));

 for(auto e : spectrum_parts) {
  params_t params(nev, e, true);
  params.sigma = sigma;
  ar(op, Bop, arpack_worker<Symmetric>::Buckling, params);
  check_eigenvectors(ar,A,M);
 }
}

// Generalized eigenproblem: Cayley transformed mode
TEST(arpack_worker_symmetric, Cayley) {
 double sigma = 1.0;
 decltype(A) inv = inverse(A - sigma*M) * (A + sigma*M);

 auto op = [&inv](vector_view<double> from, int, vector_view<double> to, int, bool) {
  to = inv * from;
 };
 auto Bop = [&](vector_const_view<double> from, int, vector_view<double> to, int) {
  to = M * from;
 };

 arpack_worker<Symmetric> ar(first_dim(A));

 for(auto e : spectrum_parts) {
  params_t params(nev, e, true);
  params.sigma = sigma;
  ar(op, Bop, arpack_worker<Symmetric>::Cayley, params);
  check_eigenvectors(ar,A,M);
 }
}

MAKE_MAIN;
