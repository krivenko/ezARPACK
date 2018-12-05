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

/////////////////////////////////////////
// Eigenproblems with complex matrices //
/////////////////////////////////////////

using params_t = arpack_worker<Complex>::params_t;

auto spectrum_parts = {params_t::LargestMagnitude,
                       params_t::SmallestMagnitude,
                       params_t::LargestReal, params_t::SmallestReal,
                       params_t::LargestImag, params_t::SmallestImag};

const int N = 100;
const double diag_coeff = 0.75;
const int offdiag_offset = 1;
const dcomplex offdiag_coeff = 1_j;
const int nev = 10;

// Hermitian matrix A
auto A = make_sparse_matrix<Complex>(N, diag_coeff, offdiag_offset, offdiag_coeff);
// Inner product matrix
auto M = make_inner_prod_matrix<Complex>(N);

TEST(arpack_worker_symmetric, InnerProduct) {
 ASSERT_GT(eigenvalues(M)(0),.0);
}

// Standard eigenproblem
TEST(arpack_worker_complex, Standard) {
 auto Aop = [](vector_const_view<dcomplex> from, vector_view<dcomplex> to) {
  to = A*from;
 };

 arpack_worker<Complex> ar(first_dim(A));

 for(auto e : spectrum_parts) {
  params_t params(nev, e, params_t::Ritz);
  params.ncv = 30;
  ar(Aop, params);
  check_eigenvectors(ar,A);
 }
}

// Generalized eigenproblem: invert mode
TEST(arpack_worker_complex, Invert) {
 decltype(A) invMA = inverse(M) * A;

 auto op = [&invMA](vector_const_view<dcomplex> from, vector_view<dcomplex> to) {
  to = invMA * from;
 };
 auto Bop = [](vector_const_view<dcomplex> from, vector_view<dcomplex> to) {
  to = M * from;
 };

 arpack_worker<Complex> ar(first_dim(A));

 for(auto e : spectrum_parts) {
  params_t params(nev, e, params_t::Ritz);
  params.ncv = 50;
  ar(op, Bop, arpack_worker<Complex>::Invert, params);
  check_eigenvectors(ar,A,M);
 }
}

// Generalized eigenproblem: Shift-and-Invert mode
TEST(arpack_worker_complex, ShiftAndInvert) {
 dcomplex sigma = 0.5 + 0.5_j;
 decltype(A) inv = inverse(A - sigma*M) * M;

 auto op = [&inv](vector_const_view<dcomplex> from, vector_view<dcomplex> to) {
  to = inv * from;
 };
 auto Bop = [](vector_const_view<dcomplex> from, vector_view<dcomplex> to) {
  to = M * from;
 };

 arpack_worker<Complex> ar(first_dim(A));

 for(auto e : spectrum_parts) {
  params_t params(nev, e, params_t::Ritz);
  params.sigma = sigma;
  params.ncv = 50;
  ar(op, Bop, arpack_worker<Complex>::ShiftAndInvert, params);
  check_eigenvectors(ar,A,M);
 }
}

MAKE_MAIN;
