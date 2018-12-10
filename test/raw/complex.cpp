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

///////////////////
// Test invert() //
///////////////////
TEST_CASE("Complex matrix is inverted", "[invert_asymmetric]") {
  const int N = 100;
  const dcomplex diag_coeff = 0.75;
  const int offdiag_offset = 1;
  const dcomplex offdiag_coeff(0, 1.0);

  auto A = make_sparse_matrix<Complex>(N, diag_coeff, offdiag_offset, offdiag_coeff);

  auto invA = make_buffer<dcomplex>(N * N);
  invert(A.get(), invA.get(), N);

  auto id = make_buffer<dcomplex>(N * N);
  for(int i = 0; i < N; ++i) {
    for(int j = 0; j < N; ++j) {
      id[i + j*N] = i == j;
    }
  }

  auto prod = make_buffer<dcomplex>(N * N);

  // A * A^{-1}
  mm_product(A.get(), invA.get(), prod.get(), N);
  CHECK_THAT(prod.get(), IsCloseTo(id.get(), N * N));

  // A^{-1} * A
  mm_product(invA.get(), A.get(), prod.get(), N);
  CHECK_THAT(prod.get(), IsCloseTo(id.get(), N * N));
}

/////////////////////////////////////////
// Eigenproblems with complex matrices //
/////////////////////////////////////////

TEST_CASE("Complex eigenproblem is solved", "[worker_complex]") {

  using worker_t = arpack_worker<Complex, raw_storage>;
  using params_t = worker_t::params_t;

  const int N = 100;
  const dcomplex diag_coeff = 0.75;
  const int offdiag_offset = 1;
  const dcomplex offdiag_coeff(0, 1.0);
  const int nev = 10;

  auto spectrum_parts = {params_t::LargestMagnitude,
                         params_t::SmallestMagnitude,
                         params_t::LargestReal, params_t::SmallestReal,
                         params_t::LargestImag, params_t::SmallestImag};

  // Hermitian matrix A
  auto A = make_sparse_matrix<Complex>(N, diag_coeff, offdiag_offset, offdiag_coeff);
  // Inner product matrix
  auto M = make_inner_prod_matrix<Complex>(N);

  SECTION("Standard eigenproblem") {
    auto Aop = [&](dcomplex const* from, dcomplex * to) {
      mv_product(A.get(), from, to, N);
    };

    worker_t ar(N);

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.ncv = 30;
      ar(Aop, params);
      check_eigenvectors(ar, A.get(), N, nev);
    }
  }

  SECTION("Generalized eigenproblem: invert mode") {
    auto invM = make_buffer<dcomplex>(N * N);
    invert(M.get(), invM.get(), N);
    auto op_matrix = make_buffer<dcomplex>(N * N);
    mm_product(invM.get(), A.get(), op_matrix.get(), N);

    auto op = [&](dcomplex const* from, dcomplex * to) {
      mv_product(op_matrix.get(), from, to, N);
    };
    auto Bop = [&](dcomplex const* from, dcomplex * to) {
      mv_product(M.get(), from, to, N);
    };

    worker_t ar(N);

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.ncv = 50;
      ar(op, Bop, worker_t::Invert, params);
      check_eigenvectors(ar, A.get(), M.get(), N, nev);
    }
  }

  SECTION("Generalized eigenproblem: Shift-and-Invert mode") {
    dcomplex sigma(0.5, 0.5);

    auto AmM = make_buffer<dcomplex>(N * N);
    for(int i = 0; i < N; ++i) {
      for(int j = 0; j < N; ++j) {
        AmM[i + j*N] = A[i + j*N] - sigma * M[i + j*N];
      }
    }
    auto invAmM = make_buffer<dcomplex>(N * N);
    invert(AmM.get(), invAmM.get(), N);
    auto op_matrix = make_buffer<dcomplex>(N * N);
    mm_product(invAmM.get(), M.get(), op_matrix.get(), N);

    auto op = [&](dcomplex const* from, dcomplex * to) {
      mv_product(op_matrix.get(), from, to, N);
    };
    auto Bop = [&](dcomplex const* from, dcomplex * to) {
      mv_product(M.get(), from, to, N);
    };

    worker_t ar(N);

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.sigma = sigma;
      params.ncv = 50;
      ar(op, Bop, worker_t::ShiftAndInvert, params);
      check_eigenvectors(ar, A.get(), M.get(), N, nev);
    }
  }
}
