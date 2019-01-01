/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2019 Igor Krivenko <igor.s.krivenko@gmail.com>
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
TEST_CASE("Asymmetric matrix is inverted", "[invert_asymmetric]") {
  const int N = 100;
  const double diag_coeff = 0.75;
  const int offdiag_offset = 3;
  const double offdiag_coeff = 1.0;

  auto A = make_sparse_matrix<Asymmetric>(N, diag_coeff, offdiag_offset, offdiag_coeff);

  auto invA = make_buffer<double>(N * N);
  invert(A.get(), invA.get(), N);

  auto id = make_buffer<double>(N * N);
  for(int i = 0; i < N; ++i) {
    for(int j = 0; j < N; ++j) {
      id[i + j*N] = i == j;
    }
  }

  auto prod = make_buffer<double>(N * N);

  // A * A^{-1}
  mm_product(A.get(), invA.get(), prod.get(), N);
  CHECK_THAT(prod.get(), IsCloseTo(id.get(), N * N));

  // A^{-1} * A
  mm_product(invA.get(), A.get(), prod.get(), N);
  CHECK_THAT(prod.get(), IsCloseTo(id.get(), N * N));
}

//////////////////////////////////////////////
// Eigenproblems with general real matrices //
//////////////////////////////////////////////

TEST_CASE("Asymmetric eigenproblem is solved", "[worker_asymmetric]") {

  using worker_t = arpack_worker<Asymmetric, raw_storage>;
  using params_t = worker_t::params_t;

  const int N = 100;
  const double diag_coeff = 0.75;
  const int offdiag_offset = 3;
  const double offdiag_coeff = 1.0;
  const int nev = 10;

  auto spectrum_parts = {params_t::LargestMagnitude,
                         params_t::SmallestMagnitude,
                         params_t::LargestReal, params_t::SmallestReal,
                         params_t::LargestImag, params_t::SmallestImag};

  // Asymmetric matrix A
  auto A = make_sparse_matrix<Asymmetric>(N, diag_coeff, offdiag_offset, offdiag_coeff);
  // Inner product matrix
  auto M = make_inner_prod_matrix<Asymmetric>(N);

  using vector_view_t = worker_t::vector_view_t;
  using vector_const_view_t = worker_t::vector_const_view_t;

  SECTION("Standard eigenproblem") {
    auto Aop = [&](vector_const_view_t from, vector_view_t to) {
      mv_product(A.get(), from, to, N);
    };

    worker_t ar(N);

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      ar(Aop, params);
      check_eigenvectors(ar, A.get(), N, nev);
    }
  }

  SECTION("Generalized eigenproblem: invert mode") {
    auto invM = make_buffer<double>(N * N);
    invert(M.get(), invM.get(), N);
    auto op_matrix = make_buffer<double>(N * N);
    mm_product(invM.get(), A.get(), op_matrix.get(), N);

    auto op = [&](vector_const_view_t from, vector_view_t to) {
      mv_product(op_matrix.get(), from, to, N);
    };
    auto Bop = [&](vector_const_view_t from, vector_view_t to) {
      mv_product(M.get(), from, to, N);
    };

    worker_t ar(N);

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.ncv = 30;
      ar(op, Bop, worker_t::Invert, params);
      check_eigenvectors(ar, A.get(), M.get(), N, nev);
    }
  }
}
