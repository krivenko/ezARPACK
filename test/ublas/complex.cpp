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

/////////////////////////////////////////
// Eigenproblems with complex matrices //
/////////////////////////////////////////

TEST_CASE("Complex eigenproblem is solved", "[worker_complex]") {

  using worker_t = arpack_worker<Complex, ublas_storage>;
  using params_t = worker_t::params_t;

  const int N = 100;
  const dcomplex diag_coeff = 1.5;
  const int offdiag_offset = 1;
  const dcomplex offdiag_coeff(0, 0.1);
  const int nev = 10;

  auto spectrum_parts = {
      params_t::LargestMagnitude, params_t::SmallestMagnitude,
      params_t::LargestReal,      params_t::SmallestReal,
      params_t::LargestImag,      params_t::SmallestImag};

  // Hermitian matrix A
  auto A =
      make_sparse_matrix<Complex>(N, diag_coeff, offdiag_offset, offdiag_coeff);
  // Inner product matrix
  auto M = make_inner_prod_matrix<Complex>(N);

  auto set_init_residual_vector = [](worker_t& ar) {
    for(int i = 0; i < N; ++i)
      ar.residual_vector()[i] = double(i) / N;
  };

  using vector_view_t = worker_t::vector_view_t;
  using vector_const_view_t = worker_t::vector_const_view_t;

  SECTION("Standard eigenproblem") {
    auto Aop = [&](vector_const_view_t from, vector_view_t to) {
      to = prod(A, from);
    };

    worker_t ar(A.size1());

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.ncv = 50;
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(Aop, params);
      check_eigenvectors(ar, A);
    }
  }

  SECTION("Generalized eigenproblem: invert mode") {
    decltype(A) op_matrix = prod(inverse(M), A);

    auto op = [&](vector_const_view_t from, vector_view_t to) {
      to = prod(op_matrix, from);
    };
    auto Bop = [&](vector_const_view_t from, vector_view_t to) {
      to = prod(M, from);
    };

    worker_t ar(A.size1());

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.ncv = 50;
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(op, Bop, worker_t::Invert, params);
      check_eigenvectors(ar, A, M);
    }
  }

  SECTION("Generalized eigenproblem: Shift-and-Invert mode") {
    dcomplex sigma(0.5, 0.5);
    decltype(A) op_matrix = prod(inverse(A - sigma * M), M);

    auto op = [&](vector_const_view_t from, vector_view_t to) {
      to = prod(op_matrix, from);
    };
    auto Bop = [&](vector_const_view_t from, vector_view_t to) {
      to = prod(M, from);
    };

    worker_t ar(A.size1());

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.sigma = sigma;
      params.ncv = 50;
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(op, Bop, worker_t::ShiftAndInvert, params);
      check_eigenvectors(ar, A, M);
    }
  }
}
