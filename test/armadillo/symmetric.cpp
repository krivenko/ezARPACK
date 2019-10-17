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

////////////////////////////////////////////////
// Eigenproblems with real symmetric matrices //
////////////////////////////////////////////////

TEST_CASE("Symmetric eigenproblem is solved", "[worker_symmetric]") {

  using worker_t = arpack_worker<Symmetric, armadillo_storage>;
  using params_t = worker_t::params_t;

  const int N = 100;
  const double diag_coeff = 1.5;
  const int offdiag_offset = 3;
  const double offdiag_coeff = 0.5;
  const int nev = 10;

  auto spectrum_parts = {params_t::Smallest,
                         params_t::Largest,
                         params_t::SmallestMagnitude,
                         params_t::LargestMagnitude,
                         params_t::BothEnds};

  // Symmetric matrix A
  auto A = make_sparse_matrix<Symmetric>(N,
                                         diag_coeff,
                                         offdiag_offset,
                                         offdiag_coeff);
  // Inner product matrix
  auto M = make_inner_prod_matrix<Symmetric>(N);

  auto set_init_residual_vector = [N](worker_t & ar) {
    for(int i = 0; i < N; ++i) ar.residual_vector()[i] = double(i) / N;
  };

  using vector_view_t = worker_t::vector_view_t;
  using vector_const_view_t = worker_t::vector_const_view_t;

  SECTION("Standard eigenproblem") {
    auto Aop = [&](vector_const_view_t from, vector_view_t to) {
      to = A * from;
    };

    worker_t ar(A.n_rows);

    for(auto e : spectrum_parts) {
      params_t params(nev, e, true);
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(Aop, params);
      check_eigenvectors(ar, A);
    }
  }

  SECTION("Generalized eigenproblem: invert mode") {
    decltype(A) invM = inv(M);

    auto op = [&](vector_view_t from, vector_view_t to) {
      from = A * from;
      to = invM * from;
    };
    auto Bop = [&](vector_const_view_t from, vector_view_t to) {
      to = M * from;
    };

    worker_t ar(A.n_rows);

    for(auto e : spectrum_parts) {
      params_t params(nev, e, true);
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(op, Bop, worker_t::Invert, params);
      check_eigenvectors(ar, A, M);
    }
  }

  SECTION("Generalized eigenproblem: Shift-and-Invert mode") {
    double sigma = 1.0;
    decltype(A) op_matrix = inv(A - sigma*M) * M;

    auto op = [&](vector_view_t from, vector_view_t to) {
      to = op_matrix * from;
    };
    auto Bop = [&](vector_const_view_t from, vector_view_t to) {
      to = M * from;
    };

    worker_t ar(A.n_rows);

    for(auto e : spectrum_parts) {
      params_t params(nev, e, true);
      params.sigma = sigma;
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(op, Bop, worker_t::ShiftAndInvert, params);
      check_eigenvectors(ar, A, M);
    }
  }

  SECTION("Generalized eigenproblem: Buckling mode") {
    double sigma = 1.0;
    decltype(A) op_matrix = inv(M - sigma*A) * M;

    auto op = [&](vector_view_t from, vector_view_t to) {
      to = op_matrix * from;
    };
    auto Bop = [&](vector_const_view_t from, vector_view_t to) {
      to = M * from;
    };

    worker_t ar(A.n_rows);

    for(auto e : spectrum_parts) {
      params_t params(nev, e, true);
      params.sigma = sigma;
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(op, Bop, worker_t::Buckling, params);
      check_eigenvectors(ar, M, A);
    }
  }

  SECTION("Generalized eigenproblem: Cayley transformed mode") {
    double sigma = 1.0;
    decltype(A) op_matrix = inv(A - sigma*M) * (A + sigma*M);

    auto op = [&](vector_view_t from, vector_view_t to) {
      to = op_matrix * from;
    };
    auto Bop = [&](vector_const_view_t from, vector_view_t to) {
      to = M * from;
    };

    worker_t ar(A.n_rows);

    for(auto e : spectrum_parts) {
      params_t params(nev, e, true);
      params.sigma = sigma;
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(op, Bop, worker_t::Cayley, params);
      check_eigenvectors(ar, A, M);
    }
  }
}
