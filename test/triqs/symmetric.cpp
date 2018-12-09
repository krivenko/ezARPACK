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

TEST_CASE("Symmetric eigenproblem is solved", "[worker_symmetric]") {

  using worker_t = arpack_worker<Symmetric, triqs_storage>;
  using params_t = worker_t::params_t;

  const int N = 100;
  const double diag_coeff = 0.5;
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

  SECTION("Standard eigenproblem") {
    auto Aop = [&](vector_const_view<double> from, vector_view<double> to) {
      to = A * from;
    };

    worker_t ar(first_dim(A));

    for(auto e : spectrum_parts) {
      params_t params(nev, e, true);
      ar(Aop, params);
      check_eigenvectors(ar, A);
    }
  }

  SECTION("Generalized eigenproblem: invert mode") {
    decltype(A) invM = inverse(M);

    auto op = [&](vector_view<double> from, vector_view<double> to) {
      from = A * from;
      to = invM * from;
    };
    auto Bop = [&](vector_const_view<double> from, vector_view<double> to) {
      to = M * from;
    };

    worker_t ar(first_dim(A));

    for(auto e : spectrum_parts) {
      params_t params(nev, e, true);
      ar(op, Bop, worker_t::Invert, params);
      check_eigenvectors(ar, A, M);
    }
  }

  SECTION("Generalized eigenproblem: Shift-and-Invert mode") {
    double sigma = 1.0;
    decltype(A) op_matrix = inverse(A - sigma*M) * M;

    auto op = [&](vector_view<double> from, vector_view<double> to) {
      to = op_matrix * from;
    };
    auto Bop = [&](vector_const_view<double> from, vector_view<double> to) {
      to = M * from;
    };

    worker_t ar(first_dim(A));

    for(auto e : spectrum_parts) {
      params_t params(nev, e, true);
      params.sigma = sigma;
      ar(op, Bop, worker_t::ShiftAndInvert, params);
      check_eigenvectors(ar, A, M);
    }
  }

  SECTION("Generalized eigenproblem: Buckling mode") {
    double sigma = 1.0;
    decltype(A) op_matrix = inverse(M - sigma*A) * M;

    auto op = [&](vector_view<double> from, vector_view<double> to) {
      to = op_matrix * from;
    };
    auto Bop = [&](vector_const_view<double> from, vector_view<double> to) {
      to = M * from;
    };

    worker_t ar(first_dim(A));

    for(auto e : spectrum_parts) {
      params_t params(nev, e, true);
      params.sigma = sigma;
      ar(op, Bop, worker_t::Buckling, params);
      check_eigenvectors(ar, M, A);
    }
  }

  SECTION("Generalized eigenproblem: Cayley transformed mode") {
    double sigma = 1.0;
    decltype(A) op_matrix = inverse(A - sigma*M) * (A + sigma*M);

    auto op = [&](vector_view<double> from, vector_view<double> to) {
      to = op_matrix * from;
    };
    auto Bop = [&](vector_const_view<double> from, vector_view<double> to) {
      to = M * from;
    };

    worker_t ar(first_dim(A));

    for(auto e : spectrum_parts) {
      params_t params(nev, e, true);
      params.sigma = sigma;
      ar(op, Bop, worker_t::Cayley, params);
      check_eigenvectors(ar, A, M);
    }
  }
}
