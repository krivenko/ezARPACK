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

//////////////////////////////////////////////
// Eigenproblems with general real matrices //
//////////////////////////////////////////////

TEST_CASE("Asymmetric eigenproblem is solved", "[worker_asymmetric]") {

  using worker_t = arpack_worker<ezarpack::Asymmetric, eigen_storage>;
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
  auto A = make_sparse_matrix<ezarpack::Asymmetric>(N, diag_coeff, offdiag_offset, offdiag_coeff);
  // Inner product matrix
  auto M = make_inner_prod_matrix<ezarpack::Asymmetric>(N);

  using vector_view_t = worker_t::vector_view_t;
  using vector_const_view_t = worker_t::vector_const_view_t;

  SECTION("Standard eigenproblem") {
    auto Aop = [&](vector_const_view_t from, vector_view_t to) {
      to = A * from;
    };

    worker_t ar(A.rows());

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      ar(Aop, params);
      check_eigenvectors(ar, A);
    }
  }

  SECTION("Generalized eigenproblem: invert mode") {
    decltype(A) op_matrix = M.inverse() * A;

    auto op = [&](vector_const_view_t from, vector_view_t to) {
      to = op_matrix * from;
    };
    auto Bop = [&](vector_const_view_t from, vector_view_t to) {
      to = M * from;
    };

    worker_t ar(A.rows());

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.ncv = 30;
      ar(op, Bop, worker_t::Invert, params);
      check_eigenvectors(ar, A, M);
    }
  }
}
