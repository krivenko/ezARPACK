/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2020 Igor Krivenko <igor.s.krivenko@gmail.com>
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

  using worker_t = arpack_worker<Complex, triqs_storage>;
  using params_t = worker_t::params_t;

  const int N = 100;
  const dcomplex diag_coeff_shift = -0.55;
  const dcomplex diag_coeff_amp = 1.0;
  const int offdiag_offset = 3;
  const dcomplex offdiag_coeff(0.5, 0.5);
  const int nev = 10;

  auto spectrum_parts = {
      params_t::LargestMagnitude, params_t::SmallestMagnitude,
      params_t::LargestReal,      params_t::SmallestReal,
      params_t::LargestImag,      params_t::SmallestImag};

  // Hermitian matrix A
  auto A = make_sparse_matrix<Complex>(N, diag_coeff_shift, diag_coeff_amp,
                                       offdiag_offset, offdiag_coeff);
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
      to = A * from;
    };

    worker_t ar(first_dim(A));

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(Aop, params);
      check_eigenvectors(ar, A);
      check_basis_vectors(ar);
    }
  }

  SECTION("Generalized eigenproblem: invert mode") {
    decltype(A) op_matrix = inverse(M) * A;

    auto op = [&](vector_const_view_t from, vector_view_t to) {
      to = op_matrix * from;
    };
    auto Bop = [&](vector_const_view_t from, vector_view_t to) {
      to = M * from;
    };

    worker_t ar(first_dim(A));

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(op, Bop, worker_t::Inverse, params);
      check_eigenvectors(ar, A, M);
      check_basis_vectors(ar, M);
    }
  }

  SECTION("Generalized eigenproblem: Shift-and-Invert mode") {
    dcomplex sigma(-1.0, 0.9);
    decltype(A) op_matrix = inverse(A - sigma * M) * M;

    auto op = [&](vector_const_view_t from, vector_view_t to) {
      to = op_matrix * from;
    };
    auto Bop = [&](vector_const_view_t from, vector_view_t to) {
      to = M * from;
    };

    worker_t ar(first_dim(A));

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.sigma = sigma;
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(op, Bop, worker_t::ShiftAndInvert, params);
      check_eigenvectors(ar, A, M);
      check_basis_vectors(ar, M);
    }
  }

  SECTION("Custom implementation of the Exact Shift Strategy") {
    std::vector<int> p;
    p.reserve(first_dim(A));
    auto shifts_f = [&](worker_t::complex_vector_const_view_t ritz_values,
                        worker_t::complex_vector_const_view_t ritz_bounds,
                        worker_t::complex_vector_view_t shifts) {
      int np = first_dim(shifts);
      if(np == 0) return;

      p.resize(np);
      std::iota(p.begin(), p.end(), 0);
      // p will contain the permutation that puts ritz_bounds in
      // the descending order of magnitude
      std::sort(p.begin(), p.end(), [&](int i, int j) {
        return std::abs(ritz_bounds(i)) > std::abs(ritz_bounds(j));
      });
      // Apply permutation p to ritz_values and use the result to fill shifts
      for(int i = 0; i < np; ++i)
        shifts(i) = ritz_values(p[i]);
    };

    SECTION("Standard eigenproblem") {
      auto Aop = [&](vector_const_view_t from, vector_view_t to) {
        to = A * from;
      };

      worker_t ar(first_dim(A));

      params_t params(nev, params_t::LargestMagnitude, params_t::Ritz);
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(Aop, params, shifts_f);
      check_eigenvectors(ar, A);
      check_basis_vectors(ar);
    }

    SECTION("Generalized eigenproblem: Shift-and-Invert mode") {
      dcomplex sigma(-1.0, 0.9);
      decltype(A) op_matrix = inverse(A - sigma * M) * M;

      auto op = [&](vector_const_view_t from, vector_view_t to) {
        to = op_matrix * from;
      };
      auto Bop = [&](vector_const_view_t from, vector_view_t to) {
        to = M * from;
      };

      worker_t ar(first_dim(A));

      params_t params(nev, params_t::LargestMagnitude, params_t::Ritz);
      params.sigma = sigma;
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(op, Bop, worker_t::ShiftAndInvert, params, shifts_f);
      check_eigenvectors(ar, A, M);
      check_basis_vectors(ar, M);
    }
  }
}
