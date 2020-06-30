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

////////////////////////////////////////////////
// Eigenproblems with real symmetric matrices //
////////////////////////////////////////////////

TEST_CASE("Symmetric eigenproblem is solved", "[worker_symmetric]") {

  using worker_t = arpack_worker<Symmetric, triqs_storage>;
  using params_t = worker_t::params_t;

  const int N = 100;
  const double diag_coeff_shift = -0.55;
  const double diag_coeff_amp = 1.0;
  const int offdiag_offset = 3;
  const double offdiag_coeff = -1.05;
  const int nev = 10;

  auto spectrum_parts = {params_t::Smallest, params_t::Largest,
                         params_t::SmallestMagnitude,
                         params_t::LargestMagnitude, params_t::BothEnds};

  // Symmetric matrix A
  auto A = make_sparse_matrix<Symmetric>(N, diag_coeff_shift, diag_coeff_amp,
                                         offdiag_offset, offdiag_coeff);
  // Inner product matrix
  auto M = make_inner_prod_matrix<Symmetric>(N);

  auto set_init_residual_vector = [](worker_t& ar) {
    for(int i = 0; i < N; ++i)
      ar.residual_vector()[i] = double(i) / N;
  };

  using vector_view_t = worker_t::vector_view_t;
  using vector_const_view_t = worker_t::vector_const_view_t;

  SECTION("Standard eigenproblem") {
    auto Aop = [&](vector_const_view_t in, vector_view_t out) {
      out = A * in;
    };

    worker_t ar(first_dim(A));

    for(auto e : spectrum_parts) {
      params_t params(nev, e, true);
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(Aop, params);
      check_eigenvectors(ar, A);
      check_basis_vectors(ar);
    }
  }

  SECTION("Generalized eigenproblem: invert mode") {
    decltype(A) invM = inverse(M);

    auto op = [&](vector_view_t in, vector_view_t out) {
      in = A * in;
      out = invM * in;
    };
    auto Bop = [&](vector_const_view_t in, vector_view_t out) {
      out = M * in;
    };

    worker_t ar(first_dim(A));

    for(auto e : spectrum_parts) {
      params_t params(nev, e, true);
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(op, Bop, worker_t::Inverse, params);
      check_eigenvectors(ar, A, M);
      check_basis_vectors(ar, M);
    }
  }

  SECTION("Generalized eigenproblem: Shift-and-Invert mode") {
    double sigma = 2.0;
    decltype(A) op_matrix = inverse(A - sigma * M) * M;

    auto op = [&](vector_view_t in, vector_view_t out) {
      out = op_matrix * in;
    };
    auto Bop = [&](vector_const_view_t in, vector_view_t out) {
      out = M * in;
    };

    worker_t ar(first_dim(A));

    for(auto e : spectrum_parts) {
      params_t params(nev, e, true);
      params.sigma = sigma;
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(op, Bop, worker_t::ShiftAndInvert, params);
      check_eigenvectors(ar, A, M);
      check_basis_vectors(ar, M);
    }
  }

  SECTION("Generalized eigenproblem: Buckling mode") {
    double sigma = 3.3;
    decltype(A) op_matrix = inverse(M - sigma * A) * M;

    auto op = [&](vector_view_t in, vector_view_t out) {
      out = op_matrix * in;
    };
    auto Bop = [&](vector_const_view_t in, vector_view_t out) {
      out = M * in;
    };

    worker_t ar(first_dim(A));

    for(auto e : spectrum_parts) {
      params_t params(nev, e, true);
      params.sigma = sigma;
      params.random_residual_vector = false;
      params.ncv = 30;
      set_init_residual_vector(ar);
      ar(op, Bop, worker_t::Buckling, params);
      check_eigenvectors(ar, M, A);
      check_basis_vectors(ar, M);
    }
  }

  SECTION("Generalized eigenproblem: Cayley transformed mode") {
    double sigma = 2.0;
    decltype(A) op_matrix = inverse(A - sigma * M) * (A + sigma * M);

    auto op = [&](vector_view_t in, vector_view_t out) {
      out = op_matrix * in;
    };
    auto Bop = [&](vector_const_view_t in, vector_view_t out) {
      out = M * in;
    };

    worker_t ar(first_dim(A));

    for(auto e : spectrum_parts) {
      params_t params(nev, e, true);
      params.sigma = sigma;
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(op, Bop, worker_t::Cayley, params);
      check_eigenvectors(ar, A, M);
      check_basis_vectors(ar, M);
    }
  }

  SECTION("Custom implementation of the Exact Shift Strategy") {
    std::vector<int> p;
    p.reserve(first_dim(A));
    auto shifts_f = [&](worker_t::real_vector_const_view_t ritz_values,
                        worker_t::real_vector_const_view_t ritz_bounds,
                        worker_t::real_vector_view_t shifts) {
      int np = first_dim(shifts);
      if(np == 0) return;

      p.resize(np);
      std::iota(p.begin(), p.end(), 0);
      // p will contain the permutation that puts ritz_bounds in
      // the descending order of magnitude
      std::sort(p.begin(), p.end(),
                [&](int i, int j) { return ritz_bounds(i) > ritz_bounds(j); });
      // Apply permutation p to ritz_values and use the result to fill shifts
      for(int i = 0; i < np; ++i)
        shifts(i) = ritz_values(p[i]);
    };

    SECTION("Standard eigenproblem") {
      auto Aop = [&](vector_const_view_t in, vector_view_t out) {
        out = A * in;
      };

      worker_t ar(first_dim(A));

      params_t params(nev, params_t::LargestMagnitude, true);
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(Aop, params, shifts_f);
      check_eigenvectors(ar, A);
      check_basis_vectors(ar);
    }

    SECTION("Generalized eigenproblem: Shift-and-Invert mode") {
      double sigma = 2.0;
      decltype(A) op_matrix = inverse(A - sigma * M) * M;

      auto op = [&](vector_view_t in, vector_view_t out) {
        out = op_matrix * in;
      };
      auto Bop = [&](vector_const_view_t in, vector_view_t out) {
        out = M * in;
      };

      worker_t ar(first_dim(A));

      params_t params(nev, params_t::LargestMagnitude, true);
      params.sigma = sigma;
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(op, Bop, worker_t::ShiftAndInvert, params, shifts_f);
      check_eigenvectors(ar, A, M);
      check_basis_vectors(ar, M);
    }
  }
}
