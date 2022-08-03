/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2022 Igor Krivenko <igor.s.krivenko@gmail.com>
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

TEST_CASE("Symmetric eigenproblem is solved", "[solver_symmetric]") {

  using solver_t = arpack_solver<Symmetric, ublas_storage>;
  using params_t = solver_t::params_t;

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

  auto set_init_residual_vector = [](solver_t& ar) {
    for(int i = 0; i < N; ++i)
      ar.residual_vector()[i] = double(i) / N;
  };

  using vector_view_t = solver_t::vector_view_t;
  using vector_const_view_t = solver_t::vector_const_view_t;

  SECTION("Standard eigenproblem") {
    auto Aop = [&](vector_const_view_t in, vector_view_t out) {
      out = prod(A, in);
    };

    solver_t ar(A.size1());
    REQUIRE(ar.dim() == A.size1());

    for(auto e : spectrum_parts) {
      params_t params(nev, e, true);
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(Aop, params);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors(ar, A);
      check_basis_vectors(ar);
    }
  }

  SECTION("Generalized eigenproblem: invert mode") {
    decltype(A) invM = inverse(M);

    auto op = [&](vector_view_t in, vector_view_t out) {
      in = prod(A, in);
      out = prod(invM, in);
    };
    auto Bop = [&](vector_const_view_t in, vector_view_t out) {
      out = prod(M, in);
    };

    solver_t ar(A.size1());
    REQUIRE(ar.dim() == A.size1());

    for(auto e : spectrum_parts) {
      params_t params(nev, e, true);
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(op, Bop, solver_t::Inverse, params);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors(ar, A, M);
      check_basis_vectors(ar, M);
    }
  }

  SECTION("Generalized eigenproblem: Shift-and-Invert mode") {
    double sigma = 2.0;
    decltype(A) op_matrix = prod(inverse(A - sigma * M), M);

    auto op = [&](vector_view_t in, vector_view_t out) {
      out = prod(op_matrix, in);
    };
    auto Bop = [&](vector_const_view_t in, vector_view_t out) {
      out = prod(M, in);
    };

    solver_t ar(A.size1());
    REQUIRE(ar.dim() == A.size1());

    for(auto e : spectrum_parts) {
      params_t params(nev, e, true);
      params.sigma = sigma;
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(op, Bop, solver_t::ShiftAndInvert, params);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors(ar, A, M);
      check_basis_vectors(ar, M);
    }
  }

  SECTION("Generalized eigenproblem: Buckling mode") {
    double sigma = 3.3;
    decltype(A) op_matrix = prod(inverse(M - sigma * A), M);

    auto op = [&](vector_view_t in, vector_view_t out) {
      out = prod(op_matrix, in);
    };
    auto Bop = [&](vector_const_view_t in, vector_view_t out) {
      out = prod(M, in);
    };

    solver_t ar(A.size1());
    REQUIRE(ar.dim() == A.size1());

    for(auto e : spectrum_parts) {
      params_t params(nev, e, true);
      params.sigma = sigma;
      params.random_residual_vector = false;
      params.ncv = 30;
      set_init_residual_vector(ar);
      ar(op, Bop, solver_t::Buckling, params);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors(ar, M, A);
      check_basis_vectors(ar, M);
    }
  }

  SECTION("Generalized eigenproblem: Cayley transformed mode") {
    double sigma = 2.0;
    decltype(A) op_matrix = prod(inverse(A - sigma * M), A + sigma * M);

    auto op = [&](vector_view_t in, vector_view_t out) {
      out = prod(op_matrix, in);
    };
    auto Bop = [&](vector_const_view_t in, vector_view_t out) {
      out = prod(M, in);
    };

    solver_t ar(A.size1());
    REQUIRE(ar.dim() == A.size1());

    for(auto e : spectrum_parts) {
      params_t params(nev, e, true);
      params.sigma = sigma;
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(op, Bop, solver_t::Cayley, params);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors(ar, A, M);
      check_basis_vectors(ar, M);
    }
  }

  SECTION("Indirect access to workspace vectors") {
    solver_t ar(A.size1());
    REQUIRE(ar.dim() == A.size1());

    auto Aop = [&](vector_const_view_t, vector_view_t) {
      auto in = ar.workspace_vector(ar.in_vector_n());
      auto out = ar.workspace_vector(ar.out_vector_n());
      out = prod(A, in);
    };

    for(auto e : spectrum_parts) {
      params_t params(nev, e, true);
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(Aop, params);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors(ar, A);
      check_basis_vectors(ar);
    }

    CHECK_THROWS(ar.workspace_vector(-1));
    CHECK_THROWS(ar.workspace_vector(3));
  }

  SECTION("Custom implementation of the Exact Shift Strategy") {
    std::vector<int> p;
    p.reserve(A.size1());
    auto shifts_f = [&](solver_t::real_vector_const_view_t ritz_values,
                        solver_t::real_vector_const_view_t ritz_bounds,
                        solver_t::real_vector_view_t shifts) {
      int np = shifts.size();
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
        out = prod(A, in);
      };

      solver_t ar(A.size1());
      REQUIRE(ar.dim() == A.size1());

      params_t params(nev, params_t::LargestMagnitude, true);
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(Aop, params, shifts_f);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors(ar, A);
      check_basis_vectors(ar);
    }

    SECTION("Generalized eigenproblem: Shift-and-Invert mode") {
      double sigma = 2.0;
      decltype(A) op_matrix = prod(inverse(A - sigma * M), M);

      auto op = [&](vector_view_t in, vector_view_t out) {
        out = prod(op_matrix, in);
      };
      auto Bop = [&](vector_const_view_t in, vector_view_t out) {
        out = prod(M, in);
      };

      solver_t ar(A.size1());
      REQUIRE(ar.dim() == A.size1());

      params_t params(nev, params_t::LargestMagnitude, true);
      params.sigma = sigma;
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(op, Bop, solver_t::ShiftAndInvert, params, shifts_f);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors(ar, A, M);
      check_basis_vectors(ar, M);
    }
  }
}
