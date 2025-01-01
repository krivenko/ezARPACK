/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2025 Igor Krivenko
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

  using solver_t = arpack_solver<ezarpack::Symmetric, ublas_storage>;

  const int N = 100;
  const double diag_coeff_mean = 1.0;
  const int offdiag_offset = 3;
  const double offdiag_coeff_mean = -0.1;
  const double offdiag_coeff_diff = 0;
  const int nev = 8;

  const double sigma = 0.102;

  // Symmetric matrix A
  auto A = make_sparse_matrix<ezarpack::Symmetric>(
      N, diag_coeff_mean, offdiag_offset, offdiag_coeff_mean,
      offdiag_coeff_diff);
  // Inner product matrix
  auto M = make_inner_prod_matrix<ezarpack::Symmetric>(N);

  // Testing helper
  auto testing = make_testing_helper<solver_t>(A, M, N, nev);

  using vv_t = solver_t::vector_view_t;
  using vcv_t = solver_t::vector_const_view_t;

  SECTION("Standard eigenproblem") {
    auto Aop = [&](vcv_t in, vv_t out) { out = prod(A, in); };

    solver_t ar(A.size1());
    testing.standard_eigenproblems(ar, Aop);
  }

  SECTION("Generalized eigenproblem: invert mode") {
    decltype(A) invM = inverse(M);

    auto op = [&](vv_t in, vv_t out) {
      in = prod(A, in);
      out = prod(invM, in);
    };
    auto Bop = [&](vcv_t in, vv_t out) { out = prod(M, in); };

    solver_t ar(A.size1());
    testing.generalized_eigenproblems(ar, solver_t::Inverse, op, Bop);
  }

  SECTION("Generalized eigenproblem: Shift-and-Invert mode") {
    decltype(A) op_mat = prod(inverse(A - sigma * M), M);

    auto op = [&](vv_t in, vv_t out) { out = prod(op_mat, in); };
    auto Bop = [&](vcv_t in, vv_t out) { out = prod(M, in); };

    solver_t ar(A.size1());
    testing.generalized_eigenproblems(ar, solver_t::ShiftAndInvert, op, Bop,
                                      sigma);
  }

  SECTION("Generalized eigenproblem: Buckling mode") {
    decltype(A) op_mat = prod(inverse(M - sigma * A), M);

    auto op = [&](vv_t in, vv_t out) { out = prod(op_mat, in); };
    auto Bop = [&](vcv_t in, vv_t out) { out = prod(M, in); };

    solver_t ar(A.size1());
    const int ncv = 30;
    testing.generalized_eigenproblems(ar, solver_t::Buckling, op, Bop, sigma,
                                      ncv);
  }

  SECTION("Generalized eigenproblem: Cayley transformed mode") {
    decltype(A) op_mat = prod(inverse(A - sigma * M), A + sigma * M);

    auto op = [&](vv_t in, vv_t out) { out = prod(op_mat, in); };
    auto Bop = [&](vcv_t in, vv_t out) { out = prod(M, in); };

    solver_t ar(A.size1());
    testing.generalized_eigenproblems(ar, solver_t::Cayley, op, Bop, sigma);
  }

  SECTION("Indirect access to workspace vectors") {
    solver_t ar(A.size1());

    auto Aop = [&](vcv_t, vv_t) {
      auto in = ar.workspace_vector(ar.in_vector_n());
      auto out = ar.workspace_vector(ar.out_vector_n());
      out = prod(A, in);
    };

    testing.standard_eigenproblems(ar, Aop);

    CHECK_THROWS(ar.workspace_vector(-1));
    CHECK_THROWS(ar.workspace_vector(3));
  }

  SECTION("Skip computation of eigenvectors") {
    solver_t ar(A.size1());

    SECTION("Standard eigenproblem") {
      auto Aop = [&](vcv_t in, vv_t out) { out = prod(A, in); };

      testing.standard_skip_eigenvectors(ar, Aop);
    }

    SECTION("Generalized eigenproblem: invert mode") {
      decltype(A) invM = inverse(M);

      auto op = [&](vv_t in, vv_t out) {
        in = prod(A, in);
        out = prod(invM, in);
      };
      auto Bop = [&](vcv_t in, vv_t out) { out = prod(M, in); };

      testing.generalized_skip_eigenvectors(ar, op, Bop);
    }
  }

  SECTION("Custom implementation of the Exact Shift Strategy") {
    using rvcv_t = solver_t::real_vector_const_view_t;
    using rvv_t = solver_t::real_vector_view_t;
    auto size_f = [](rvv_t shifts) -> int { return shifts.size(); };
    exact_shift_strategy<ezarpack::Symmetric, rvcv_t, rvv_t> shifts_f(size_f);

    SECTION("Standard eigenproblem") {
      auto Aop = [&](vcv_t in, vv_t out) { out = prod(A, in); };

      solver_t ar(A.size1());
      testing.standard_custom_exact_shifts(ar, Aop, shifts_f);
    }

    SECTION("Generalized eigenproblem: Shift-and-Invert mode") {
      decltype(A) op_mat = prod(inverse(A - sigma * M), M);

      auto op = [&](vv_t in, vv_t out) { out = prod(op_mat, in); };
      auto Bop = [&](vcv_t in, vv_t out) { out = prod(M, in); };

      solver_t ar(A.size1());
      testing.generalized_custom_exact_shifts(ar, op, Bop, shifts_f, sigma);
    }
  }
}
