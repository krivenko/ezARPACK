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

  using solver_t = arpack_solver<ezarpack::Symmetric, eigen_storage>;

  const int N = 100;
  const double diag_coeff_shift = -0.55;
  const double diag_coeff_amp = 1.0;
  const int offdiag_offset = 3;
  const double offdiag_coeff = -1.05;
  const int nev = 10;

  // Symmetric matrix A
  auto A = make_sparse_matrix<ezarpack::Symmetric>(
      N, diag_coeff_shift, diag_coeff_amp, offdiag_offset, offdiag_coeff);
  // Inner product matrix
  auto M = make_inner_prod_matrix<ezarpack::Symmetric>(N);

  // Testing helper
  auto testing = make_testing_helper<solver_t>(A, M, N, nev);

  using vv_t = solver_t::vector_view_t;
  using vcv_t = solver_t::vector_const_view_t;

  SECTION("Standard eigenproblem") {
    auto Aop = [&](vcv_t in, vv_t out) { out = A * in; };

    solver_t ar(A.rows());
    testing.standard_eigenproblems(ar, Aop);
  }

  SECTION("Generalized eigenproblem: invert mode") {
    decltype(A) invM = M.inverse();

    auto op = [&](vv_t in, vv_t out) {
      in = A * in;
      out = invM * in;
    };
    auto Bop = [&](vcv_t in, vv_t out) { out = M * in; };

    solver_t ar(A.rows());
    testing.generalized_eigenproblems(ar, solver_t::Inverse, op, Bop);
  }

  SECTION("Generalized eigenproblem: Shift-and-Invert mode") {
    double sigma = 2.0;
    decltype(A) op_mat = (A - sigma * M).inverse() * M;

    auto op = [&](vv_t in, vv_t out) { out = op_mat * in; };
    auto Bop = [&](vcv_t in, vv_t out) { out = M * in; };

    solver_t ar(A.rows());
    testing.generalized_eigenproblems(ar, solver_t::ShiftAndInvert, op, Bop,
                                      sigma);
  }

  SECTION("Generalized eigenproblem: Buckling mode") {
    double sigma = 3.3;
    decltype(A) op_mat = (M - sigma * A).inverse() * M;

    auto op = [&](vv_t in, vv_t out) { out = op_mat * in; };
    auto Bop = [&](vcv_t in, vv_t out) { out = M * in; };

    solver_t ar(A.rows());
    const int ncv = 30;
    testing.generalized_eigenproblems(ar, solver_t::Buckling, op, Bop, sigma,
                                      ncv);
  }

  SECTION("Generalized eigenproblem: Cayley transformed mode") {
    double sigma = 2.0;
    decltype(A) op_mat = (A - sigma * M).inverse() * (A + sigma * M);

    auto op = [&](vv_t in, vv_t out) { out = op_mat * in; };
    auto Bop = [&](vcv_t in, vv_t out) { out = M * in; };

    solver_t ar(A.rows());
    testing.generalized_eigenproblems(ar, solver_t::Cayley, op, Bop, sigma);
  }

  SECTION("Indirect access to workspace vectors") {
    solver_t ar(A.rows());

    auto Aop = [&](vcv_t, vv_t) {
      auto in = ar.workspace_vector(ar.in_vector_n());
      auto out = ar.workspace_vector(ar.out_vector_n());
      out = A * in;
    };

    testing.standard_eigenproblems(ar, Aop);

    CHECK_THROWS(ar.workspace_vector(-1));
    CHECK_THROWS(ar.workspace_vector(3));
  }

  SECTION("Skip computation of eigenvectors") {
    solver_t ar(A.rows());

    SECTION("Standard eigenproblem") {
      auto Aop = [&](vcv_t in, vv_t out) { out = A * in; };

      testing.standard_skip_eigenvectors(ar, Aop);
    }

    SECTION("Generalized eigenproblem: invert mode") {
      decltype(A) invM = M.inverse();

      auto op = [&](vv_t in, vv_t out) {
        in = A * in;
        out = invM * in;
      };
      auto Bop = [&](vcv_t in, vv_t out) { out = M * in; };

      testing.generalized_skip_eigenvectors(ar, op, Bop);
    }
  }

  SECTION("Custom implementation of the Exact Shift Strategy") {
    using rvcv_t = solver_t::real_vector_const_view_t;
    using rvv_t = solver_t::real_vector_view_t;
    auto size_f = [](rvv_t shifts) -> int { return shifts.size(); };
    exact_shift_strategy<ezarpack::Symmetric, rvcv_t, rvv_t> shifts_f(size_f);

    SECTION("Standard eigenproblem") {
      auto Aop = [&](vcv_t in, vv_t out) { out = A * in; };

      solver_t ar(A.rows());
      testing.standard_custom_exact_shifts(ar, Aop, shifts_f);
    }

    SECTION("Generalized eigenproblem: Shift-and-Invert mode") {
      double sigma = 2.0;
      decltype(A) op_mat = (A - sigma * M).inverse() * M;

      auto op = [&](vv_t in, vv_t out) { out = op_mat * in; };
      auto Bop = [&](vcv_t in, vv_t out) { out = M * in; };

      solver_t ar(A.rows());
      testing.generalized_custom_exact_shifts(ar, op, Bop, shifts_f, sigma);
    }
  }
}
