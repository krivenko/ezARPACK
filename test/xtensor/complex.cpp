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

/////////////////////////////////////////
// Eigenproblems with complex matrices //
/////////////////////////////////////////

TEST_CASE("Complex eigenproblem is solved", "[solver_complex]") {

  using xt::linalg::dot;
  using xt::linalg::inv;

  using solver_t = arpack_solver<ezarpack::Complex, xtensor_storage>;

  const int N = 100;
  const dcomplex diag_coeff_shift = -0.55;
  const dcomplex diag_coeff_amp = 1.0;
  const int offdiag_offset = 3;
  const dcomplex offdiag_coeff(0.5, 0.5);
  const int nev = 10;

  // Hermitian matrix A
  auto A = make_sparse_matrix<ezarpack::Complex>(
      N, diag_coeff_shift, diag_coeff_amp, offdiag_offset, offdiag_coeff);
  // Inner product matrix
  auto M = make_inner_prod_matrix<ezarpack::Complex>(N);

  // Testing helper
  auto testing = make_testing_helper<solver_t>(A, M, N, nev);

  using vv_t = solver_t::vector_view_t;
  using vcv_t = solver_t::vector_const_view_t;

  SECTION("Standard eigenproblem") {
    auto Aop = [&](vcv_t in, vv_t out) { out = dot(A, in); };

    solver_t ar(A.shape(0));
    testing.standard_eigenproblems(ar, Aop);
  }

  SECTION("Generalized eigenproblem: invert mode") {
    decltype(A) op_mat = dot(inv(M), A);

    auto op = [&](vcv_t in, vv_t out) { out = dot(op_mat, in); };
    auto Bop = [&](vcv_t in, vv_t out) { out = dot(M, in); };

    solver_t ar(A.shape(0));
    testing.generalized_eigenproblems(ar, solver_t::Inverse, op, Bop);
  }

  SECTION("Generalized eigenproblem: Shift-and-Invert mode") {
    dcomplex sigma(-1.0, 0.9);
    decltype(A) op_mat = dot(inv(eval(A - sigma * M)), M);

    auto op = [&](vcv_t in, vv_t out) { out = dot(op_mat, in); };
    auto Bop = [&](vcv_t in, vv_t out) { out = dot(M, in); };

    solver_t ar(A.shape(0));
    testing.generalized_eigenproblems(ar, solver_t::ShiftAndInvert, op, Bop,
                                      sigma);
  }

  SECTION("Indirect access to workspace vectors") {
    solver_t ar(A.shape(0));

    auto Aop = [&](vcv_t, vv_t) {
      auto in = ar.workspace_vector(ar.in_vector_n());
      auto out = ar.workspace_vector(ar.out_vector_n());
      out = dot(A, in);
    };

    testing.standard_eigenproblems(ar, Aop);

    CHECK_THROWS(ar.workspace_vector(-1));
    CHECK_THROWS(ar.workspace_vector(3));
  }

  SECTION("Various compute_vectors") {
    solver_t ar(A.shape(0));

    SECTION("Standard eigenproblem") {
      auto Aop = [&](vcv_t in, vv_t out) { out = dot(A, in); };

      testing.standard_compute_vectors(ar, Aop);
    }

    SECTION("Generalized eigenproblem: invert mode") {
      decltype(A) op_mat = dot(inv(M), A);

      auto op = [&](vcv_t in, vv_t out) { out = dot(op_mat, in); };
      auto Bop = [&](vcv_t in, vv_t out) { out = dot(M, in); };

      testing.generalized_compute_vectors(ar, op, Bop);
    }
  }

  SECTION("Custom implementation of the Exact Shift Strategy") {
    using cvcv_t = solver_t::complex_vector_const_view_t;
    using cvv_t = solver_t::complex_vector_view_t;
    auto size_f = [](cvv_t shifts) -> int { return shifts.size(); };
    exact_shift_strategy<ezarpack::Complex, cvcv_t, cvv_t> shifts_f(size_f);

    SECTION("Standard eigenproblem") {
      auto Aop = [&](vcv_t in, vv_t out) { out = dot(A, in); };

      solver_t ar(A.shape(0));
      testing.standard_custom_exact_shifts(ar, Aop, shifts_f);
    }

    SECTION("Generalized eigenproblem: Shift-and-Invert mode") {
      dcomplex sigma(-1.0, 0.9);
      decltype(A) op_mat = dot(inv(eval(A - sigma * M)), M);

      auto op = [&](vcv_t in, vv_t out) { out = dot(op_mat, in); };
      auto Bop = [&](vcv_t in, vv_t out) { out = dot(M, in); };

      solver_t ar(A.shape(0));
      testing.generalized_custom_exact_shifts(ar, op, Bop, shifts_f, sigma);
    }
  }
}
