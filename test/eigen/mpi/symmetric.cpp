/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2024 Igor Krivenko
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

  using solver_t = mpi::arpack_solver<ezarpack::Symmetric, eigen_storage>;

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

  // Matrix-distributed vector multiplication
  auto mat_vec = mpi_mat_vec<false>(N, MPI_COMM_WORLD);

  using vv_t = solver_t::vector_view_t;
  using vcv_t = solver_t::vector_const_view_t;

  SECTION("Constructors") { test_mpi_arpack_solver_ctor<solver_t>(); }

  SECTION("Standard eigenproblem") {
    auto Aop = [&](vcv_t in, vv_t out) { mat_vec(A, in, out); };

    solver_t ar(A.rows(), MPI_COMM_WORLD);
    testing.standard_eigenproblems(ar, Aop);
  }

  SECTION("Generalized eigenproblem: invert mode") {
    decltype(A) invM = M.inverse();

    auto op = [&](vv_t in, vv_t out) {
      mat_vec(A, in, in);
      mat_vec(invM, in, out);
    };
    auto Bop = [&](vcv_t in, vv_t out) { mat_vec(M, in, out); };

    solver_t ar(A.rows(), MPI_COMM_WORLD);
    testing.generalized_eigenproblems(ar, solver_t::Inverse, op, Bop);
  }

  SECTION("Generalized eigenproblem: Shift-and-Invert mode") {
    decltype(A) op_mat = (A - sigma * M).inverse() * M;

    auto op = [&](vv_t in, vv_t out) { mat_vec(op_mat, in, out); };
    auto Bop = [&](vcv_t in, vv_t out) { mat_vec(M, in, out); };

    solver_t ar(A.rows(), MPI_COMM_WORLD);
    testing.generalized_eigenproblems(ar, solver_t::ShiftAndInvert, op, Bop,
                                      sigma);
  }

  SECTION("Generalized eigenproblem: Buckling mode") {
    decltype(A) op_mat = (M - sigma * A).inverse() * M;

    auto op = [&](vv_t in, vv_t out) { mat_vec(op_mat, in, out); };
    auto Bop = [&](vcv_t in, vv_t out) { mat_vec(M, in, out); };

    solver_t ar(A.rows(), MPI_COMM_WORLD);
    const int ncv = 30;
    testing.generalized_eigenproblems(ar, solver_t::Buckling, op, Bop, sigma,
                                      ncv);
  }

  SECTION("Generalized eigenproblem: Cayley transformed mode") {
    decltype(A) op_mat = (A - sigma * M).inverse() * (A + sigma * M);

    auto op = [&](vv_t in, vv_t out) { mat_vec(op_mat, in, out); };
    auto Bop = [&](vcv_t in, vv_t out) { mat_vec(M, in, out); };

    solver_t ar(A.rows(), MPI_COMM_WORLD);
    testing.generalized_eigenproblems(ar, solver_t::Cayley, op, Bop, sigma);
  }

  SECTION("Indirect access to workspace vectors") {
    solver_t ar(A.rows(), MPI_COMM_WORLD);

    auto Aop = [&](vcv_t, vv_t) {
      auto in = ar.workspace_vector(ar.in_vector_n());
      auto out = ar.workspace_vector(ar.out_vector_n());
      mat_vec(A, in, out);
    };

    testing.standard_eigenproblems(ar, Aop);

    CHECK_THROWS(ar.workspace_vector(-1));
    CHECK_THROWS(ar.workspace_vector(3));
  }

  SECTION("Skip computation of eigenvectors") {
    solver_t ar(A.rows(), MPI_COMM_WORLD);

    SECTION("Standard eigenproblem") {
      auto Aop = [&](vcv_t in, vv_t out) { mat_vec(A, in, out); };

      testing.standard_skip_eigenvectors(ar, Aop);
    }

    SECTION("Generalized eigenproblem: invert mode") {
      decltype(A) invM = M.inverse();

      auto op = [&](vv_t in, vv_t out) {
        mat_vec(A, in, in);
        mat_vec(invM, in, out);
      };
      auto Bop = [&](vcv_t in, vv_t out) { mat_vec(M, in, out); };

      testing.generalized_skip_eigenvectors(ar, op, Bop);
    }
  }

  SECTION("Custom implementation of the Exact Shift Strategy") {
    using rvcv_t = solver_t::real_vector_const_view_t;
    using rvv_t = solver_t::real_vector_view_t;
    auto size_f = [](rvv_t shifts) -> int { return shifts.size(); };
    exact_shift_strategy<ezarpack::Symmetric, rvcv_t, rvv_t> shifts_f(size_f);

    SECTION("Standard eigenproblem") {
      auto Aop = [&](vcv_t in, vv_t out) { mat_vec(A, in, out); };

      solver_t ar(A.rows(), MPI_COMM_WORLD);
      testing.standard_custom_exact_shifts(ar, Aop, shifts_f);
    }

    SECTION("Generalized eigenproblem: Shift-and-Invert mode") {
      decltype(A) op_mat = (A - sigma * M).inverse() * M;

      auto op = [&](vv_t in, vv_t out) { mat_vec(op_mat, in, out); };
      auto Bop = [&](vcv_t in, vv_t out) { mat_vec(M, in, out); };

      solver_t ar(A.rows(), MPI_COMM_WORLD);
      testing.generalized_custom_exact_shifts(ar, op, Bop, shifts_f, sigma);
    }
  }
}
