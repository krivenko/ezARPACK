/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2026 Igor Krivenko
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

TEST_CASE("Asymmetric eigenproblem is solved", "[solver_asymmetric]") {

  using solver_t =
      ezarpack::mpi::arpack_solver<ezarpack::Asymmetric, triqs_storage>;

  const int N = 100;
  const double diag_coeff_mean = 1.0;
  const int offdiag_offset = 3;
  const double offdiag_coeff_mean = -1.0;
  const double offdiag_coeff_diff = 0.1;
  const int nev = 8;

  const dcomplex sigma(0.5, 0.5);

  // Asymmetric matrix A
  auto A = make_sparse_matrix<ezarpack::Asymmetric>(
      N, diag_coeff_mean, offdiag_offset, offdiag_coeff_mean,
      offdiag_coeff_diff);
  // Inner product matrix
  auto M = make_inner_prod_matrix<ezarpack::Asymmetric>(N);

  // Testing helper
  auto testing = make_testing_helper<solver_t>(A, M, N, nev);

  // Matrix-distributed vector multiplication
  auto mat_vec = mpi_mat_vec<false>(N, MPI_COMM_WORLD);

  using vv_t = solver_t::vector_view_t;
  using vcv_t = solver_t::vector_const_view_t;

  SECTION("Constructors") { test_mpi_arpack_solver_ctor<solver_t>(); }

  SECTION("Standard eigenproblem") {
    auto Aop = [&](vcv_t in, vv_t out) { mat_vec(A, in, out); };

    solver_t ar(first_dim(A), MPI_COMM_WORLD);
    const int ncv = 30;
    testing.standard_eigenproblems(ar, Aop, ncv);
  }

  SECTION("Generalized eigenproblem: invert mode") {
    decltype(A) op_mat = inverse(M) * A;

    auto op = [&](vcv_t in, vv_t out) { mat_vec(op_mat, in, out); };
    auto Bop = [&](vcv_t in, vv_t out) { mat_vec(M, in, out); };

    solver_t ar(first_dim(A), MPI_COMM_WORLD);
    testing.generalized_eigenproblems(ar, solver_t::Inverse, op, Bop);
  }

  SECTION("Generalized eigenproblem: Shift-and-Invert mode (real part)") {
    decltype(A) op_mat = real(inverse(A - sigma * M) * M);

    auto op = [&](vcv_t in, vv_t out) { mat_vec(op_mat, in, out); };
    auto Bop = [&](vcv_t in, vv_t out) { mat_vec(M, in, out); };

    solver_t ar(first_dim(A), MPI_COMM_WORLD);
    testing.generalized_eigenproblems(ar, solver_t::ShiftAndInvertReal, op, Bop,
                                      sigma);
  }

  SECTION("Generalized eigenproblem: Shift-and-Invert mode (imaginary part)") {
    decltype(A) op_mat = imag(inverse(A - sigma * M) * M);

    auto op = [&](vcv_t in, vv_t out) { mat_vec(op_mat, in, out); };
    auto Bop = [&](vcv_t in, vv_t out) { mat_vec(M, in, out); };

    solver_t ar(first_dim(A), MPI_COMM_WORLD);
    testing.generalized_eigenproblems(ar, solver_t::ShiftAndInvertImag, op, Bop,
                                      sigma);
  }

  SECTION("Indirect access to workspace vectors") {
    solver_t ar(first_dim(A), MPI_COMM_WORLD);

    auto Aop = [&](vcv_t, vv_t) {
      auto in = ar.workspace_vector(ar.in_vector_n());
      auto out = ar.workspace_vector(ar.out_vector_n());
      mat_vec(A, in, out);
    };

    const int ncv = 30;
    testing.standard_eigenproblems(ar, Aop, ncv);

    CHECK_THROWS(ar.workspace_vector(-1));
    CHECK_THROWS(ar.workspace_vector(3));
  }

  SECTION("Various compute_vectors") {
    solver_t ar(first_dim(A), MPI_COMM_WORLD);

    SECTION("Standard eigenproblem") {
      auto Aop = [&](vcv_t in, vv_t out) { mat_vec(A, in, out); };

      testing.standard_compute_vectors(ar, Aop);
    }

    SECTION("Generalized eigenproblem: invert mode") {
      decltype(A) op_mat = inverse(M) * A;

      auto op = [&](vcv_t in, vv_t out) { mat_vec(op_mat, in, out); };
      auto Bop = [&](vcv_t in, vv_t out) { mat_vec(M, in, out); };

      testing.generalized_compute_vectors(ar, op, Bop);
    }
  }

  SECTION("Custom implementation of the Exact Shift Strategy") {
    using rvcv_t = solver_t::real_vector_const_view_t;
    using rvv_t = solver_t::real_vector_view_t;
    auto size_f = [](rvv_t shifts) -> int { return shifts.size(); };
    exact_shift_strategy<ezarpack::Asymmetric, rvcv_t, rvv_t> shifts_f(size_f);

    SECTION("Standard eigenproblem") {
      auto Aop = [&](vcv_t in, vv_t out) { mat_vec(A, in, out); };

      solver_t ar(first_dim(A), MPI_COMM_WORLD);
      testing.standard_custom_exact_shifts(ar, Aop, shifts_f);
    }

    SECTION("Generalized eigenproblem: Shift-and-Invert mode (real part)") {
      decltype(A) op_mat = real(inverse(A - sigma * M) * M);

      auto op = [&](vcv_t in, vv_t out) { mat_vec(op_mat, in, out); };
      auto Bop = [&](vcv_t in, vv_t out) { mat_vec(M, in, out); };

      solver_t ar(first_dim(A), MPI_COMM_WORLD);
      testing.generalized_custom_exact_shifts(ar, op, Bop, shifts_f, sigma);
    }
  }
}
