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

///////////////////
// Test invert() //
///////////////////
TEST_CASE("Asymmetric matrix is inverted", "[invert_asymmetric]") {
  const int N = 100;
  const double diag_coeff_mean = 1.0;
  const int offdiag_offset = 3;
  const double offdiag_coeff_mean = -1.0;
  const double offdiag_coeff_diff = 0.1;

  // Asymmetric matrix A
  auto A = make_sparse_matrix<ezarpack::Asymmetric>(
      N, diag_coeff_mean, offdiag_offset, offdiag_coeff_mean,
      offdiag_coeff_diff);

  auto invA = make_buffer<double>(N * N);
  invert(A.get(), invA.get(), N);

  auto id = make_buffer<double>(N * N);
  for(int i = 0; i < N; ++i) {
    for(int j = 0; j < N; ++j) {
      id[i + j * N] = i == j;
    }
  }

  auto prod = make_buffer<double>(N * N);

  // A * A^{-1}
  mm_prod(A.get(), invA.get(), prod.get(), N);
  CHECK_THAT(prod.get(), IsCloseTo(id.get(), N * N));

  // A^{-1} * A
  mm_prod(invA.get(), A.get(), prod.get(), N);
  CHECK_THAT(prod.get(), IsCloseTo(id.get(), N * N));
}

//////////////////////////////////////////////
// Eigenproblems with general real matrices //
//////////////////////////////////////////////

TEST_CASE("Asymmetric eigenproblem is solved", "[solver_asymmetric]") {

  using solver_t = mpi::arpack_solver<ezarpack::Asymmetric, raw_storage>;

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
    auto Aop = [&](vcv_t in, vv_t out) { mat_vec(A.get(), in, out); };

    solver_t ar(N, MPI_COMM_WORLD);
    const int ncv = 30;
    testing.standard_eigenproblems(ar, Aop, ncv);
  }

  SECTION("Generalized eigenproblem: invert mode") {
    auto invM = make_buffer<double>(N * N);
    invert(M.get(), invM.get(), N);
    auto op_mat = make_buffer<double>(N * N);
    mm_prod(invM.get(), A.get(), op_mat.get(), N);

    auto op = [&](vcv_t in, vv_t out) { mat_vec(op_mat.get(), in, out); };
    auto Bop = [&](vcv_t in, vv_t out) { mat_vec(M.get(), in, out); };

    solver_t ar(N, MPI_COMM_WORLD);
    testing.generalized_eigenproblems(ar, solver_t::Inverse, op, Bop);
  }

  SECTION("Generalized eigenproblem: Shift-and-Invert mode (real part)") {

    auto AmM = make_buffer<dcomplex>(N * N);
    for(int i = 0; i < N; ++i) {
      for(int j = 0; j < N; ++j) {
        AmM[i + j * N] = A[i + j * N] - sigma * M[i + j * N];
      }
    }
    auto invAmM = make_buffer<dcomplex>(N * N);
    invert(AmM.get(), invAmM.get(), N);
    auto op_mat = make_buffer<dcomplex>(N * N);
    mm_prod(invAmM.get(), M.get(), op_mat.get(), N);
    auto op_mat_re = make_buffer<double>(N * N);
    for(int i = 0; i < N; ++i) {
      for(int j = 0; j < N; ++j) {
        op_mat_re[i + j * N] = op_mat[i + j * N].real();
      }
    }

    auto op = [&](vcv_t in, vv_t out) { mat_vec(op_mat_re.get(), in, out); };
    auto Bop = [&](vcv_t in, vv_t out) { mat_vec(M.get(), in, out); };

    solver_t ar(N, MPI_COMM_WORLD);
    testing.generalized_eigenproblems(ar, solver_t::ShiftAndInvertReal, op, Bop,
                                      sigma);
  }

  SECTION("Generalized eigenproblem: Shift-and-Invert mode (imaginary part)") {

    auto AmM = make_buffer<dcomplex>(N * N);
    for(int i = 0; i < N; ++i) {
      for(int j = 0; j < N; ++j) {
        AmM[i + j * N] = A[i + j * N] - sigma * M[i + j * N];
      }
    }
    auto invAmM = make_buffer<dcomplex>(N * N);
    invert(AmM.get(), invAmM.get(), N);
    auto op_mat = make_buffer<dcomplex>(N * N);
    mm_prod(invAmM.get(), M.get(), op_mat.get(), N);
    auto op_mat_im = make_buffer<double>(N * N);
    for(int i = 0; i < N; ++i) {
      for(int j = 0; j < N; ++j) {
        op_mat_im[i + j * N] = op_mat[i + j * N].imag();
      }
    }

    auto op = [&](vcv_t in, vv_t out) { mat_vec(op_mat_im.get(), in, out); };
    auto Bop = [&](vcv_t in, vv_t out) { mat_vec(M.get(), in, out); };

    solver_t ar(N, MPI_COMM_WORLD);
    testing.generalized_eigenproblems(ar, solver_t::ShiftAndInvertImag, op, Bop,
                                      sigma);
  }

  SECTION("Indirect access to workspace vectors") {
    solver_t ar(N, MPI_COMM_WORLD);

    auto Aop = [&](vcv_t, vv_t) {
      auto in = ar.workspace_vector(ar.in_vector_n());
      auto out = ar.workspace_vector(ar.out_vector_n());
      mat_vec(A.get(), in, out);
    };

    const int ncv = 30;
    testing.standard_eigenproblems(ar, Aop, ncv);

    CHECK_THROWS(ar.workspace_vector(-1));
    CHECK_THROWS(ar.workspace_vector(3));
  }

  SECTION("Various compute_vectors") {
    solver_t ar(N, MPI_COMM_WORLD);

    SECTION("Standard eigenproblem") {
      auto Aop = [&](vcv_t in, vv_t out) { mat_vec(A.get(), in, out); };

      testing.standard_compute_vectors(ar, Aop);
    }

    SECTION("Generalized eigenproblem: invert mode") {
      auto invM = make_buffer<double>(N * N);
      invert(M.get(), invM.get(), N);
      auto op_mat = make_buffer<double>(N * N);
      mm_prod(invM.get(), A.get(), op_mat.get(), N);

      auto op = [&](vcv_t in, vv_t out) { mat_vec(op_mat.get(), in, out); };
      auto Bop = [&](vcv_t in, vv_t out) { mat_vec(M.get(), in, out); };

      testing.generalized_compute_vectors(ar, op, Bop);
    }
  }
}
