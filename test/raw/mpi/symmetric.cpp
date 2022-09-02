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

///////////////////
// Test invert() //
///////////////////
TEST_CASE("Symmetric matrix is inverted", "[invert_symmetric]") {
  const int N = 100;
  const double diag_coeff_mean = 1.0;
  const int offdiag_offset = 3;
  const double offdiag_coeff_mean = -0.1;
  const double offdiag_coeff_diff = 0;

  // Symmetric matrix A
  auto A = make_sparse_matrix<ezarpack::Symmetric>(
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

////////////////////////////////////////////////
// Eigenproblems with real symmetric matrices //
////////////////////////////////////////////////

TEST_CASE("Symmetric eigenproblem is solved", "[solver_symmetric]") {

  using solver_t = mpi::arpack_solver<ezarpack::Symmetric, raw_storage>;

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
    auto Aop = [&](vcv_t in, vv_t out) { mat_vec(A.get(), in, out); };

    solver_t ar(N, MPI_COMM_WORLD);
    testing.standard_eigenproblems(ar, Aop);
  }

  SECTION("Generalized eigenproblem: invert mode") {
    auto invM = make_buffer<double>(N * N);
    invert(M.get(), invM.get(), N);

    solver_t ar(N, MPI_COMM_WORLD);

    auto tmp = make_buffer<double>(ar.local_block_size());
    auto op = [&](vv_t in, vv_t out) {
      mat_vec(A.get(), in, tmp.get());
      std::copy(tmp.get(), tmp.get() + ar.local_block_size(), in);
      mat_vec(invM.get(), in, out);
    };
    auto Bop = [&](vcv_t in, vv_t out) { mat_vec(M.get(), in, out); };

    testing.generalized_eigenproblems(ar, solver_t::Inverse, op, Bop);
  }

  SECTION("Generalized eigenproblem: Shift-and-Invert mode") {
    auto AmM = make_buffer<double>(N * N);
    for(int i = 0; i < N; ++i) {
      for(int j = 0; j < N; ++j) {
        AmM[i + j * N] = A[i + j * N] - sigma * M[i + j * N];
      }
    }
    auto invAmM = make_buffer<double>(N * N);
    invert(AmM.get(), invAmM.get(), N);
    auto op_mat = make_buffer<double>(N * N);
    mm_prod(invAmM.get(), M.get(), op_mat.get(), N);

    auto op = [&](vv_t in, vv_t out) { mat_vec(op_mat.get(), in, out); };
    auto Bop = [&](vcv_t in, vv_t out) { mat_vec(M.get(), in, out); };

    solver_t ar(N, MPI_COMM_WORLD);
    testing.generalized_eigenproblems(ar, solver_t::ShiftAndInvert, op, Bop,
                                      sigma);
  }

  SECTION("Generalized eigenproblem: Buckling mode") {

    auto MmA = make_buffer<double>(N * N);
    for(int i = 0; i < N; ++i) {
      for(int j = 0; j < N; ++j) {
        MmA[i + j * N] = M[i + j * N] - sigma * A[i + j * N];
      }
    }
    auto invMmA = make_buffer<double>(N * N);
    invert(MmA.get(), invMmA.get(), N);
    auto op_mat = make_buffer<double>(N * N);
    mm_prod(invMmA.get(), M.get(), op_mat.get(), N);

    auto op = [&](vv_t in, vv_t out) { mat_vec(op_mat.get(), in, out); };
    auto Bop = [&](vcv_t in, vv_t out) { mat_vec(M.get(), in, out); };

    solver_t ar(N, MPI_COMM_WORLD);
    const int ncv = 30;
    testing.generalized_eigenproblems(ar, solver_t::Buckling, op, Bop, sigma,
                                      ncv);
  }

  SECTION("Generalized eigenproblem: Cayley transformed mode") {
    auto AmM = make_buffer<double>(N * N);
    auto ApM = make_buffer<double>(N * N);
    for(int i = 0; i < N; ++i) {
      for(int j = 0; j < N; ++j) {
        AmM[i + j * N] = A[i + j * N] - sigma * M[i + j * N];
        ApM[i + j * N] = A[i + j * N] + sigma * M[i + j * N];
      }
    }
    auto invAmM = make_buffer<double>(N * N);
    invert(AmM.get(), invAmM.get(), N);
    auto op_mat = make_buffer<double>(N * N);
    mm_prod(invAmM.get(), ApM.get(), op_mat.get(), N);

    auto op = [&](vv_t in, vv_t out) { mat_vec(op_mat.get(), in, out); };
    auto Bop = [&](vcv_t in, vv_t out) { mat_vec(M.get(), in, out); };

    solver_t ar(N, MPI_COMM_WORLD);
    testing.generalized_eigenproblems(ar, solver_t::Cayley, op, Bop, sigma);
  }

  SECTION("Indirect access to workspace vectors") {
    solver_t ar(N, MPI_COMM_WORLD);

    auto Aop = [&](vcv_t, vv_t) {
      auto in = ar.workspace_vector(ar.in_vector_n());
      auto out = ar.workspace_vector(ar.out_vector_n());
      mat_vec(A.get(), in, out);
    };

    testing.standard_eigenproblems(ar, Aop);

    CHECK_THROWS(ar.workspace_vector(-1));
    CHECK_THROWS(ar.workspace_vector(3));
  }

  SECTION("Skip computation of eigenvectors") {
    solver_t ar(N, MPI_COMM_WORLD);

    SECTION("Standard eigenproblem") {
      auto Aop = [&](vcv_t in, vv_t out) { mat_vec(A.get(), in, out); };

      testing.standard_skip_eigenvectors(ar, Aop);
    }

    SECTION("Generalized eigenproblem: invert mode") {
      auto invM = make_buffer<double>(N * N);
      invert(M.get(), invM.get(), N);
      auto tmp = make_buffer<double>(ar.local_block_size());
      auto op = [&](vv_t in, vv_t out) {
        mat_vec(A.get(), in, tmp.get());
        std::copy(tmp.get(), tmp.get() + ar.local_block_size(), in);
        mat_vec(invM.get(), in, out);
      };
      auto Bop = [&](vcv_t in, vv_t out) { mat_vec(M.get(), in, out); };

      testing.generalized_skip_eigenvectors(ar, op, Bop);
    }
  }
}
