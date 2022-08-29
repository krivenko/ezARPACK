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
TEST_CASE("Complex matrix is inverted", "[invert_asymmetric]") {
  const int N = 100;
  const dcomplex diag_coeff_shift = -0.55;
  const dcomplex diag_coeff_amp = 1.0;
  const int offdiag_offset = 3;
  const dcomplex offdiag_coeff(0.5, 0.5);

  auto A = make_sparse_matrix<Complex>(N, diag_coeff_shift, diag_coeff_amp,
                                       offdiag_offset, offdiag_coeff);

  auto invA = make_buffer<dcomplex>(N * N);
  invert(A.get(), invA.get(), N);

  auto id = make_buffer<dcomplex>(N * N);
  for(int i = 0; i < N; ++i) {
    for(int j = 0; j < N; ++j) {
      id[i + j * N] = i == j;
    }
  }

  auto prod = make_buffer<dcomplex>(N * N);

  // A * A^{-1}
  mm_product(A.get(), invA.get(), prod.get(), N);
  CHECK_THAT(prod.get(), IsCloseTo(id.get(), N * N));

  // A^{-1} * A
  mm_product(invA.get(), A.get(), prod.get(), N);
  CHECK_THAT(prod.get(), IsCloseTo(id.get(), N * N));
}

/////////////////////////////////////////
// Eigenproblems with complex matrices //
/////////////////////////////////////////

TEST_CASE("Complex eigenproblem is solved", "[solver_complex]") {

  using solver_t = mpi::arpack_solver<Complex, raw_storage>;
  using params_t = solver_t::params_t;

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

  auto set_init_residual_vector = [](solver_t& ar) {
    int const block_start = ar.local_block_start();
    int const block_size = ar.local_block_size();
    for(int i = 0; i < block_size; ++i)
      ar.residual_vector()[i] = double(i + block_start) / N;
  };

  // Matrix-distributed vector multiplication
  auto mat_vec = mpi_mat_vec<true>(N, MPI_COMM_WORLD);

  using vector_view_t = solver_t::vector_view_t;
  using vector_const_view_t = solver_t::vector_const_view_t;

  SECTION("Standard eigenproblem") {
    auto Aop = [&](vector_const_view_t in, vector_view_t out) {
      mat_vec(A.get(), in, out);
    };

    solver_t ar(N, MPI_COMM_WORLD);
    REQUIRE(ar.dim() == N);

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(Aop, params);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors(ar, A.get(), nev);
      check_basis_vectors(ar, nev);
    }
  }

  SECTION("Generalized eigenproblem: invert mode") {
    auto invM = make_buffer<dcomplex>(N * N);
    invert(M.get(), invM.get(), N);
    auto op_matrix = make_buffer<dcomplex>(N * N);
    mm_product(invM.get(), A.get(), op_matrix.get(), N);

    auto op = [&](vector_const_view_t in, vector_view_t out) {
      mat_vec(op_matrix.get(), in, out);
    };
    auto Bop = [&](vector_const_view_t in, vector_view_t out) {
      mat_vec(M.get(), in, out);
    };

    solver_t ar(N, MPI_COMM_WORLD);
    REQUIRE(ar.dim() == N);

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(op, Bop, solver_t::Inverse, params);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors(ar, A.get(), M.get(), nev);
      check_basis_vectors(ar, M.get(), nev);
    }
  }

  SECTION("Generalized eigenproblem: Shift-and-Invert mode") {
    dcomplex sigma(-1.0, 0.9);

    auto AmM = make_buffer<dcomplex>(N * N);
    for(int i = 0; i < N; ++i) {
      for(int j = 0; j < N; ++j) {
        AmM[i + j * N] = A[i + j * N] - sigma * M[i + j * N];
      }
    }
    auto invAmM = make_buffer<dcomplex>(N * N);
    invert(AmM.get(), invAmM.get(), N);
    auto op_matrix = make_buffer<dcomplex>(N * N);
    mm_product(invAmM.get(), M.get(), op_matrix.get(), N);

    auto op = [&](vector_const_view_t in, vector_view_t out) {
      mat_vec(op_matrix.get(), in, out);
    };
    auto Bop = [&](vector_const_view_t in, vector_view_t out) {
      mat_vec(M.get(), in, out);
    };

    solver_t ar(N, MPI_COMM_WORLD);
    REQUIRE(ar.dim() == N);

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.sigma = sigma;
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(op, Bop, solver_t::ShiftAndInvert, params);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors(ar, A.get(), M.get(), nev);
      check_basis_vectors(ar, M.get(), nev);
    }
  }

  SECTION("Indirect access to workspace vectors") {
    solver_t ar(N, MPI_COMM_WORLD);
    REQUIRE(ar.dim() == N);

    auto Aop = [&](vector_const_view_t, vector_view_t) {
      auto in = ar.workspace_vector(ar.in_vector_n());
      auto out = ar.workspace_vector(ar.out_vector_n());
      mat_vec(A.get(), in, out);
    };

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(Aop, params);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors(ar, A.get(), nev);
      check_basis_vectors(ar, nev);
    }

    CHECK_THROWS(ar.workspace_vector(-1));
    CHECK_THROWS(ar.workspace_vector(3));
  }

  SECTION("Various compute_vectors") {
    solver_t ar(N, MPI_COMM_WORLD);
    REQUIRE(ar.dim() == N);

    SECTION("Standard eigenproblem") {
      auto Aop = [&](vector_const_view_t in, vector_view_t out) {
        mat_vec(A.get(), in, out);
      };

      params_t params(nev, params_t::LargestMagnitude, params_t::None);
      params.random_residual_vector = false;

      set_init_residual_vector(ar);
      ar(Aop, params);
      CHECK_THROWS_AS(ar.schur_vectors(), std::runtime_error);
      CHECK_THROWS_AS(ar.eigenvectors(), std::runtime_error);

      params.compute_vectors = params_t::Schur;
      set_init_residual_vector(ar);
      ar(Aop, params);
      check_basis_vectors(ar, nev);
      CHECK_THROWS_AS(ar.eigenvectors(), std::runtime_error);

      params.compute_vectors = params_t::Ritz;
      set_init_residual_vector(ar);
      ar(Aop, params);
      check_basis_vectors(ar, nev);
      check_eigenvectors(ar, A.get(), nev);
    }

    SECTION("Generalized eigenproblem: invert mode") {
      auto invM = make_buffer<dcomplex>(N * N);
      invert(M.get(), invM.get(), N);
      auto op_matrix = make_buffer<dcomplex>(N * N);
      mm_product(invM.get(), A.get(), op_matrix.get(), N);

      auto op = [&](vector_const_view_t in, vector_view_t out) {
        mat_vec(op_matrix.get(), in, out);
      };
      auto Bop = [&](vector_const_view_t in, vector_view_t out) {
        mat_vec(M.get(), in, out);
      };

      params_t params(nev, params_t::LargestMagnitude, params_t::None);
      params.random_residual_vector = false;

      set_init_residual_vector(ar);
      ar(op, Bop, solver_t::Inverse, params);
      CHECK_THROWS_AS(ar.schur_vectors(), std::runtime_error);
      CHECK_THROWS_AS(ar.eigenvectors(), std::runtime_error);

      params.compute_vectors = params_t::Schur;
      set_init_residual_vector(ar);
      ar(op, Bop, solver_t::Inverse, params);
      check_basis_vectors(ar, M.get(), nev);
      CHECK_THROWS_AS(ar.eigenvectors(), std::runtime_error);

      params.compute_vectors = params_t::Ritz;
      set_init_residual_vector(ar);
      ar(op, Bop, solver_t::Inverse, params);
      check_basis_vectors(ar, M.get(), nev);
      check_eigenvectors(ar, A.get(), M.get(), nev);
    }
  }
}
