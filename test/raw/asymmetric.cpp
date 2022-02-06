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
TEST_CASE("Asymmetric matrix is inverted", "[invert_asymmetric]") {
  const int N = 100;
  const double diag_coeff_shift = -0.55;
  const double diag_coeff_amp = 1.0;
  const int offdiag_offset = 3;
  const double offdiag_coeff = -1.05;

  auto A = make_sparse_matrix<Asymmetric>(N, diag_coeff_shift, diag_coeff_amp,
                                          offdiag_offset, offdiag_coeff);

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
  mm_product(A.get(), invA.get(), prod.get(), N);
  CHECK_THAT(prod.get(), IsCloseTo(id.get(), N * N));

  // A^{-1} * A
  mm_product(invA.get(), A.get(), prod.get(), N);
  CHECK_THAT(prod.get(), IsCloseTo(id.get(), N * N));
}

//////////////////////////////////////////////
// Eigenproblems with general real matrices //
//////////////////////////////////////////////

TEST_CASE("Asymmetric eigenproblem is solved", "[solver_asymmetric]") {

  using solver_t = arpack_solver<Asymmetric, raw_storage>;
  using params_t = solver_t::params_t;

  const int N = 100;
  const double diag_coeff_shift = -0.55;
  const double diag_coeff_amp = 1.0;
  const int offdiag_offset = 3;
  const double offdiag_coeff = -1.05;
  const int nev = 10;

  auto spectrum_parts = {
      params_t::LargestMagnitude, params_t::SmallestMagnitude,
      params_t::LargestReal,      params_t::SmallestReal,
      params_t::LargestImag,      params_t::SmallestImag};

  // Asymmetric matrix A
  auto A = make_sparse_matrix<Asymmetric>(N, diag_coeff_shift, diag_coeff_amp,
                                          offdiag_offset, offdiag_coeff);
  // Inner product matrix
  auto M = make_inner_prod_matrix<Asymmetric>(N);

  auto set_init_residual_vector = [](solver_t& ar) {
    for(int i = 0; i < N; ++i)
      ar.residual_vector()[i] = double(i) / N;
  };

  using vector_view_t = solver_t::vector_view_t;
  using vector_const_view_t = solver_t::vector_const_view_t;

  SECTION("Standard eigenproblem") {
    auto Aop = [&](vector_const_view_t in, vector_view_t out) {
      mv_product(A.get(), in, out, N);
    };

    solver_t ar(N);

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(Aop, params);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors(ar, A.get(), N, nev);
      check_basis_vectors(ar, N, nev);
    }
  }

  SECTION("Generalized eigenproblem: invert mode") {
    auto invM = make_buffer<double>(N * N);
    invert(M.get(), invM.get(), N);
    auto op_matrix = make_buffer<double>(N * N);
    mm_product(invM.get(), A.get(), op_matrix.get(), N);

    auto op = [&](vector_const_view_t in, vector_view_t out) {
      mv_product(op_matrix.get(), in, out, N);
    };
    auto Bop = [&](vector_const_view_t in, vector_view_t out) {
      mv_product(M.get(), in, out, N);
    };

    solver_t ar(N);

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(op, Bop, solver_t::Inverse, params);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors(ar, A.get(), M.get(), N, nev);
      check_basis_vectors(ar, M.get(), N, nev);
    }
  }

  SECTION("Generalized eigenproblem: Shift-and-Invert mode (real part)") {
    dcomplex sigma(1.0, -0.1);

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
    auto op_matrix_re = make_buffer<double>(N * N);
    for(int i = 0; i < N; ++i) {
      for(int j = 0; j < N; ++j) {
        op_matrix_re[i + j * N] = op_matrix[i + j * N].real();
      }
    }

    auto op = [&](vector_const_view_t in, vector_view_t out) {
      mv_product(op_matrix_re.get(), in, out, N);
    };
    auto Bop = [&](vector_const_view_t in, vector_view_t out) {
      mv_product(M.get(), in, out, N);
    };

    solver_t ar(N);

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.random_residual_vector = false;
      params.sigma = sigma;
      set_init_residual_vector(ar);
      ar(op, Bop, solver_t::ShiftAndInvertReal, params);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors_shift_and_invert(ar, A.get(), M.get(), N, nev);
      check_basis_vectors(ar, M.get(), N, nev);
    }
  }

  SECTION("Generalized eigenproblem: Shift-and-Invert mode (imaginary part)") {
    dcomplex sigma(-0.1, 1.0);

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
    auto op_matrix_im = make_buffer<double>(N * N);
    for(int i = 0; i < N; ++i) {
      for(int j = 0; j < N; ++j) {
        op_matrix_im[i + j * N] = op_matrix[i + j * N].imag();
      }
    }

    auto op = [&](vector_const_view_t in, vector_view_t out) {
      mv_product(op_matrix_im.get(), in, out, N);
    };
    auto Bop = [&](vector_const_view_t in, vector_view_t out) {
      mv_product(M.get(), in, out, N);
    };

    solver_t ar(N);

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.random_residual_vector = false;
      params.sigma = sigma;
      set_init_residual_vector(ar);
      ar(op, Bop, solver_t::ShiftAndInvertImag, params);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors_shift_and_invert(ar, A.get(), M.get(), N, nev);
      check_basis_vectors(ar, M.get(), N, nev);
    }
  }

  SECTION("Indirect access to workspace vectors") {
    solver_t ar(N);

    auto Aop = [&](vector_const_view_t, vector_view_t) {
      auto in = ar.workspace_vector(ar.in_vector_n());
      auto out = ar.workspace_vector(ar.out_vector_n());
      mv_product(A.get(), in, out, N);
    };

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(Aop, params);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors(ar, A.get(), N, nev);
      check_basis_vectors(ar, N, nev);
    }

    CHECK_THROWS(ar.workspace_vector(-1));
    CHECK_THROWS(ar.workspace_vector(3));
  }
}
