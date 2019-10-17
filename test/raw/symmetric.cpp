/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2019 Igor Krivenko <igor.s.krivenko@gmail.com>
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
  const double diag_coeff = 1.5;
  const int offdiag_offset = 3;
  const double offdiag_coeff = 0.5;

  auto A = make_sparse_matrix<Symmetric>(N,
                                         diag_coeff,
                                         offdiag_offset,
                                         offdiag_coeff);

  auto invA = make_buffer<double>(N * N);
  invert(A.get(), invA.get(), N);

  auto id = make_buffer<double>(N * N);
  for(int i = 0; i < N; ++i) {
    for(int j = 0; j < N; ++j) {
      id[i + j*N] = i == j;
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

////////////////////////////////////////////////
// Eigenproblems with real symmetric matrices //
////////////////////////////////////////////////

TEST_CASE("Symmetric eigenproblem is solved", "[worker_symmetric]") {

  using worker_t = arpack_worker<Symmetric, raw_storage>;
  using params_t = worker_t::params_t;

  const int N = 100;
  const double diag_coeff = 1.5;
  const int offdiag_offset = 3;
  const double offdiag_coeff = 0.5;
  const int nev = 10;

  auto spectrum_parts = {params_t::Smallest,
                         params_t::Largest,
                         params_t::SmallestMagnitude,
                         params_t::LargestMagnitude,
                         params_t::BothEnds};

  // Symmetric matrix A
  auto A = make_sparse_matrix<Symmetric>(N,
                                         diag_coeff,
                                         offdiag_offset,
                                         offdiag_coeff);
  // Inner product matrix
  auto M = make_inner_prod_matrix<Symmetric>(N);

  auto set_init_residual_vector = [N](worker_t & ar) {
    for(int i = 0; i < N; ++i) ar.residual_vector()[i] = double(i) / N;
  };

  using vector_view_t = worker_t::vector_view_t;
  using vector_const_view_t = worker_t::vector_const_view_t;

  SECTION("Standard eigenproblem") {
    auto Aop = [&](vector_const_view_t from, vector_view_t to) {
      mv_product(A.get(), from, to, N);
    };

    worker_t ar(N);

    for(auto e : spectrum_parts) {
      params_t params(nev, e, true);
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(Aop, params);
      check_eigenvectors(ar, A.get(), N, nev);
    }
  }

  SECTION("Generalized eigenproblem: invert mode") {
    auto invM = make_buffer<double>(N * N);
    invert(M.get(), invM.get(), N);

    auto tmp = make_buffer<double>(N);
    auto op = [&](vector_view_t from, vector_view_t to) {
      mv_product(A.get(), from, tmp.get(), N);
      std::copy(tmp.get(), tmp.get() + N, from);

      mv_product(invM.get(), from, to, N);
    };
    auto Bop = [&](vector_const_view_t from, vector_view_t to) {
      mv_product(M.get(), from, to, N);
    };

    worker_t ar(N);

    for(auto e : spectrum_parts) {
      params_t params(nev, e, true);
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(op, Bop, worker_t::Invert, params);
      check_eigenvectors(ar, A.get(), M.get(), N, nev);
    }
  }

  SECTION("Generalized eigenproblem: Shift-and-Invert mode") {
    double sigma = 1.0;

    auto AmM = make_buffer<double>(N * N);
    for(int i = 0; i < N; ++i) {
      for(int j = 0; j < N; ++j) {
        AmM[i + j*N] = A[i + j*N] - sigma * M[i + j*N];
      }
    }
    auto invAmM = make_buffer<double>(N * N);
    invert(AmM.get(), invAmM.get(), N);
    auto op_matrix = make_buffer<double>(N * N);
    mm_product(invAmM.get(), M.get(), op_matrix.get(), N);

    auto op = [&](vector_view_t from, vector_view_t to) {
      mv_product(op_matrix.get(), from, to, N);
    };
    auto Bop = [&](vector_const_view_t from, vector_view_t to) {
      mv_product(M.get(), from, to, N);
    };

    worker_t ar(N);

    for(auto e : spectrum_parts) {
      params_t params(nev, e, true);
      params.sigma = sigma;
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(op, Bop, worker_t::ShiftAndInvert, params);
      check_eigenvectors(ar, A.get(), M.get(), N, nev);
    }
  }

  SECTION("Generalized eigenproblem: Buckling mode") {
    double sigma = 1.0;

    auto MmA = make_buffer<double>(N * N);
    for(int i = 0; i < N; ++i) {
      for(int j = 0; j < N; ++j) {
        MmA[i + j*N] = M[i + j*N] - sigma * A[i + j*N];
      }
    }
    auto invMmA = make_buffer<double>(N * N);
    invert(MmA.get(), invMmA.get(), N);
    auto op_matrix = make_buffer<double>(N * N);
    mm_product(invMmA.get(), M.get(), op_matrix.get(), N);

    auto op = [&](vector_view_t from, vector_view_t to) {
      mv_product(op_matrix.get(), from, to, N);
    };
    auto Bop = [&](vector_const_view_t from, vector_view_t to) {
      mv_product(M.get(), from, to, N);
    };

    worker_t ar(N);

    for(auto e : spectrum_parts) {
      params_t params(nev, e, true);
      params.sigma = sigma;
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(op, Bop, worker_t::Buckling, params);
      check_eigenvectors(ar, M.get(), A.get(), N, nev);
    }
  }

  SECTION("Generalized eigenproblem: Cayley transformed mode") {
    double sigma = 1.0;

    auto AmM = make_buffer<double>(N * N);
    auto ApM = make_buffer<double>(N * N);
    for(int i = 0; i < N; ++i) {
      for(int j = 0; j < N; ++j) {
        AmM[i + j*N] = A[i + j*N] - sigma * M[i + j*N];
        ApM[i + j*N] = A[i + j*N] + sigma * M[i + j*N];
      }
    }
    auto invAmM = make_buffer<double>(N * N);
    invert(AmM.get(), invAmM.get(), N);
    auto op_matrix = make_buffer<double>(N * N);
    mm_product(invAmM.get(), ApM.get(), op_matrix.get(), N);

    auto op = [&](vector_view_t from, vector_view_t to) {
      mv_product(op_matrix.get(), from, to, N);
    };
    auto Bop = [&](vector_const_view_t from, vector_view_t to) {
      mv_product(M.get(), from, to, N);
    };

    worker_t ar(N);

    for(auto e : spectrum_parts) {
      params_t params(nev, e, true);
      params.sigma = sigma;
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(op, Bop, worker_t::Cayley, params);
      check_eigenvectors(ar, A.get(), M.get(), N, nev);
    }
  }
}
