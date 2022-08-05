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

  using solver_t = mpi::arpack_solver<ezarpack::Complex, eigen_storage>;
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
  auto A = make_sparse_matrix<ezarpack::Complex>(
      N, diag_coeff_shift, diag_coeff_amp, offdiag_offset, offdiag_coeff);
  // Inner product matrix
  auto M = make_inner_prod_matrix<ezarpack::Complex>(N);

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

  SECTION("Constructors") {
    std::vector<std::vector<unsigned int>> block_sizes = {
        {100}, {50, 50}, {34, 33, 33}, {25, 25, 25, 25}};
    std::vector<std::vector<unsigned int>> block_starts = {
        {0}, {0, 50}, {0, 34, 67}, {0, 25, 50, 75}};

    auto comm_size = mpi::size(MPI_COMM_WORLD);
    auto comm_rank = mpi::rank(MPI_COMM_WORLD);

    if(comm_size <= 4) {
      solver_t ar1(N, MPI_COMM_WORLD);
      CHECK(ar1.local_block_size() == block_sizes[comm_size - 1][comm_rank]);
      CHECK(ar1.local_block_start() == block_starts[comm_size - 1][comm_rank]);

      solver_t ar2(block_sizes[comm_size - 1], MPI_COMM_WORLD);
      CHECK(ar2.local_block_size() == block_sizes[comm_size - 1][comm_rank]);
      CHECK(ar2.local_block_start() == block_starts[comm_size - 1][comm_rank]);
    }
  }

  SECTION("Standard eigenproblem") {
    auto Aop = [&](vector_const_view_t in, vector_view_t out) {
      mat_vec(A, in, out);
    };

    solver_t ar(A.rows(), MPI_COMM_WORLD);
    REQUIRE(ar.dim() == A.rows());

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(Aop, params);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors(ar, A);
      check_basis_vectors(ar);
    }
  }

  SECTION("Generalized eigenproblem: invert mode") {
    decltype(A) op_matrix = M.inverse() * A;

    auto op = [&](vector_const_view_t in, vector_view_t out) {
      mat_vec(op_matrix, in, out);
    };
    auto Bop = [&](vector_const_view_t in, vector_view_t out) {
      mat_vec(M, in, out);
    };

    solver_t ar(A.rows(), MPI_COMM_WORLD);
    REQUIRE(ar.dim() == A.rows());

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(op, Bop, solver_t::Inverse, params);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors(ar, A, M);
      check_basis_vectors(ar, M);
    }
  }

  SECTION("Generalized eigenproblem: Shift-and-Invert mode") {
    dcomplex sigma(-1.0, 0.9);
    decltype(A) op_matrix = (A - sigma * M).inverse() * M;

    auto op = [&](vector_const_view_t in, vector_view_t out) {
      mat_vec(op_matrix, in, out);
    };
    auto Bop = [&](vector_const_view_t in, vector_view_t out) {
      mat_vec(M, in, out);
    };

    solver_t ar(A.rows(), MPI_COMM_WORLD);
    REQUIRE(ar.dim() == A.rows());

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.sigma = sigma;
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(op, Bop, solver_t::ShiftAndInvert, params);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors(ar, A, M);
      check_basis_vectors(ar, M);
    }
  }

  SECTION("Indirect access to workspace vectors") {
    solver_t ar(A.rows(), MPI_COMM_WORLD);
    REQUIRE(ar.dim() == A.rows());

    auto Aop = [&](vector_const_view_t, vector_view_t) {
      auto in = ar.workspace_vector(ar.in_vector_n());
      auto out = ar.workspace_vector(ar.out_vector_n());
      mat_vec(A, in, out);
    };

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(Aop, params);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors(ar, A);
      check_basis_vectors(ar);
    }

    CHECK_THROWS(ar.workspace_vector(-1));
    CHECK_THROWS(ar.workspace_vector(3));
  }

  SECTION("Custom implementation of the Exact Shift Strategy") {
    std::vector<int> p;
    p.reserve(A.rows());
    auto shifts_f = [&](solver_t::complex_vector_const_view_t ritz_values,
                        solver_t::complex_vector_const_view_t ritz_bounds,
                        solver_t::complex_vector_view_t shifts) {
      int np = shifts.size();
      if(np == 0) return;

      p.resize(np);
      std::iota(p.begin(), p.end(), 0);
      // p will contain the permutation that puts ritz_bounds in
      // the descending order of magnitude
      std::sort(p.begin(), p.end(), [&](int i, int j) {
        return std::abs(ritz_bounds(i)) > std::abs(ritz_bounds(j));
      });
      // Apply permutation p to ritz_values and use the result to fill shifts
      for(int i = 0; i < np; ++i)
        shifts(i) = ritz_values(p[i]);
    };

    SECTION("Standard eigenproblem") {
      auto Aop = [&](vector_const_view_t in, vector_view_t out) {
        mat_vec(A, in, out);
      };

      solver_t ar(A.rows(), MPI_COMM_WORLD);
      REQUIRE(ar.dim() == A.rows());

      params_t params(nev, params_t::LargestMagnitude, params_t::Ritz);
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(Aop, params, shifts_f);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors(ar, A);
      check_basis_vectors(ar);
    }

    SECTION("Generalized eigenproblem: Shift-and-Invert mode") {
      dcomplex sigma(-1.0, 0.9);
      decltype(A) op_matrix = (A - sigma * M).inverse() * M;

      auto op = [&](vector_const_view_t in, vector_view_t out) {
        mat_vec(op_matrix, in, out);
      };
      auto Bop = [&](vector_const_view_t in, vector_view_t out) {
        mat_vec(M, in, out);
      };

      solver_t ar(A.rows(), MPI_COMM_WORLD);
      REQUIRE(ar.dim() == A.rows());

      params_t params(nev, params_t::LargestMagnitude, params_t::Ritz);
      params.sigma = sigma;
      params.random_residual_vector = false;
      set_init_residual_vector(ar);
      ar(op, Bop, solver_t::ShiftAndInvert, params, shifts_f);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors(ar, A, M);
      check_basis_vectors(ar, M);
    }
  }
}
