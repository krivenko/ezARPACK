/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2024 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/
#pragma once

#include <stdexcept>

#include "ezarpack/common.hpp"

// Helper class used to test solvers
template<typename SolverType, typename MatrixType> class testing_helper;

//
// Specialization of testing_helper for symmetric eigenproblems
//

template<template<ezarpack::operator_kind, typename> class SolverTemplate,
         typename Backend,
         typename MatrixType>
class testing_helper<SolverTemplate<ezarpack::Symmetric, Backend>, MatrixType> {
  // Matrix of the eigenproblem
  MatrixType const& A;
  // Inner product matrix
  MatrixType const& M;
  // Dimension of the eigenproblem
  int N;
  // Number of requested eigenvalues
  int nev;

  using solver_t = SolverTemplate<ezarpack::Symmetric, Backend>;
  using params_t = typename solver_t::params_t;

  const typename params_t::eigenvalues_select_t spectrum_parts[5] = {
      params_t::Smallest, params_t::Largest, params_t::SmallestMagnitude,
      params_t::LargestMagnitude, params_t::BothEnds};

public:
  testing_helper(MatrixType const& A, MatrixType const& M, int N, int nev)
      : A(A), M(M), N(N), nev(nev) {}

  // Test standard eigenproblems
  template<typename AOpType>
  void standard_eigenproblems(solver_t& ar, AOpType const& Aop, int ncv = 0) {
    REQUIRE(ar.dim() == N);

    for(auto e : spectrum_parts) {
      params_t params(nev, e, true);
      params.random_residual_vector = false;
      if(ncv > 0) params.ncv = ncv;
      set_init_residual_vector(ar);
      ar(Aop, params);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors(ar, A);
      check_basis_vectors(ar);
    }
  }

  // Test generalized eigenproblems
  template<typename OpType, typename BOpType>
  void generalized_eigenproblems(solver_t& ar,
                                 typename solver_t::Mode mode,
                                 OpType const& op,
                                 BOpType const& Bop,
                                 double sigma = 0,
                                 int ncv = 0) {
    REQUIRE(ar.dim() == N);

    for(auto e : spectrum_parts) {
      params_t params(nev, e, true);
      params.random_residual_vector = false;
      if(mode == solver_t::ShiftAndInvert || mode == solver_t::Buckling ||
         mode == solver_t::Cayley)
        params.sigma = sigma;

      if(ncv > 0) params.ncv = ncv;
      set_init_residual_vector(ar);
      ar(op, Bop, mode, params);
      CHECK(ar.nconv() >= nev);
      if(mode == solver_t::Buckling)
        check_eigenvectors(ar, M, A);
      else
        check_eigenvectors(ar, A, M);
      check_basis_vectors(ar, M);
    }
  }

  // Skip computation of eigenvectors: Standard eigenproblem
  template<typename OpType>
  void standard_skip_eigenvectors(solver_t& ar, OpType const& op) {
    REQUIRE(ar.dim() == N);

    params_t params(nev, params_t::LargestMagnitude, false);
    params.random_residual_vector = false;

    set_init_residual_vector(ar);
    ar(op, params);
    CHECK(ar.nconv() >= nev);
    CHECK_THROWS_AS(ar.eigenvectors(), std::runtime_error);
  }

  // Skip computation of eigenvectors: Generalized eigenproblem
  template<typename OpType, typename BOpType>
  void generalized_skip_eigenvectors(solver_t& ar,
                                     OpType const& op,
                                     BOpType const& Bop) {
    REQUIRE(ar.dim() == N);

    params_t params(nev, params_t::LargestMagnitude, false);
    params.random_residual_vector = false;

    set_init_residual_vector(ar);
    ar(op, Bop, solver_t::Inverse, params);
    CHECK(ar.nconv() >= nev);
    CHECK_THROWS_AS(ar.eigenvectors(), std::runtime_error);
  }

  // Custom implementation of the Exact Shift Strategy (standard eigenproblem)
  template<typename AOpType, typename ShiftsF>
  void standard_custom_exact_shifts(solver_t& ar,
                                    AOpType const& Aop,
                                    ShiftsF&& shifts_f) {
    REQUIRE(ar.dim() == N);

    params_t params(nev, params_t::LargestMagnitude, true);
    params.random_residual_vector = false;
    set_init_residual_vector(ar);
    ar(Aop, params, shifts_f);
    CHECK(ar.nconv() >= nev);
    check_eigenvectors(ar, A);
    check_basis_vectors(ar);
  }

  // Custom implementation of the Exact Shift Strategy (generalized
  // eigenproblem)
  template<typename OpType, typename BOpType, typename ShiftsF>
  void generalized_custom_exact_shifts(solver_t& ar,
                                       OpType const& op,
                                       BOpType const& Bop,
                                       ShiftsF&& shifts_f,
                                       double sigma) {
    REQUIRE(ar.dim() == N);

    params_t params(nev, params_t::LargestMagnitude, true);
    params.sigma = sigma;
    params.random_residual_vector = false;
    set_init_residual_vector(ar);
    ar(op, Bop, solver_t::ShiftAndInvert, params, shifts_f);
    CHECK(ar.nconv() >= nev);
    check_eigenvectors(ar, A, M);
    check_basis_vectors(ar, M);
  }
};

//
// Specialization of testing_helper for asymmetric eigenproblems
//

// Various compute_vectors (standard eigenproblems)
template<typename SolverType, typename OpType, typename MatrixType>
void standard_compute_vectors_impl(SolverType& ar,
                                   MatrixType const& A,
                                   OpType const& op,
                                   int N,
                                   int nev) {
  using params_t = typename SolverType::params_t;
  REQUIRE(ar.dim() == N);

  params_t params(nev, params_t::LargestMagnitude, params_t::None);
  params.random_residual_vector = false;

  set_init_residual_vector(ar);
  ar(op, params);
  CHECK_THROWS_AS(ar.schur_vectors(), std::runtime_error);
  CHECK_THROWS_AS(ar.eigenvectors(), std::runtime_error);

  params.compute_vectors = params_t::Schur;
  set_init_residual_vector(ar);
  ar(op, params);
  check_basis_vectors(ar);
  CHECK_THROWS_AS(ar.eigenvectors(), std::runtime_error);

  params.compute_vectors = params_t::Ritz;
  set_init_residual_vector(ar);
  ar(op, params);
  check_basis_vectors(ar);
  check_eigenvectors(ar, A);
}

template<typename SolverType,
         typename OpType,
         typename BOpType,
         typename MatrixType>
void generalized_compute_vectors_impl(SolverType& ar,
                                      MatrixType const& A,
                                      MatrixType const& M,
                                      OpType const& op,
                                      BOpType const& Bop,
                                      int N,
                                      int nev) {
  using params_t = typename SolverType::params_t;
  REQUIRE(ar.dim() == N);

  params_t params(nev, params_t::LargestMagnitude, params_t::None);
  params.random_residual_vector = false;

  set_init_residual_vector(ar);
  ar(op, Bop, SolverType::Inverse, params);
  CHECK_THROWS_AS(ar.schur_vectors(), std::runtime_error);
  CHECK_THROWS_AS(ar.eigenvectors(), std::runtime_error);

  params.compute_vectors = params_t::Schur;
  set_init_residual_vector(ar);
  ar(op, Bop, SolverType::Inverse, params);
  check_basis_vectors(ar, M);
  CHECK_THROWS_AS(ar.eigenvectors(), std::runtime_error);

  params.compute_vectors = params_t::Ritz;
  set_init_residual_vector(ar);
  ar(op, Bop, SolverType::Inverse, params);
  check_basis_vectors(ar, M);
  check_eigenvectors(ar, A, M);
}

template<template<ezarpack::operator_kind, typename> class SolverTemplate,
         typename Backend,
         typename MatrixType>
class testing_helper<SolverTemplate<ezarpack::Asymmetric, Backend>,
                     MatrixType> {
  // Matrix of the eigenproblem
  MatrixType const& A;
  // Inner product matrix
  MatrixType const& M;
  // Dimension of the eigenproblem
  int N;
  // Number of requested eigenvalues
  int nev;

  using solver_t = SolverTemplate<ezarpack::Asymmetric, Backend>;
  using params_t = typename solver_t::params_t;

  const typename params_t::eigenvalues_select_t spectrum_parts[6] = {
      params_t::LargestMagnitude, params_t::SmallestMagnitude,
      params_t::LargestReal,      params_t::SmallestReal,
      params_t::LargestImag,      params_t::SmallestImag};

public:
  testing_helper(MatrixType const& A, MatrixType const& M, int N, int nev)
      : A(A), M(M), N(N), nev(nev) {}

  // Test standard eigenproblems
  template<typename AOpType>
  void standard_eigenproblems(solver_t& ar, AOpType const& Aop, int ncv = 0) {
    REQUIRE(ar.dim() == N);

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.random_residual_vector = false;
      if(ncv > 0) params.ncv = ncv;
      set_init_residual_vector(ar);
      ar(Aop, params);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors(ar, A);
      check_basis_vectors(ar);
    }
  }

  // Test generalized eigenproblems
  template<typename OpType, typename BOpType>
  void generalized_eigenproblems(solver_t& ar,
                                 typename solver_t::Mode mode,
                                 OpType const& op,
                                 BOpType const& Bop,
                                 dcomplex sigma = 0,
                                 int ncv = 0) {
    REQUIRE(ar.dim() == N);

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.random_residual_vector = false;
      if(mode == solver_t::ShiftAndInvertReal ||
         mode == solver_t::ShiftAndInvertImag)
        params.sigma = sigma;

      if(ncv > 0) params.ncv = ncv;
      set_init_residual_vector(ar);
      ar(op, Bop, mode, params);
      CHECK(ar.nconv() >= nev);
      if(mode == solver_t::ShiftAndInvertReal ||
         mode == solver_t::ShiftAndInvertImag)
        check_eigenvectors_shift_and_invert(ar, A, M);
      else
        check_eigenvectors(ar, A, M);
      check_basis_vectors(ar, M);
    }
  }

  // Various compute_vectors (standard eigenproblems)
  template<typename OpType>
  void standard_compute_vectors(solver_t& ar, OpType const& op) {
    standard_compute_vectors_impl(ar, A, op, N, nev);
  }

  // Various compute_vectors (generalized eigenproblems)
  template<typename OpType, typename BOpType>
  void generalized_compute_vectors(solver_t& ar,
                                   OpType const& op,
                                   BOpType const& Bop) {
    generalized_compute_vectors_impl(ar, A, M, op, Bop, N, nev);
  }

  // Custom implementation of the Exact Shift Strategy (standard eigenproblem)
  template<typename AOpType, typename ShiftsF>
  void standard_custom_exact_shifts(solver_t& ar,
                                    AOpType const& Aop,
                                    ShiftsF&& shifts_f) {
    REQUIRE(ar.dim() == N);

    params_t params(nev, params_t::LargestMagnitude, params_t::Ritz);
    params.random_residual_vector = false;
    set_init_residual_vector(ar);
    ar(Aop, params, shifts_f);
    CHECK(ar.nconv() >= nev);
    check_eigenvectors(ar, A);
    check_basis_vectors(ar);
  }

  // Custom implementation of the Exact Shift Strategy (generalized
  // eigenproblem)
  template<typename OpType, typename BOpType, typename ShiftsF>
  void generalized_custom_exact_shifts(solver_t& ar,
                                       OpType const& op,
                                       BOpType const& Bop,
                                       ShiftsF&& shifts_f,
                                       dcomplex sigma) {
    REQUIRE(ar.dim() == N);

    params_t params(nev, params_t::LargestMagnitude, params_t::Ritz);
    params.random_residual_vector = false;
    params.sigma = sigma;
    set_init_residual_vector(ar);
    ar(op, Bop, solver_t::ShiftAndInvertReal, params, shifts_f);
    CHECK(ar.nconv() >= nev);
    check_eigenvectors_shift_and_invert(ar, A, M);
    check_basis_vectors(ar, M);
  }
};

//
// Specialization of testing_helper for complex eigenproblems
//

template<template<ezarpack::operator_kind, typename> class SolverTemplate,
         typename Backend,
         typename MatrixType>
class testing_helper<SolverTemplate<ezarpack::Complex, Backend>, MatrixType> {
  // Matrix of the eigenproblem
  MatrixType const& A;
  // Inner product matrix
  MatrixType const& M;
  // Dimension of the eigenproblem
  int N;
  // Number of requested eigenvalues
  int nev;

  using solver_t = SolverTemplate<ezarpack::Complex, Backend>;
  using params_t = typename solver_t::params_t;

  const typename params_t::eigenvalues_select_t spectrum_parts[6] = {
      params_t::LargestMagnitude, params_t::SmallestMagnitude,
      params_t::LargestReal,      params_t::SmallestReal,
      params_t::LargestImag,      params_t::SmallestImag};

public:
  testing_helper(MatrixType const& A, MatrixType const& M, int N, int nev)
      : A(A), M(M), N(N), nev(nev) {}

  // Test standard eigenproblems
  template<typename AOpType>
  void standard_eigenproblems(solver_t& ar, AOpType const& Aop, int ncv = 0) {
    REQUIRE(ar.dim() == N);

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.random_residual_vector = false;
      if(ncv > 0) params.ncv = ncv;
      set_init_residual_vector(ar);
      ar(Aop, params);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors(ar, A);
      check_basis_vectors(ar);
    }
  }

  // Test generalized eigenproblems
  template<typename OpType, typename BOpType>
  void generalized_eigenproblems(solver_t& ar,
                                 typename solver_t::Mode mode,
                                 OpType const& op,
                                 BOpType const& Bop,
                                 dcomplex sigma = 0,
                                 int ncv = 0) {
    REQUIRE(ar.dim() == N);

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.random_residual_vector = false;
      if(mode == solver_t::ShiftAndInvert) params.sigma = sigma;

      if(ncv > 0) params.ncv = ncv;
      set_init_residual_vector(ar);
      ar(op, Bop, mode, params);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors(ar, A, M);
      check_basis_vectors(ar, M);
    }
  }

  // Various compute_vectors (standard eigenproblems)
  template<typename OpType>
  void standard_compute_vectors(solver_t& ar, OpType const& op) {
    standard_compute_vectors_impl(ar, A, op, N, nev);
  }

  // Various compute_vectors (generalized eigenproblems)
  template<typename OpType, typename BOpType>
  void generalized_compute_vectors(solver_t& ar,
                                   OpType const& op,
                                   BOpType const& Bop) {
    generalized_compute_vectors_impl(ar, A, M, op, Bop, N, nev);
  }

  // Custom implementation of the Exact Shift Strategy (standard eigenproblem)
  template<typename AOpType, typename ShiftsF>
  void standard_custom_exact_shifts(solver_t& ar,
                                    AOpType const& Aop,
                                    ShiftsF&& shifts_f) {
    REQUIRE(ar.dim() == N);

    params_t params(nev, params_t::LargestMagnitude, params_t::Ritz);
    params.random_residual_vector = false;
    set_init_residual_vector(ar);
    ar(Aop, params, shifts_f);
    CHECK(ar.nconv() >= nev);
    check_eigenvectors(ar, A);
    check_basis_vectors(ar);
  }

  // Custom implementation of the Exact Shift Strategy (generalized
  // eigenproblem)
  template<typename OpType, typename BOpType, typename ShiftsF>
  void generalized_custom_exact_shifts(solver_t& ar,
                                       OpType const& op,
                                       BOpType const& Bop,
                                       ShiftsF&& shifts_f,
                                       dcomplex sigma) {
    REQUIRE(ar.dim() == N);

    params_t params(nev, params_t::LargestMagnitude, params_t::Ritz);
    params.sigma = sigma;
    params.random_residual_vector = false;
    set_init_residual_vector(ar);
    ar(op, Bop, solver_t::ShiftAndInvert, params, shifts_f);
    CHECK(ar.nconv() >= nev);
    check_eigenvectors(ar, A, M);
    check_basis_vectors(ar, M);
  }
};

template<typename SolverType, typename MatrixType>
testing_helper<SolverType, MatrixType>
make_testing_helper(MatrixType const& A, MatrixType const& M, int N, int nev) {
  return testing_helper<SolverType, MatrixType>(A, M, N, nev);
}
