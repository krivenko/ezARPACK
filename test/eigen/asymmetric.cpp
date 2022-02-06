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

//////////////////////////////////////////////
// Eigenproblems with general real matrices //
//////////////////////////////////////////////

TEST_CASE("Asymmetric eigenproblem is solved", "[solver_asymmetric]") {

  using solver_t = arpack_solver<ezarpack::Asymmetric, eigen_storage>;
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
  auto A = make_sparse_matrix<ezarpack::Asymmetric>(
      N, diag_coeff_shift, diag_coeff_amp, offdiag_offset, offdiag_coeff);
  // Inner product matrix
  auto M = make_inner_prod_matrix<ezarpack::Asymmetric>(N);

  auto set_init_residual_vector = [](solver_t& ar) {
    for(int i = 0; i < N; ++i)
      ar.residual_vector()[i] = double(i) / N;
  };

  using vector_view_t = solver_t::vector_view_t;
  using vector_const_view_t = solver_t::vector_const_view_t;

  SECTION("Standard eigenproblem") {
    auto Aop = [&](vector_const_view_t in, vector_view_t out) {
      out = A * in;
    };

    solver_t ar(A.rows());

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
      out = op_matrix * in;
    };
    auto Bop = [&](vector_const_view_t in, vector_view_t out) {
      out = M * in;
    };

    solver_t ar(A.rows());

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

  SECTION("Generalized eigenproblem: Shift-and-Invert mode (real part)") {
    dcomplex sigma(1.0, -0.1);
#ifdef EIGEN_CAN_MIX_REAL_COMPLEX_EXPR
    decltype(A) op_matrix = ((A - sigma * M).inverse() * M).real();
#else
    decltype(A) op_matrix =
        ((A.cast<dcomplex>() - sigma * M.cast<dcomplex>()).inverse() *
         M.cast<dcomplex>())
            .real();
#endif

    auto op = [&](vector_const_view_t in, vector_view_t out) {
      out = op_matrix * in;
    };
    auto Bop = [&](vector_const_view_t in, vector_view_t out) {
      out = M * in;
    };

    solver_t ar(A.rows());

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.random_residual_vector = false;
      params.sigma = sigma;
      set_init_residual_vector(ar);
      ar(op, Bop, solver_t::ShiftAndInvertReal, params);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors_shift_and_invert(ar, A, M);
      check_basis_vectors(ar, M);
    }
  }

  SECTION("Generalized eigenproblem: Shift-and-Invert mode (imaginary part)") {
    dcomplex sigma(-0.1, 1.0);
#ifdef EIGEN_CAN_MIX_REAL_COMPLEX_EXPR
    decltype(A) op_matrix = ((A - sigma * M).inverse() * M).imag();
#else
    decltype(A) op_matrix =
        ((A.cast<dcomplex>() - sigma * M.cast<dcomplex>()).inverse() *
         M.cast<dcomplex>())
            .imag();
#endif

    auto op = [&](vector_const_view_t in, vector_view_t out) {
      out = op_matrix * in;
    };
    auto Bop = [&](vector_const_view_t in, vector_view_t out) {
      out = M * in;
    };

    solver_t ar(A.rows());

    for(auto e : spectrum_parts) {
      params_t params(nev, e, params_t::Ritz);
      params.random_residual_vector = false;
      params.sigma = sigma;
      set_init_residual_vector(ar);
      ar(op, Bop, solver_t::ShiftAndInvertImag, params);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors_shift_and_invert(ar, A, M);
      check_basis_vectors(ar, M);
    }
  }

  SECTION("Indirect access to workspace vectors") {
    solver_t ar(A.rows());

    auto Aop = [&](vector_const_view_t, vector_view_t) {
      auto in = ar.workspace_vector(ar.in_vector_n());
      auto out = ar.workspace_vector(ar.out_vector_n());
      out = A * in;
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
    auto shifts_f = [&](solver_t::real_vector_const_view_t ritz_values_re,
                        solver_t::real_vector_const_view_t ritz_values_im,
                        solver_t::real_vector_const_view_t ritz_bounds,
                        solver_t::real_vector_view_t shifts_re,
                        solver_t::real_vector_view_t shifts_im) {
      int np = shifts_re.size();
      if(np == 0) return;

      p.resize(np);
      std::iota(p.begin(), p.end(), 0);
      // p will contain the permutation that puts ritz_bounds in
      // the descending order of magnitude
      std::sort(p.begin(), p.end(),
                [&](int i, int j) { return ritz_bounds(i) > ritz_bounds(j); });
      // Apply permutation p to ritz_values_* and use the result to fill shifts
      for(int i = 0; i < np; ++i) {
        shifts_re(i) = ritz_values_re(p[i]);
        shifts_im(i) = ritz_values_im(p[i]);
      }
    };

    SECTION("Standard eigenproblem") {
      auto Aop = [&](vector_const_view_t in, vector_view_t out) {
        out = A * in;
      };

      solver_t ar(A.rows());

      params_t params(nev, params_t::LargestMagnitude, params_t::Ritz);
      params.random_residual_vector = false;
      params.tolerance = 1e-10;
      set_init_residual_vector(ar);
      ar(Aop, params, shifts_f);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors(ar, A);
      check_basis_vectors(ar);
    }

    SECTION("Generalized eigenproblem: Shift-and-Invert mode (real part)") {
      dcomplex sigma(1.0, -0.1);
#ifdef EIGEN_CAN_MIX_REAL_COMPLEX_EXPR
      decltype(A) op_matrix = ((A - sigma * M).inverse() * M).real();
#else
      decltype(A) op_matrix =
          ((A.cast<dcomplex>() - sigma * M.cast<dcomplex>()).inverse() *
           M.cast<dcomplex>())
              .real();
#endif

      auto op = [&](vector_const_view_t in, vector_view_t out) {
        out = op_matrix * in;
      };
      auto Bop = [&](vector_const_view_t in, vector_view_t out) {
        out = M * in;
      };

      solver_t ar(A.rows());

      params_t params(nev, params_t::LargestMagnitude, params_t::Ritz);
      params.random_residual_vector = false;
      params.tolerance = 1e-10;
      params.sigma = sigma;
      set_init_residual_vector(ar);
      ar(op, Bop, solver_t::ShiftAndInvertReal, params, shifts_f);
      CHECK(ar.nconv() >= nev);
      check_eigenvectors_shift_and_invert(ar, A, M);
      check_basis_vectors(ar, M);
    }
  }
}
