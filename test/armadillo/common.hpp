/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2020 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/
#pragma once

#include <algorithm>
#include <cmath>
#include <numeric>
#include <type_traits>
#include <vector>

#include <catch2/catch.hpp>

#include "ezarpack/arpack_worker.hpp"
#include "ezarpack/storages/armadillo.hpp"

using namespace ezarpack;
using namespace arma;

////////////////////////////////////////////////////////////////////////////////

template<operator_kind MKind>
using scalar_t =
    typename std::conditional<MKind == Complex, dcomplex, double>::type;

template<operator_kind MKind> scalar_t<MKind> reflect_coeff(scalar_t<MKind> x);
template<> double reflect_coeff<Symmetric>(double x) { return x; }
template<> double reflect_coeff<Asymmetric>(double x) { return -x; }
template<> dcomplex reflect_coeff<Complex>(dcomplex x) { return -x; }

// Make a test sparse matrix
template<operator_kind MKind, typename T = scalar_t<MKind>>
Mat<T> make_sparse_matrix(int N,
                          T diag_coeff_shift,
                          T diag_coeff_amp,
                          int offdiag_offset,
                          T offdiag_coeff) {
  auto refl_offdiag_coeff = reflect_coeff<MKind>(offdiag_coeff);
  Mat<T> M(N, N);
  for(int i = 0; i < N; ++i) {
    for(int j = 0; j < N; ++j) {
      if(i == j)
        M(i, j) = diag_coeff_amp * T(i % 2) + diag_coeff_shift;
      else if(j - i == offdiag_offset)
        M(i, j) = offdiag_coeff;
      else if(i - j == offdiag_offset)
        M(i, j) = refl_offdiag_coeff;
      else
        M(i, j) = T(0);
    }
  }
  return M;
}

// Make a test inner product matrix
template<operator_kind MKind, typename T = scalar_t<MKind>>
Mat<T> make_inner_prod_matrix(int N) {
  Mat<T> M(N, N);
  for(int i = 0; i < N; ++i) {
    for(int j = 0; j < N; ++j) {
      M(i, j) = (i == j) ? 1.5 : (std::abs(i - j) == 1 ? 0.25 : 0);
    }
  }
  return M;
}

////////////////////////////////////////////////////////////////////////////////

// Catch2 Matcher class that checks proximity of two Armadillo vectors
template<typename Scalar>
class IsCloseToMatcher : public Catch::MatcherBase<Col<Scalar>> {
  Col<Scalar> ref;
  double tol;

public:
  template<typename T>
  IsCloseToMatcher(T&& ref, double tol) : ref(ref), tol(tol) {}

  virtual bool match(Col<Scalar> const& x) const override {
    return max(abs(x - ref)) < tol;
  }

  virtual std::string describe() const override {
    std::ostringstream ss;
    ss << "is close to " << ref << "(tol = " << tol << ")";
    return ss.str();
  }
};

template<typename T>
inline IsCloseToMatcher<typename T::elem_type> IsCloseTo(T&& ref,
                                                         double tol = 1e-10) {
  return IsCloseToMatcher<typename T::elem_type>(ref, tol);
}

////////////////////////////////////////////////////////////////////////////////

// Check that 'ar' contains the correct solution of a standard eigenproblem
template<typename AR, typename M>
void check_eigenvectors(AR const& ar, M const& A) {
  auto lambda = ar.eigenvalues();
  auto vecs = ar.eigenvectors();
  for(int i = 0; i < int(lambda.n_elem); ++i) {
    auto vec = vecs.col(i);
    CHECK_THAT(A * vec, IsCloseTo(lambda[i] * vec));
  }
}

// Check that 'ar' contains the correct solution of a generalized eigenproblem
template<typename AR, typename MT>
void check_eigenvectors(AR const& ar, MT const& A, MT const& M) {
  auto lambda = ar.eigenvalues();
  auto vecs = ar.eigenvectors();
  for(int i = 0; i < int(lambda.n_elem); ++i) {
    auto vec = vecs.col(i);
    CHECK_THAT(A * vec, IsCloseTo(lambda[i] * M * vec));
  }
}

// Check that 'ar' contains the correct solution of a generalized eigenproblem
// (Asymmetric Shift-and-Invert modes)
template<typename MT>
void check_eigenvectors_shift_and_invert(
    arpack_worker<Asymmetric, armadillo_storage> const& ar,
    MT const& A,
    MT const& M) {
  using worker_t = arpack_worker<Asymmetric, armadillo_storage>;
  using vector_view_t = worker_t::vector_view_t;
  using vector_const_view_t = worker_t::vector_const_view_t;
  auto Aop = [&](vector_const_view_t from, vector_view_t to) { to = A * from; };
  auto lambda = ar.eigenvalues(Aop);
  auto vecs = ar.eigenvectors();
  for(int i = 0; i < int(lambda.size()); ++i) {
    auto vec = vecs.col(i);
    CHECK_THAT(A * vec, IsCloseTo(lambda[i] * M * vec, 1e-9));
  }
}

////////////////////////////////////////////////////////////////////////////////

// In the real symmetric case, eigenvectors form an orthonormal basis
auto get_basis_vectors(arpack_worker<Symmetric, armadillo_storage> const& ar)
    -> decltype(ar.eigenvectors()) {
  return ar.eigenvectors();
}
// In the other two cases we must call schur_vectors()
template<typename AR>
auto get_basis_vectors(AR const& ar) -> decltype(ar.schur_vectors()) {
  return ar.schur_vectors();
}

// Check orthogonality of basis vectors (standard eigenproblem)
template<typename AR> void check_basis_vectors(AR const& ar) {
  auto vecs = get_basis_vectors(ar);
  for(int i = 0; i < int(vecs.n_cols); ++i) {
    auto vi = vecs.col(i);
    for(int j = 0; j < int(vecs.n_cols); ++j) {
      auto vj = vecs.col(j);
      CHECK(std::abs(dot(conj(vi), vj) - double(i == j)) < 1e-10);
    }
  }
}
// Check orthogonality of basis vectors (generalized eigenproblem)
template<typename AR, typename MT>
void check_basis_vectors(AR const& ar, MT const& B) {
  auto vecs = get_basis_vectors(ar);
  for(int i = 0; i < int(vecs.n_cols); ++i) {
    auto vi = vecs.col(i);
    for(int j = 0; j < int(vecs.n_cols); ++j) {
      auto vj = vecs.col(j);
      CHECK(std::abs(dot(conj(vi), B * vj) - double(i == j)) < 1e-10);
    }
  }
}
