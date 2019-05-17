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
#pragma once

#include <type_traits>

#include <catch2/catch.hpp>

#include "ezarpack/storages/armadillo.hpp"
#include "ezarpack/arpack_worker.hpp"

using namespace ezarpack;
using namespace arma;

//////////////////////////////////////////////////////////////////////////////////////////

template<operator_kind MKind>
using scalar_t = typename std::conditional<MKind==Complex, dcomplex, double>::type;

template<operator_kind MKind> scalar_t<MKind> reflect_coeff(scalar_t<MKind> x);
template<> double reflect_coeff<Symmetric>(double x) { return x; }
template<> double reflect_coeff<Asymmetric>(double x) { return -x; }
template<> dcomplex reflect_coeff<Complex>(dcomplex x) { return std::conj(x); }

// Make a test sparse matrix
template<operator_kind MKind, typename T = scalar_t<MKind>>
Mat<T> make_sparse_matrix(int N, T diag_coeff, int offdiag_offset, T offdiag_coeff) {
  auto refl_offdiag_coeff = reflect_coeff<MKind>(offdiag_coeff);
  Mat<T> M(N, N);
  for(int i = 0; i < N; ++i) {
    for(int j = 0; j < N; ++j) {
      if(i == j)
        M(i, j) = diag_coeff / T(i + 1);
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
      M(i, j) = std::exp(-(i-j)*(i-j)/2.0);
    }
  }
  return M;
}

//////////////////////////////////////////////////////////////////////////////////////////

// Catch2 Matcher class that checks proximity of two Armadillo vectors
template<typename Scalar>
class IsCloseToMatcher : public Catch::MatcherBase<Col<Scalar>> {
  Col<Scalar> ref;
  double tol;

public:
  template<typename T>
  IsCloseToMatcher(T && ref, double tol) : ref(ref), tol(tol) {}

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
inline IsCloseToMatcher<typename T::elem_type> IsCloseTo(T && ref, double tol = 1e-10) {
  return IsCloseToMatcher<typename T::elem_type>(ref, tol);
}

//////////////////////////////////////////////////////////////////////////////////////////

// Check that 'ar' contains the correct solution of a standard eigenproblem
template<typename AR, typename M> void check_eigenvectors(AR const& ar, M const& A) {
  auto lambda = ar.eigenvalues();
  auto vecs = ar.eigenvectors();
  for(int i = 0; i < int(lambda.n_elem); ++i) {
    auto vec = vecs.col(i);
    CHECK_THAT(A * vec, IsCloseTo(lambda[i] * vec));
  }
}

// Check that 'ar' contains the correct solution of a generalized eigenproblem
template<typename AR, typename MT> void check_eigenvectors(AR const& ar, MT const& A, MT const& M) {
  auto lambda = ar.eigenvalues();
  auto vecs = ar.eigenvectors();
  for(int i = 0; i < int(lambda.n_elem); ++i) {
    auto vec = vecs.col(i);
    CHECK_THAT(A * vec, IsCloseTo(lambda[i] * M * vec));
  }
}
