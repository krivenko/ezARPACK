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

#include "ezarpack/storages/xtensor.hpp"
#include "ezarpack/arpack_worker.hpp"

#include <xtensor/xio.hpp>
#include <xtensor/xnorm.hpp>
#include <xtensor/xeval.hpp>
#include <xtensor-blas/xlinalg.hpp>

using namespace ezarpack;
using namespace xt;

////////////////////////////////////////////////////////////////////////////////

template<operator_kind MKind>
using scalar_t =
  typename std::conditional<MKind==Complex, dcomplex, double>::type;
template<typename T> using vector = xtensor<T, 1>;
template<typename T> using matrix = xtensor<T, 2, layout_type::column_major>;

template<operator_kind MKind> scalar_t<MKind> reflect_coeff(scalar_t<MKind> x);
template<> double reflect_coeff<ezarpack::Symmetric>(double x) { return x; }
template<> double reflect_coeff<ezarpack::Asymmetric>(double x) { return -x; }
template<> dcomplex reflect_coeff<ezarpack::Complex>(dcomplex x) {
  return std::conj(x);
}

// Make a test sparse matrix
template<operator_kind MKind, typename T = scalar_t<MKind>>
matrix<T> make_sparse_matrix(int N,
                             T diag_coeff,
                             int offdiag_offset,
                             T offdiag_coeff) {
  auto refl_offdiag_coeff = reflect_coeff<MKind>(offdiag_coeff);
  auto M = matrix<T>::from_shape({size_t(N), size_t(N)});
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
matrix<T> make_inner_prod_matrix(int N) {
  auto M = matrix<T>::from_shape({size_t(N), size_t(N)});
  for(int i = 0; i < N; ++i) {
    for(int j = 0; j < N; ++j) {
      M(i, j) = std::exp(-(i-j)*(i-j)/2.0);
    }
  }
  return M;
}

////////////////////////////////////////////////////////////////////////////////

// Catch2 Matcher class that checks proximity of two xtensor vectors
template<typename Scalar>
class IsCloseToMatcher : public Catch::MatcherBase<vector<Scalar>> {
  vector<Scalar> ref;
  double tol;

public:
  template<typename T>
  IsCloseToMatcher(T && ref, double tol) : ref(ref), tol(tol) {}

  virtual bool match(vector<Scalar> const& x) const override {
    return norm_linf(x - ref, {0})[0] < tol;
  }

  virtual std::string describe() const override {
    std::ostringstream ss;
    ss << "is close to " << ref << "(tol = " << tol << ")";
    return ss.str();
  }
};

template<typename T>
inline IsCloseToMatcher<typename T::value_type> IsCloseTo(T && ref,
                                                          double tol = 1e-10) {
  return IsCloseToMatcher<typename T::value_type>(ref, tol);
}

////////////////////////////////////////////////////////////////////////////////

// Check that 'ar' contains the correct solution of a standard eigenproblem
template<typename AR, typename M> void check_eigenvectors(AR const& ar,
                                                          M const& A) {
  auto lambda = ar.eigenvalues();
  auto const vecs = ar.eigenvectors();
  for(int i = 0; i < int(lambda.size()); ++i) {
    auto vec = view(vecs, all(), i);
    CHECK_THAT(linalg::dot(A, vec), IsCloseTo(lambda(i) * vec));
  }
}

// Check that 'ar' contains the correct solution of a generalized eigenproblem
template<typename AR, typename MT> void check_eigenvectors(AR const& ar,
                                                           MT const& A,
                                                           MT const& M) {
  auto lambda = ar.eigenvalues();
  auto const vecs = ar.eigenvectors();
  for(int i = 0; i < int(lambda.size()); ++i) {
    auto vec = view(vecs, all(), i);
    CHECK_THAT(linalg::dot(A, vec), IsCloseTo(lambda(i) * linalg::dot(M, vec)));
  }
}
