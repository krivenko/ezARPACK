/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2025 Igor Krivenko
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

#include "ezarpack/arpack_solver.hpp"
#include "ezarpack/storages/xtensor.hpp"

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/core/xeval.hpp>
#include <xtensor/io/xio.hpp>
#include <xtensor/reducers/xnorm.hpp>

#include "../common.hpp"
#include "../tests.hpp"

using namespace ezarpack;
using namespace xt;

template<typename T> using vector = xtensor<T, 1>;
template<typename T> using matrix = xtensor<T, 2, layout_type::column_major>;

// Make a test sparse matrix
template<operator_kind MKind, typename T = scalar_t<MKind>>
matrix<T> make_sparse_matrix(int N,
                             T diag_coeff_mean,
                             int offdiag_offset,
                             T offdiag_coeff_mean,
                             T offdiag_coeff_diff) {
  auto M = matrix<T>::from_shape({size_t(N), size_t(N)});
  for(int i = 0; i < N; ++i) {
    for(int j = 0; j < N; ++j) {
      if(i == j)
        M(i, j) = diag_coeff_mean;
      else if(j - i == offdiag_offset)
        M(i, j) = offdiag_coeff_mean + offdiag_coeff_diff;
      else if(i - j == offdiag_offset)
        M(i, j) = offdiag_coeff_mean - offdiag_coeff_diff;
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
      M(i, j) = (i == j) ? 1.0 : (std::abs(i - j) == 1 ? 0.1 : 0);
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
  IsCloseToMatcher(T&& ref, double tol) : ref(ref), tol(tol) {}

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
inline IsCloseToMatcher<typename std::remove_reference<T>::type::value_type>
IsCloseTo(T&& ref, double tol = 1e-10) {
  return IsCloseToMatcher<typename std::remove_reference<T>::type::value_type>(
      ref, tol);
}

////////////////////////////////////////////////////////////////////////////////

namespace ezarpack { // To make ADL for the following functions work

// Check that 'ar' contains the correct solution of a standard eigenproblem
template<operator_kind MKind, typename M>
void check_eigenvectors(arpack_solver<MKind, xtensor_storage> const& ar,
                        M const& A) {
  auto lambda = ar.eigenvalues();
  auto const vecs = ar.eigenvectors();
  for(int i = 0; i < int(lambda.size()); ++i) {
    auto vec = view(vecs, all(), i);
    CHECK_THAT(linalg::dot(A, vec), IsCloseTo(lambda(i) * vec));
  }
}

// Check that 'ar' contains the correct solution of a generalized eigenproblem
template<operator_kind MKind, typename MT>
void check_eigenvectors(arpack_solver<MKind, xtensor_storage> const& ar,
                        MT const& A,
                        MT const& M) {
  auto lambda = ar.eigenvalues();
  auto const vecs = ar.eigenvectors();
  for(int i = 0; i < int(lambda.size()); ++i) {
    auto vec = view(vecs, all(), i);
    using linalg::dot;
    CHECK_THAT(dot(A, vec), IsCloseTo(lambda(i) * dot(M, vec)));
  }
}

// Check that 'ar' contains the correct solution of a generalized eigenproblem
// (Asymmetric Shift-and-Invert modes)
template<typename MT>
void check_eigenvectors_shift_and_invert(
    arpack_solver<Asymmetric, xtensor_storage> const& ar,
    MT const& A,
    MT const& M) {
  using solver_t = arpack_solver<Asymmetric, xtensor_storage>;
  using vector_view_t = solver_t::vector_view_t;
  using vector_const_view_t = solver_t::vector_const_view_t;
  using linalg::dot;
  auto Aop = [&](vector_const_view_t in, vector_view_t out) {
    out = dot(A, in);
  };
  auto lambda = ar.eigenvalues(Aop);
  auto vecs = ar.eigenvectors();
  for(int i = 0; i < int(lambda.size()); ++i) {
    auto vec = view(vecs, all(), i);
    CHECK_THAT(dot(A, vec), IsCloseTo(lambda(i) * dot(M, vec), 1e-9));
  }
}

////////////////////////////////////////////////////////////////////////////////

// In the real symmetric case, eigenvectors form an orthonormal basis
auto get_basis_vectors(arpack_solver<Symmetric, xtensor_storage> const& ar)
    -> decltype(ar.eigenvectors()) {
  return ar.eigenvectors();
}
// In the other two cases we must call schur_vectors()
template<operator_kind MKind>
auto get_basis_vectors(arpack_solver<MKind, xtensor_storage> const& ar)
    -> decltype(ar.schur_vectors()) {
  return ar.schur_vectors();
}

// Check orthogonality of basis vectors (standard eigenproblem)
template<operator_kind MKind>
void check_basis_vectors(arpack_solver<MKind, xtensor_storage> const& ar) {
  auto vecs = get_basis_vectors(ar);
  for(int i = 0; i < int(vecs.shape(1)); ++i) {
    auto vi = view(vecs, all(), i);
    for(int j = 0; j < int(vecs.shape(1)); ++j) {
      auto vj = view(vecs, all(), j);
      CHECK(std::abs(linalg::vdot(vi, vj) - double(i == j)) < 1e-10);
    }
  }
}
// Check orthogonality of basis vectors (generalized eigenproblem)
template<operator_kind MKind, typename MT>
void check_basis_vectors(arpack_solver<MKind, xtensor_storage> const& ar,
                         MT const& B) {
  auto vecs = get_basis_vectors(ar);
  for(int i = 0; i < int(vecs.shape(1)); ++i) {
    auto vi = view(vecs, all(), i);
    for(int j = 0; j < int(vecs.shape(1)); ++j) {
      auto vj = view(vecs, all(), j);
      CHECK(std::abs(linalg::vdot(vi, linalg::dot(B, vj)) - double(i == j)) <
            1e-10);
    }
  }
}

} // namespace ezarpack
