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
#pragma once

#include <algorithm>
#include <cmath>
#include <numeric>
#include <type_traits>
#include <vector>

#include <catch2/catch.hpp>

#include "ezarpack/arpack_solver.hpp"
#include "ezarpack/storages/triqs.hpp"

#include "../common.hpp"

using namespace ezarpack;
using namespace triqs::arrays;

// Make a test sparse matrix
template<operator_kind MKind, typename T = scalar_t<MKind>>
matrix<T> make_sparse_matrix(int N,
                             T diag_coeff_shift,
                             T diag_coeff_amp,
                             int offdiag_offset,
                             T offdiag_coeff) {
  auto refl_offdiag_coeff = reflect_coeff<MKind>(offdiag_coeff);
  matrix<T> M(N, N);
  assign_foreach(M, [&](int i, int j) {
    if(i == j)
      return diag_coeff_amp * T(i % 2) + diag_coeff_shift;
    else if(j - i == offdiag_offset)
      return offdiag_coeff;
    else if(i - j == offdiag_offset)
      return refl_offdiag_coeff;
    else
      return T(0);
  });
  return M;
}

// Make a test inner product matrix
template<operator_kind MKind, typename T = scalar_t<MKind>>
matrix<T> make_inner_prod_matrix(int N) {
  matrix<T> M(N, N);
  assign_foreach(M, [](int i, int j) {
    return (i == j) ? 1.5 : (std::abs(i - j) == 1 ? 0.25 : 0);
  });
  return M;
}

////////////////////////////////////////////////////////////////////////////////

// Catch2 Matcher class that checks proximity of two TRIQS vectors
template<typename Scalar>
class IsCloseToMatcher : public Catch::MatcherBase<vector<Scalar>> {
  vector<Scalar> ref;
  double tol;

public:
  template<typename T>
  IsCloseToMatcher(T&& ref, double tol)
      : ref(triqs::clef::eval(ref)), tol(tol) {}

  virtual bool match(vector<Scalar> const& x) const override {
    return max_element(abs(x - ref)) < tol;
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

// Check that 'ar' contains the correct solution of a standard eigenproblem
template<operator_kind MKind, typename M>
void check_eigenvectors(arpack_solver<MKind, triqs_storage> const& ar,
                        M const& A) {
  auto lambda = ar.eigenvalues();
  auto vecs = ar.eigenvectors();
  for(int i : range(lambda.size())) {
    auto vec = vecs(range(), i);
    CHECK_THAT(A * vec, IsCloseTo(lambda(i) * vec));
  }
}

// Check that 'ar' contains the correct solution of a generalized eigenproblem
template<operator_kind MKind, typename MT>
void check_eigenvectors(arpack_solver<MKind, triqs_storage> const& ar,
                        MT const& A,
                        MT const& M) {
  auto lambda = ar.eigenvalues();
  auto vecs = ar.eigenvectors();
  for(int i : range(lambda.size())) {
    auto vec = vecs(range(), i);
    CHECK_THAT(A * vec, IsCloseTo(lambda(i) * M * vec));
  }
}

// Check that 'ar' contains the correct solution of a generalized eigenproblem
// (Asymmetric Shift-and-Invert modes)
template<typename MT>
void check_eigenvectors_shift_and_invert(
    arpack_solver<Asymmetric, triqs_storage> const& ar,
    MT const& A,
    MT const& M) {
  using solver_t = arpack_solver<Asymmetric, triqs_storage>;
  using vector_view_t = solver_t::vector_view_t;
  using vector_const_view_t = solver_t::vector_const_view_t;
  auto Aop = [&](vector_const_view_t in, vector_view_t out) { out = A * in; };
  auto lambda = ar.eigenvalues(Aop);
  auto vecs = ar.eigenvectors();
  for(int i = 0; i < int(lambda.size()); ++i) {
    auto vec = vecs(range(), i);
    CHECK_THAT(A * vec, IsCloseTo(lambda(i) * M * vec, 1e-9));
  }
}

////////////////////////////////////////////////////////////////////////////////

// In the real symmetric case, eigenvectors form an orthonormal basis
auto get_basis_vectors(arpack_solver<Symmetric, triqs_storage> const& ar)
    -> decltype(ar.eigenvectors()) {
  return ar.eigenvectors();
}
// In the other two cases we must call schur_vectors()
template<operator_kind MKind>
auto get_basis_vectors(arpack_solver<MKind, triqs_storage> const& ar)
    -> decltype(ar.schur_vectors()) {
  return ar.schur_vectors();
}

// Check orthogonality of basis vectors (standard eigenproblem)
template<operator_kind MKind>
void check_basis_vectors(arpack_solver<MKind, triqs_storage> const& ar) {
  auto vecs = get_basis_vectors(ar);
  for(int i : range(second_dim(vecs))) {
    auto vi = vecs(range(), i);
    for(int j : range(second_dim(vecs))) {
      auto vj = vecs(range(), j);
      CHECK(std::abs(dotc(vi, vj) - double(i == j)) < 1e-10);
    }
  }
}
// Check orthogonality of basis vectors (generalized eigenproblem)
template<operator_kind MKind, typename MT>
void check_basis_vectors(arpack_solver<MKind, triqs_storage> const& ar,
                         MT const& B) {
  auto vecs = get_basis_vectors(ar);
  for(int i : range(second_dim(vecs))) {
    auto vi = vecs(range(), i);
    for(int j : range(second_dim(vecs))) {
      auto vj = vecs(range(), j);
      CHECK(std::abs(dotc(vi, B * vj) - double(i == j)) < 1e-10);
    }
  }
}
