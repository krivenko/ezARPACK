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
#include <utility>
#include <vector>

#include <catch2/catch.hpp>

#include "ezarpack/arpack_solver.hpp"
#include "ezarpack/storages/raw.hpp"

using namespace ezarpack;

////////////////////////////////////////////////////////////////////////////////

template<operator_kind MKind>
using scalar_t =
    typename std::conditional<MKind == Complex, dcomplex, double>::type;

double conj(double x) { return x; }
dcomplex conj(dcomplex x) { return std::conj(x); }

template<operator_kind MKind> scalar_t<MKind> reflect_coeff(scalar_t<MKind> x);
template<> double reflect_coeff<Symmetric>(double x) { return x; }
template<> double reflect_coeff<Asymmetric>(double x) { return -x; }
template<> dcomplex reflect_coeff<Complex>(dcomplex x) { return -x; }

// Make memory buffer ro accommodate N elements of type T
template<typename T> std::unique_ptr<T[]> make_buffer(int N) {
  return std::unique_ptr<T[]>(new T[N]);
}

// Make a test sparse matrix
template<operator_kind MKind, typename T = scalar_t<MKind>>
std::unique_ptr<T[]> make_sparse_matrix(int N,
                                        T diag_coeff_shift,
                                        T diag_coeff_amp,
                                        int offdiag_offset,
                                        T offdiag_coeff) {
  auto refl_offdiag_coeff = reflect_coeff<MKind>(offdiag_coeff);
  auto M = make_buffer<T>(N * N);
  for(int i = 0; i < N; ++i) {
    for(int j = 0; j < N; ++j) {
      int n = i + j * N;
      if(i == j)
        M[n] = diag_coeff_amp * T(i % 2) + diag_coeff_shift;
      else if(j - i == offdiag_offset)
        M[n] = offdiag_coeff;
      else if(i - j == offdiag_offset)
        M[n] = refl_offdiag_coeff;
      else
        M[n] = T(0);
    }
  }
  return M;
}

// Make a test inner product matrix
template<operator_kind MKind, typename T = scalar_t<MKind>>
std::unique_ptr<T[]> make_inner_prod_matrix(int N) {
  auto M = make_buffer<T>(N * N);
  for(int i = 0; i < N; ++i) {
    for(int j = 0; j < N; ++j) {
      M[i + j * N] = (i == j) ? 1.5 : (std::abs(i - j) == 1 ? 0.25 : 0);
    }
  }
  return M;
}

////////////////////////////////////////////////////////////////////////////////

// Matrix-vector product m * v
template<typename M, typename V, typename O>
void mv_product(M const* m, V const* v, O* out, int N) {
  for(int i = 0; i < N; ++i) {
    out[i] = .0;
    for(int j = 0; j < N; ++j) {
      out[i] += m[i + j * N] * v[j];
    }
  }
}

// Matrix-matrix product m1 * m2
template<typename M1, typename M2, typename O>
void mm_product(M1 const* m1, M2 const* m2, O* out, int N) {
  for(int i = 0; i < N; ++i) {
    for(int j = 0; j < N; ++j) {
      out[i + j * N] = .0;
      for(int k = 0; k < N; ++k) {
        out[i + j * N] += m1[i + k * N] * m2[k + j * N];
      }
    }
  }
}

// Multiplication of a vector by a scalar
template<typename V, typename T, typename O>
void scale(V const* v, T scalar, O* out, int N) {
  for(int i = 0; i < N; ++i) {
    out[i] = v[i] * scalar;
  }
}

// Based on https://en.wikipedia.org/wiki/LU_decomposition
template<typename T> void invert(T const* m, T* out, int N) {

  // Construct LU decomposition
  auto A = make_buffer<T>(N * N);
  std::copy(m, m + N * N, A.get());
  auto P = make_buffer<int>(N + 1);

  double absA;

  for(int i = 0; i <= N; ++i)
    P[i] = i;

  for(int i = 0; i < N; ++i) {
    double maxA = 0.0;
    int imax = i;

    for(int k = i; k < N; ++k) {
      if((absA = std::abs(A[k + i * N])) > maxA) {
        maxA = absA;
        imax = k;
      }
    }
    if(maxA < 1e-14) throw std::runtime_error("invert: matrix is degenerate");

    if(imax != i) {
      // pivoting P
      std::swap(P[i], P[imax]);

      // pivoting rows of A
      for(int j = 0; j < N; ++j)
        std::swap(A[i + N * j], A[imax + N * j]);

      // counting pivots starting from N (for determinant)
      ++P[N];
    }

    for(int j = i + 1; j < N; ++j) {
      A[j + i * N] /= A[i + i * N];

      for(int k = i + 1; k < N; ++k)
        A[j + k * N] -= A[j + i * N] * A[i + k * N];
    }
  }

  // Use the LU decomposition to compute the inverse
  for(int j = 0; j < N; ++j) {
    for(int i = 0; i < N; ++i) {
      if(P[i] == j)
        out[i + j * N] = 1.0;
      else
        out[i + j * N] = 0.0;

      for(int k = 0; k < i; ++k)
        out[i + j * N] -= A[i + k * N] * out[k + j * N];
    }

    for(int i = N - 1; i >= 0; --i) {
      for(int k = i + 1; k < N; ++k)
        out[i + j * N] -= A[i + k * N] * out[k + j * N];

      out[i + j * N] = out[i + j * N] / A[i + i * N];
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

// Catch2 matcher class that checks proximity of two raw memory arrays
template<typename Scalar>
class IsCloseToMatcher : public Catch::MatcherBase<Scalar const*> {
  Scalar const* ref;
  int N;
  double tol;

public:
  IsCloseToMatcher(Scalar const* ref, int N, double tol)
      : ref(ref), N(N), tol(tol) {}

  virtual bool match(Scalar const* const& x) const override {
    double max_diff = 0;
    for(int i = 0; i < N; ++i) {
      max_diff = std::max(max_diff, std::abs(x[i] - ref[i]));
    }
    return max_diff < tol;
  }

  virtual std::string describe() const override {
    std::ostringstream ss;
    ss << "is close to [";
    std::for_each(ref, ref + N, [&](Scalar x) { ss << x << ", "; });
    ss << "] (tol = " << tol << ")";
    return ss.str();
  }
};

template<typename Scalar>
inline IsCloseToMatcher<Scalar>
IsCloseTo(Scalar const* ref, int N, double tol = 1e-10) {
  return IsCloseToMatcher<Scalar>(ref, N, tol);
}

////////////////////////////////////////////////////////////////////////////////

template<typename T> T const* get_ptr(T const* x) { return x; }
template<typename T> T const* get_ptr(std::unique_ptr<T[]> const& x) {
  return x.get();
}

// Check that 'ar' contains the correct solution of a standard eigenproblem
template<
    operator_kind MKind,
    typename M,
    typename T =
        typename std::conditional<MKind == Symmetric, double, dcomplex>::type>
void check_eigenvectors(arpack_solver<MKind, raw_storage> const& ar,
                        M const* m,
                        int N,
                        int nev) {
  auto eigenvalues = ar.eigenvalues();
  auto eigenvectors = ar.eigenvectors();
  for(int i = 0; i < nev; ++i) {
    // RHS
    auto rhs = make_buffer<T>(N);
    scale(get_ptr(eigenvectors) + i * N, eigenvalues[i], rhs.get(), N);
    // LHS
    auto lhs = make_buffer<T>(N);
    mv_product(m, get_ptr(eigenvectors) + i * N, lhs.get(), N);

    CHECK_THAT(rhs.get(), IsCloseTo(lhs.get(), N));
  }
}

// Check that 'ar' contains the correct solution of a generalized eigenproblem
template<
    operator_kind MKind,
    typename M,
    typename T =
        typename std::conditional<MKind == Symmetric, double, dcomplex>::type>
void check_eigenvectors(arpack_solver<MKind, raw_storage> const& ar,
                        M const* a,
                        M const* m,
                        int N,
                        int nev) {
  auto eigenvalues = ar.eigenvalues();
  auto eigenvectors = ar.eigenvectors();
  for(int i = 0; i < nev; ++i) {
    // RHS
    auto rhs = make_buffer<T>(N);
    mv_product(m, get_ptr(eigenvectors) + i * N, rhs.get(), N);
    scale(rhs.get(), eigenvalues[i], rhs.get(), N);
    // LHS
    auto lhs = make_buffer<T>(N);
    mv_product(a, get_ptr(eigenvectors) + i * N, lhs.get(), N);

    CHECK_THAT(rhs.get(), IsCloseTo(lhs.get(), N));
  }
}

// Check that 'ar' contains the correct solution of a generalized eigenproblem
// (Asymmetric Shift-and-Invert modes)
template<typename M>
void check_eigenvectors_shift_and_invert(
    arpack_solver<Asymmetric, raw_storage> const& ar,
    M const* a,
    M const* m,
    int N,
    int nev) {
  using solver_t = arpack_solver<Asymmetric, raw_storage>;
  using vector_view_t = solver_t::vector_view_t;
  using vector_const_view_t = solver_t::vector_const_view_t;
  auto Aop = [&](vector_const_view_t in, vector_view_t out) {
    mv_product(a, in, out, N);
  };
  auto eigenvalues = ar.eigenvalues(Aop);
  auto eigenvectors = ar.eigenvectors();
  for(int i = 0; i < nev; ++i) {
    // RHS
    auto rhs = make_buffer<dcomplex>(N);
    mv_product(m, get_ptr(eigenvectors) + i * N, rhs.get(), N);
    scale(rhs.get(), eigenvalues[i], rhs.get(), N);
    // LHS
    auto lhs = make_buffer<dcomplex>(N);
    mv_product(a, get_ptr(eigenvectors) + i * N, lhs.get(), N);

    CHECK_THAT(rhs.get(), IsCloseTo(lhs.get(), N));
  }
}

////////////////////////////////////////////////////////////////////////////////

// In the real symmetric case, eigenvectors form an orthonormal basis
auto get_basis_vectors(arpack_solver<Symmetric, raw_storage> const& ar)
    -> decltype(ar.eigenvectors()) {
  return ar.eigenvectors();
}
// In the other two cases we must call schur_vectors()
template<typename AR>
auto get_basis_vectors(AR const& ar) -> decltype(ar.schur_vectors()) {
  return ar.schur_vectors();
}

// Check orthogonality of basis vectors (standard eigenproblem)
template<operator_kind MKind>
void check_basis_vectors(arpack_solver<MKind, raw_storage> const& ar,
                         int N,
                         int nev) {
  auto vecs = get_basis_vectors(ar);
  for(int i = 0; i < nev; ++i) {
    auto iptr = get_ptr(vecs) + i * N;
    for(int j = 0; j < nev; ++j) {
      auto jptr = get_ptr(vecs) + j * N;
      scalar_t<MKind> prod = {};
      for(int k = 0; k < N; ++k)
        prod += conj(*(iptr + k)) * *(jptr + k);
      CHECK(std::abs(prod - double(i == j)) < 1e-10);
    }
  }
}
// Check orthogonality of basis vectors (generalized eigenproblem)
template<operator_kind MKind, typename M>
void check_basis_vectors(arpack_solver<MKind, raw_storage> const& ar,
                         M const* b,
                         int N,
                         int nev) {
  auto vecs = get_basis_vectors(ar);
  auto bx = make_buffer<scalar_t<MKind>>(N);
  for(int i = 0; i < nev; ++i) {
    auto iptr = get_ptr(vecs) + i * N;
    for(int j = 0; j < nev; ++j) {
      auto jptr = get_ptr(vecs) + j * N;
      mv_product(b, jptr, bx.get(), N);
      scalar_t<MKind> prod = {};
      for(int k = 0; k < N; ++k)
        prod += conj(*(iptr + k)) * *(bx.get() + k);
      CHECK(std::abs(prod - double(i == j)) < 1e-10);
    }
  }
}
