/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2018 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/
#pragma once

#include <type_traits>

#include <triqs/test_tools/arrays.hpp>

#include <triqs/arrays/linalg/eigenelements.hpp>

#include "ezarpack/storages/triqs.hpp"
#include "ezarpack/arpack_worker.hpp"

using namespace ezarpack;
using namespace triqs::arrays;
using triqs::arrays::linalg::eigenvalues;

template<operator_kind MKind>
using matrix_t = matrix<typename std::conditional<MKind==Complex,dcomplex,double>::type>;

// Fill diagonal
template<operator_kind MKind> void fill_diag(matrix_t<MKind> & M, double coeff) {
 for(int n : range(first_dim(M))) M(n,n) = coeff / (n+1);
}

// Fill off-diagonal elements
template<operator_kind MKind> typename std::enable_if<MKind==Symmetric,void>::type
fill_offdiag(matrix_t<MKind> & M, int offdiag_offset, double coeff) {
 int N = first_dim(M);
 for(int n : range(N-offdiag_offset)) {
  M(n,n+offdiag_offset) = coeff;
  M(n+offdiag_offset,n) = coeff;
 }
}

template<operator_kind MKind> typename std::enable_if<MKind==Asymmetric,void>::type
fill_offdiag(matrix_t<MKind> & M, int offdiag_offset, double coeff) {
 int N = first_dim(M);
 for(int n : range(N-offdiag_offset)) {
  M(n,n+offdiag_offset) = coeff;
  M(n+offdiag_offset,n) = -coeff;
 }
}

template<operator_kind MKind> typename std::enable_if<MKind==Complex,void>::type
fill_offdiag(matrix_t<MKind> & M, int offdiag_offset, dcomplex coeff) {
 int N = first_dim(M);
 for(int n : range(N-offdiag_offset)) {
  M(n,n+offdiag_offset) = coeff;
  M(n+offdiag_offset,n) = std::conj(coeff);
 }
}

template<operator_kind MKind> matrix_t<MKind> make_sparse_matrix(int N, double diag_coeff,
                          int offdiag_offset, typename matrix_t<MKind>::value_type offdiag_coeff) {
 matrix_t<MKind> M(N,N);
 M() = 0;
 fill_diag   <MKind>(M, diag_coeff);
 fill_offdiag<MKind>(M, offdiag_offset, offdiag_coeff);
 return M;
}

template<operator_kind MKind> matrix_t<MKind> make_inner_prod_matrix(int N) {
 matrix_t<MKind> M(N,N);
 assign_foreach(M,[](int i, int j){ return std::exp(-(i-j)*(i-j)/2.0); });
 return M;
}

template<typename AR, typename M> void check_eigenvectors(AR const& ar, M const& A) {
 auto lambda = ar.eigenvalues();
 for(int i : range(lambda.size())) {
  auto vec = ar.eigenvectors()(range(), i);
  EXPECT_ARRAY_NEAR(A*vec, lambda(i)*vec);
 }
}

template<typename AR, typename MT> void check_eigenvectors(AR const& ar, MT const& A, MT const& M) {
 auto lambda = ar.eigenvalues();
 for(int i : range(lambda.size())) {
  auto vec = ar.eigenvectors()(range(), i);
  EXPECT_ARRAY_NEAR(A*vec, lambda(i)*M*vec);
 }
}

template<typename AR, typename MT>
void check_eigenvectors_rayleigh(AR const& ar, MT const& A, MT const& M) {
 auto Aop = [&A](vector_const_view<double> from, vector_view<double> to) {
  to = A * from;
 };
 auto lambda = ar.eigenvalues(Aop);
 for(int i : range(lambda.size())) {
  auto vec = ar.eigenvectors()(range(), i);
  EXPECT_ARRAY_NEAR(A*vec, lambda(i)*M*vec);
 }
}
