/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2016 I. Krivenko
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once

#include <triqs/test_tools/arrays.hpp>

#include <triqs/utility/c14.hpp>
#include <triqs/arrays/linalg/eigenelements.hpp>

#include "arpack_worker.hpp"

using namespace triqs::arrays;
using namespace triqs::arrays::arpack;
using triqs::arrays::linalg::eigenvalues;

template<operator_kind MKind>
using matrix_t = matrix<std14::conditional_t<MKind==Complex,dcomplex,double>>;

// Fill diagonal
template<operator_kind MKind> void fill_diag(matrix_t<MKind> & M, double coeff) {
 for(int n : range(first_dim(M))) M(n,n) = coeff / (n+1);
}

// Fill off-diagonal elements
template<operator_kind MKind> std14::enable_if_t<MKind==Symmetric,void>
fill_offdiag(matrix_t<MKind> & M, int offdiag_offset, double coeff) {
 int N = first_dim(M);
 for(int n : range(N-offdiag_offset)) {
  M(n,n+offdiag_offset) = coeff;
  M(n+offdiag_offset,n) = coeff;
 }
}

template<operator_kind MKind> std14::enable_if_t<MKind==Asymmetric,void>
fill_offdiag(matrix_t<MKind> & M, int offdiag_offset, double coeff) {
 int N = first_dim(M);
 for(int n : range(N-offdiag_offset)) {
  M(n,n+offdiag_offset) = coeff;
  M(n+offdiag_offset,n) = -coeff;
 }
}

template<operator_kind MKind> std14::enable_if_t<MKind==Complex,void>
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
  EXPECT_CLOSE_ARRAY(A*vec, lambda(i)*vec);
 }
}

template<typename AR, typename MT> void check_eigenvectors(AR const& ar, MT const& A, MT const& M) {
 auto lambda = ar.eigenvalues();
 for(int i : range(lambda.size())) {
  auto vec = ar.eigenvectors()(range(), i);
  EXPECT_CLOSE_ARRAY(A*vec, lambda(i)*M*vec);
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
  EXPECT_CLOSE_ARRAY(A*vec, lambda(i)*M*vec);
 }
}
