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

#include "common.hpp"

//////////////////////////////////////////////
// Eigenproblems with general real matrices //
//////////////////////////////////////////////

using params_t = arpack_worker<Asymmetric>::params_t;

const int N = 100;
const double diag_coeff = 0.75;
const int offdiag_offset = 3;
const double offdiag_coeff = 1.0;
const int nev = 10;

// Symmetric matrix A
auto A = make_sparse_matrix<Asymmetric>(N, diag_coeff, offdiag_offset, offdiag_coeff);
// Inner product matrix
auto M = make_inner_prod_matrix<Asymmetric>(N);

auto spectrum_parts = {params_t::LargestMagnitude,
                       params_t::SmallestMagnitude,
                       params_t::LargestReal, params_t::SmallestReal,
                       params_t::LargestImag, params_t::SmallestImag };

TEST(arpack_worker_symmetric, InnerProduct) {
 ASSERT_GT(eigenvalues(M)(0),.0);
}

// Standard eigenproblem
TEST(arpack_worker_asymmetric, Standard) {
 auto Aop = [](vector_const_view<double> from, int, vector_view<double> to, int) {
  to = A*from;
 };

 arpack_worker<Asymmetric> ar(first_dim(A));

 for(auto e : spectrum_parts) {
  params_t params(nev, e, params_t::Ritz);
  ar(Aop, params);
  check_eigenvectors(ar,A);
 }
}

// Generalized eigenproblem: invert mode
TEST(arpack_worker_asymmetric, Invert) {
 decltype(A) invMA = inverse(M) * A;

 auto op = [&invMA](vector_const_view<double> from, int, vector_view<double> to, int, bool) {
  to = invMA * from;
 };
 auto Bop = [](vector_const_view<double> from, int, vector_view<double> to, int) {
  to = M * from;
 };

 arpack_worker<Asymmetric> ar(first_dim(A));

 for(auto e : spectrum_parts) {
  params_t params(nev, e, params_t::Ritz);
  params.ncv = 30;
  ar(op, Bop, arpack_worker<Asymmetric>::Invert, params);
  check_eigenvectors(ar,A,M);
 }
}

MAKE_MAIN;
