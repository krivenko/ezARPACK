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

#include "arpack.hpp"

using params_t = arpack_worker<Symmetric>::params_t;

const int N = 100;
const double diag_coeff = 0.75;
const int offdiag_offset = 3;
const double offdiag_coeff = 1.0;
const int nev = 10;

// Symmetric matrix A
auto A = make_sparse_matrix<Symmetric>(N, diag_coeff, offdiag_offset, offdiag_coeff);
// Inner product matrix
auto M = make_inner_prod_matrix<Symmetric>(N);

auto spectrum_parts = {params_t::Smallest, params_t::Largest,
                       params_t::SmallestMagnitude, params_t::LargestMagnitude,
                       params_t::BothEnds};

TEST(arpack_worker_symmetric, InnerProduct) {
 ASSERT_GT(eigenvalues(M)(0),.0);
}

TEST(arpack_worker_symmetric, Standard) {
 auto Aop = [](vector_const_view<double> from, int, vector_view<double> to, int) {
  to = A*from;
 };

 arpack_worker<Symmetric> ar(first_dim(A));

 for(auto e : spectrum_parts) {
  params_t params(nev, e, true);
  ar(Aop, params);
  check_eigenvectors(ar,A);
 }
}

TEST(arpack_worker_symmetric, Invert) {
 decltype(A) invM = inverse(M);

 auto op = [&invM](vector_view<double> from, int, vector_view<double> to, int, bool) {
  from = A * from;
  to = invM * from;
 };
 auto Bop = [](vector_const_view<double> from, int, vector_view<double> to, int) {
  to = M * from;
 };

 arpack_worker<Symmetric> ar(first_dim(A));

 for(auto e : spectrum_parts) {
  params_t params(nev, e, true);
  ar(op, Bop, arpack_worker<Symmetric>::Invert, params);
  check_eigenvectors(ar,A,M);
 }
}

TEST(arpack_worker_symmetric, ShiftAndInvert) {
 double sigma = 1.0;
 decltype(A) inv = inverse(A - sigma*M) * M;

 auto op = [&inv](vector_view<double> from, int, vector_view<double> to, int, bool) {
  to = inv * from;
 };
 auto Bop = [](vector_const_view<double> from, int, vector_view<double> to, int) {
  to = M * from;
 };

 arpack_worker<Symmetric> ar(first_dim(A));

 for(auto e : spectrum_parts) {
  params_t params(nev, e, true);
  params.sigma = sigma;
  ar(op, Bop, arpack_worker<Symmetric>::ShiftAndInvert, params);
  check_eigenvectors(ar,A,M);
 }
}

TEST(arpack_worker_symmetric, Buckling) {
 double sigma = 1.0;
 decltype(A) inv = inverse(A - sigma*M) * A;

 auto op = [&inv](vector_view<double> from, int, vector_view<double> to, int, bool) {
  to = inv * from;
 };
 auto Bop = [](vector_const_view<double> from, int, vector_view<double> to, int) {
  to = M * from;
 };

 arpack_worker<Symmetric> ar(first_dim(A));

 for(auto e : spectrum_parts) {
  params_t params(nev, e, true);
  params.sigma = sigma;
  ar(op, Bop, arpack_worker<Symmetric>::Buckling, params);
  check_eigenvectors(ar,A,M);
 }
}

TEST(arpack_worker_symmetric, Cayley) {
 double sigma = 1.0;
 decltype(A) inv = inverse(A - sigma*M) * (A + sigma*M);

 auto op = [&inv](vector_view<double> from, int, vector_view<double> to, int, bool) {
  to = inv * from;
 };
 auto Bop = [&](vector_const_view<double> from, int, vector_view<double> to, int) {
  to = M * from;
 };

 arpack_worker<Symmetric> ar(first_dim(A));

 for(auto e : spectrum_parts) {
  params_t params(nev, e, true);
  params.sigma = sigma;
  ar(op, Bop, arpack_worker<Symmetric>::Cayley, params);
  check_eigenvectors(ar,A,M);
 }
}

MAKE_MAIN;
