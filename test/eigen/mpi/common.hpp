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

#include "ezarpack/mpi/arpack_solver.hpp"

#include "../../mpi_util.hpp"
#include "../common.hpp"

// Compute matrix-vector product u = M v, where blocks of vectors v and u are
// distributed among MPI ranks.
template<operator_kind MKind> class mpi_mat_vec : public mpi_dist_op_base {

  using scalar_t = scalar_t<MKind>;

  // M multiplied by a block of v
  mutable vector<scalar_t> M_v_block;

public:
  mpi_mat_vec(int N, MPI_Comm const& comm)
      : mpi_dist_op_base(N, comm), M_v_block(N) {}

  template<typename MT, typename VT, typename REST>
  void operator()(MT const& M, VT const& v, REST& res) const {
    M_v_block = M.middleCols(local_block_start, local_block_size) * v;

    MPI_Datatype datatype = std::is_same<scalar_t, double>::value
                                ? MPI_DOUBLE
                                : MPI_CXX_DOUBLE_COMPLEX;

    MPI_Reduce_scatter(M_v_block.data(), res.data(), block_sizes.data(),
                       datatype, MPI_SUM, comm);
  }
};

// Compute dot product dot(v1, v2), where blocks of vectors v1 and v2 are
// distributed among MPI ranks.
template<operator_kind MKind> class mpi_dot {
  using scalar_t = scalar_t<MKind>;

  MPI_Comm comm;

public:
  mpi_dot(MPI_Comm const& comm) : comm(comm) {}

  template<typename VT1, typename VT2>
  scalar_t operator()(VT1 const& v1, VT2 const& v2) const {
    scalar_t partial_dot_prod = v1.dot(v2);

    MPI_Datatype datatype = std::is_same<scalar_t, double>::value
                                ? MPI_DOUBLE
                                : MPI_CXX_DOUBLE_COMPLEX;
    scalar_t dot_prod;
    MPI_Allreduce(&partial_dot_prod, &dot_prod, 1, datatype, MPI_SUM, comm);

    return dot_prod;
  }
};

////////////////////////////////////////////////////////////////////////////////

// Check that 'ar' contains the correct solution of a standard eigenproblem
template<operator_kind MKind, typename M>
void check_eigenvectors(mpi::arpack_solver<MKind, eigen_storage> const& ar,
                        M const& A) {
  auto lambda = ar.eigenvalues();
  auto vecs = ar.eigenvectors();
  mpi_mat_vec<MKind> prod(ar.dim(), MPI_COMM_WORLD);
  for(int i = 0; i < lambda.size(); ++i) {
    auto vec = vecs.col(i);
    vector<scalar_t<MKind>> lhs(vec.size());
    prod(A, vec, lhs);
    CHECK_THAT(lhs, IsCloseTo(lambda(i) * vec));
  }
}

// Check that 'ar' contains the correct solution of a generalized eigenproblem
template<operator_kind MKind, typename MT>
void check_eigenvectors(mpi::arpack_solver<MKind, eigen_storage> const& ar,
                        MT const& A,
                        MT const& M) {
  auto lambda = ar.eigenvalues();
  auto vecs = ar.eigenvectors();
  mpi_mat_vec<MKind> prod(ar.dim(), MPI_COMM_WORLD);
  for(int i = 0; i < lambda.size(); ++i) {
    auto vec = vecs.col(i);
    vector<scalar_t<MKind>> lhs(vec.size());
    vector<scalar_t<MKind>> rhs(vec.size());
    prod(A, vec, lhs);
    prod(M, vec, rhs);
    rhs *= lambda(i);
    CHECK_THAT(lhs, IsCloseTo(rhs));
  }
}

////////////////////////////////////////////////////////////////////////////////

// In the real symmetric case, eigenvectors form an orthonormal basis
auto get_basis_vectors(
    mpi::arpack_solver<ezarpack::Symmetric, eigen_storage> const& ar)
    -> decltype(ar.eigenvectors()) {
  return ar.eigenvectors();
}

// Check orthogonality of basis vectors (standard eigenproblem)
template<operator_kind MKind>
void check_basis_vectors(mpi::arpack_solver<MKind, eigen_storage> const& ar) {
  auto vecs = get_basis_vectors(ar);
  mpi_dot<MKind> dot(MPI_COMM_WORLD);
  for(int i = 0; i < vecs.cols(); ++i) {
    auto vi = vecs.col(i);
    for(int j = 0; j < vecs.cols(); ++j) {
      auto vj = vecs.col(j);
      CHECK(std::abs(dot(vi, vj) - double(i == j)) < 1e-10);
    }
  }
}

// Check orthogonality of basis vectors (generalized eigenproblem)
template<operator_kind MKind, typename MT>
void check_basis_vectors(mpi::arpack_solver<MKind, eigen_storage> const& ar,
                         MT const& B) {
  auto vecs = get_basis_vectors(ar);
  mpi_mat_vec<MKind> prod(ar.dim(), MPI_COMM_WORLD);
  mpi_dot<MKind> dot(MPI_COMM_WORLD);
  vector<scalar_t<MKind>> Bvj(ar.local_block_size());
  for(int i = 0; i < vecs.cols(); ++i) {
    auto vi = vecs.col(i);
    for(int j = 0; j < vecs.cols(); ++j) {
      auto vj = vecs.col(j);
      prod(B, vj, Bvj);
      CHECK(std::abs(dot(vi, Bvj) - double(i == j)) < 1e-10);
    }
  }
}
