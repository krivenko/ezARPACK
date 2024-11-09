/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2024 Igor Krivenko
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/
#pragma once

#include "ezarpack/mpi/arpack_solver.hpp"

#include "../../common_mpi.hpp"
#include "../common.hpp"

// Compute matrix-vector product u = M v, where blocks of vectors v and u are
// distributed among MPI ranks.
template<bool Complex> class mpi_mat_vec : public mpi_dist_op_base {
  using scalar_t = typename std::conditional<Complex, dcomplex, double>::type;

  // M multiplied by a block of v
  mutable vector<scalar_t> M_v_block;

public:
  mpi_mat_vec(int N, MPI_Comm const& comm)
      : mpi_dist_op_base(N, comm), M_v_block(N) {}

  template<typename MT, typename VT, typename REST>
  void operator()(MT const& M, VT const& v, REST& res) const {
    M_v_block = M.middleCols(local_block_start, local_block_size) * v;

    MPI_Datatype datatype = Complex ? MPI_CXX_DOUBLE_COMPLEX : MPI_DOUBLE;

    MPI_Reduce_scatter(M_v_block.data(), res.data(), block_sizes.data(),
                       datatype, MPI_SUM, comm);
  }
};

// Compute dot product dot(v1, v2), where blocks of vectors v1 and v2 are
// distributed among MPI ranks.
template<bool Complex> class mpi_dot {
  using scalar_t = typename std::conditional<Complex, dcomplex, double>::type;

  MPI_Comm comm;

public:
  mpi_dot(MPI_Comm const& comm) : comm(comm) {}

  template<typename VT1, typename VT2>
  scalar_t operator()(VT1 const& v1, VT2 const& v2) const {
    scalar_t partial_dot_prod = v1.dot(v2);

    MPI_Datatype datatype = Complex ? MPI_CXX_DOUBLE_COMPLEX : MPI_DOUBLE;
    scalar_t dot_prod;
    MPI_Allreduce(&partial_dot_prod, &dot_prod, 1, datatype, MPI_SUM, comm);

    return dot_prod;
  }
};

////////////////////////////////////////////////////////////////////////////////

// To make ADL for the following functions work
namespace ezarpack {
namespace mpi {

// Check that 'ar' contains the correct solution of a standard eigenproblem
template<operator_kind MKind, typename M>
void check_eigenvectors(mpi::arpack_solver<MKind, eigen_storage> const& ar,
                        M const& A) {
  auto lambda = ar.eigenvalues();
  auto vecs = ar.eigenvectors();

  constexpr bool const ComplexEigenVecs = (MKind != ezarpack::Symmetric);
  mpi_mat_vec<ComplexEigenVecs> prod(ar.dim(), MPI_COMM_WORLD);

  vector<typename decltype(vecs)::value_type> lhs(ar.local_block_size());
  for(int i = 0; i < lambda.size(); ++i) {
    auto vec = vecs.col(i);
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

  constexpr bool const ComplexEigenVecs = (MKind != ezarpack::Symmetric);
  mpi_mat_vec<ComplexEigenVecs> prod(ar.dim(), MPI_COMM_WORLD);

  vector<typename decltype(vecs)::value_type> lhs(ar.local_block_size());
  vector<typename decltype(vecs)::value_type> rhs(ar.local_block_size());
  for(int i = 0; i < lambda.size(); ++i) {
    auto vec = vecs.col(i);
    prod(A, vec, lhs);
    prod(M, vec, rhs);
    rhs *= lambda(i);
    CHECK_THAT(lhs, IsCloseTo(rhs));
  }
}

// Check that 'ar' contains the correct solution of a generalized eigenproblem
// (Asymmetric Shift-and-Invert modes)
template<typename MT>
void check_eigenvectors_shift_and_invert(
    mpi::arpack_solver<ezarpack::Asymmetric, eigen_storage> const& ar,
    MT const& A,
    MT const& M) {
  using solver_t = arpack_solver<ezarpack::Asymmetric, eigen_storage>;
  using vector_view_t = solver_t::vector_view_t;
  using vector_const_view_t = solver_t::vector_const_view_t;

  mpi_mat_vec<false> prod(ar.dim(), MPI_COMM_WORLD);

  auto Aop = [&](vector_const_view_t in, vector_view_t out) {
    prod(A, in, out);
  };
  auto lambda = ar.eigenvalues(Aop);
  auto vecs = ar.eigenvectors();

  mpi_mat_vec<true> prod_complex(ar.dim(), MPI_COMM_WORLD);

  vector<dcomplex> lhs(ar.local_block_size());
  vector<dcomplex> rhs(ar.local_block_size());
  for(int i = 0; i < lambda.size(); ++i) {
    auto vec = vecs.col(i);
    prod_complex(A, vec, lhs);
    prod_complex(M, vec, rhs);
    rhs *= lambda(i);
    CHECK_THAT(lhs, IsCloseTo(rhs, 1e-9));
  }
}

////////////////////////////////////////////////////////////////////////////////

// In the real symmetric case, eigenvectors form an orthonormal basis
auto get_basis_vectors(
    mpi::arpack_solver<ezarpack::Symmetric, eigen_storage> const& ar)
    -> decltype(ar.eigenvectors()) {
  return ar.eigenvectors();
}
// In the other two cases we must call schur_vectors()
template<operator_kind MKind>
auto get_basis_vectors(mpi::arpack_solver<MKind, eigen_storage> const& ar)
    -> decltype(ar.schur_vectors()) {
  return ar.schur_vectors();
}

// Check orthogonality of basis vectors (standard eigenproblem)
template<operator_kind MKind>
void check_basis_vectors(mpi::arpack_solver<MKind, eigen_storage> const& ar) {
  auto vecs = get_basis_vectors(ar);

  constexpr bool const ComplexBasisVecs = MKind == ezarpack::Complex;
  mpi_dot<ComplexBasisVecs> dot(MPI_COMM_WORLD);

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

  constexpr bool const ComplexBasisVecs = MKind == ezarpack::Complex;
  mpi_mat_vec<ComplexBasisVecs> prod(ar.dim(), MPI_COMM_WORLD);
  mpi_dot<ComplexBasisVecs> dot(MPI_COMM_WORLD);

  vector<typename decltype(vecs)::value_type> Bvj(ar.local_block_size());
  for(int i = 0; i < vecs.cols(); ++i) {
    auto vi = vecs.col(i);
    for(int j = 0; j < vecs.cols(); ++j) {
      auto vj = vecs.col(j);
      prod(B, vj, Bvj);
      CHECK(std::abs(dot(vi, Bvj) - double(i == j)) < 1e-10);
    }
  }
}

} // namespace mpi
} // namespace ezarpack
