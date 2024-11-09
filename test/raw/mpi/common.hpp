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

  int const N;

  // M multiplied by a block of v
  mutable std::unique_ptr<scalar_t[]> M_v_block;

public:
  mpi_mat_vec(int N, MPI_Comm const& comm)
      : mpi_dist_op_base(N, comm), N(N), M_v_block(make_buffer<scalar_t>(N)) {}

  template<typename MT, typename VT, typename REST>
  void operator()(MT const* M, VT const* v, REST* res) const {
    for(int i = 0; i < N; ++i) {
      M_v_block[i] = .0;
      for(int j = 0; j < local_block_size; ++j) {
        M_v_block[i] += M[i + (local_block_start + j) * N] * v[j];
      }
    }

    MPI_Datatype datatype = Complex ? MPI_CXX_DOUBLE_COMPLEX : MPI_DOUBLE;

    MPI_Reduce_scatter(M_v_block.get(), res, block_sizes.data(), datatype,
                       MPI_SUM, comm);
  }
};

// Compute dot product dot(v1, v2), where blocks of vectors v1 and v2 are
// distributed among MPI ranks.
template<bool Complex> class mpi_dot : public mpi_dist_op_base {
  using scalar_t = typename std::conditional<Complex, dcomplex, double>::type;

public:
  mpi_dot(int N, MPI_Comm const& comm) : mpi_dist_op_base(N, comm) {}

  template<typename VT1, typename VT2>
  scalar_t operator()(VT1 const* v1, VT2 const* v2) const {

    scalar_t partial_dot_prod = {};
    for(int k = 0; k < local_block_size; ++k)
      partial_dot_prod += conj(*(v1 + k)) * *(v2 + k);

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
void check_eigenvectors(mpi::arpack_solver<MKind, raw_storage> const& ar,
                        M&& m) {
  using scalar_t =
      typename std::conditional<MKind == Symmetric, double, dcomplex>::type;

  auto eigenvalues = ar.eigenvalues();
  auto eigenvectors = ar.eigenvectors();

  constexpr bool const ComplexEigenVecs = (MKind != ezarpack::Symmetric);
  mpi_mat_vec<ComplexEigenVecs> prod(ar.dim(), MPI_COMM_WORLD);

  auto block_size = ar.local_block_size();
  int nev = ar.nconv();

  auto lhs = make_buffer<scalar_t>(block_size);
  auto rhs = make_buffer<scalar_t>(block_size);
  for(int i = 0; i < nev; ++i) {
    // LHS
    prod(get_ptr(m), get_ptr(eigenvectors) + i * block_size, lhs.get());
    // RHS
    scale(get_ptr(eigenvectors) + i * block_size, eigenvalues[i], rhs.get(),
          block_size);

    CHECK_THAT(lhs.get(), IsCloseTo(rhs.get(), block_size));
  }
}

// Check that 'ar' contains the correct solution of a generalized eigenproblem
template<operator_kind MKind, typename M>
void check_eigenvectors(mpi::arpack_solver<MKind, raw_storage> const& ar,
                        M&& a,
                        M&& m) {
  using scalar_t =
      typename std::conditional<MKind == Symmetric, double, dcomplex>::type;

  auto eigenvalues = ar.eigenvalues();
  auto eigenvectors = ar.eigenvectors();

  constexpr bool const ComplexEigenVecs = (MKind != ezarpack::Symmetric);
  mpi_mat_vec<ComplexEigenVecs> prod(ar.dim(), MPI_COMM_WORLD);

  auto block_size = ar.local_block_size();
  int nev = ar.nconv();

  auto lhs = make_buffer<scalar_t>(block_size);
  auto rhs = make_buffer<scalar_t>(block_size);
  for(int i = 0; i < nev; ++i) {
    // LHS
    prod(get_ptr(a), get_ptr(eigenvectors) + i * block_size, lhs.get());
    // RHS
    prod(get_ptr(m), get_ptr(eigenvectors) + i * block_size, rhs.get());
    scale(rhs.get(), eigenvalues[i], rhs.get(), block_size);

    CHECK_THAT(lhs.get(), IsCloseTo(rhs.get(), block_size));
  }
}

// Check that 'ar' contains the correct solution of a generalized eigenproblem
// (Asymmetric Shift-and-Invert modes)
template<typename M>
void check_eigenvectors_shift_and_invert(
    mpi::arpack_solver<Asymmetric, raw_storage> const& ar,
    M&& a,
    M&& m) {
  using solver_t = mpi::arpack_solver<Asymmetric, raw_storage>;
  using vector_view_t = solver_t::vector_view_t;
  using vector_const_view_t = solver_t::vector_const_view_t;

  mpi_mat_vec<false> prod(ar.dim(), MPI_COMM_WORLD);

  auto block_size = ar.local_block_size();
  int nev = ar.nconv();

  auto Aop = [&](vector_const_view_t in, vector_view_t out) {
    prod(get_ptr(a), in, out);
  };
  auto eigenvalues = ar.eigenvalues(Aop);
  auto eigenvectors = ar.eigenvectors();

  mpi_mat_vec<true> prod_complex(ar.dim(), MPI_COMM_WORLD);

  auto lhs = make_buffer<dcomplex>(block_size);
  auto rhs = make_buffer<dcomplex>(block_size);
  for(int i = 0; i < nev; ++i) {
    // LHS
    prod_complex(get_ptr(a), get_ptr(eigenvectors) + i * block_size, lhs.get());
    // RHS
    prod_complex(get_ptr(m), get_ptr(eigenvectors) + i * block_size, rhs.get());
    scale(rhs.get(), eigenvalues[i], rhs.get(), block_size);

    CHECK_THAT(lhs.get(), IsCloseTo(rhs.get(), block_size));
  }
}

////////////////////////////////////////////////////////////////////////////////

// In the real symmetric case, eigenvectors form an orthonormal basis
auto get_basis_vectors(mpi::arpack_solver<Symmetric, raw_storage> const& ar)
    -> decltype(ar.eigenvectors()) {
  return ar.eigenvectors();
}
// In the other two cases we must call schur_vectors()
template<operator_kind MKind>
auto get_basis_vectors(mpi::arpack_solver<MKind, raw_storage> const& ar)
    -> decltype(ar.schur_vectors()) {
  return ar.schur_vectors();
}

// Check orthogonality of basis vectors (standard eigenproblem)
template<operator_kind MKind>
void check_basis_vectors(mpi::arpack_solver<MKind, raw_storage> const& ar) {
  auto vecs = get_basis_vectors(ar);

  constexpr bool const ComplexBasisVecs = MKind == ezarpack::Complex;
  mpi_dot<ComplexBasisVecs> dot(ar.dim(), MPI_COMM_WORLD);

  auto block_size = ar.local_block_size();
  int nev = ar.nconv();

  for(int i = 0; i < nev; ++i) {
    auto iptr = get_ptr(vecs) + i * block_size;
    for(int j = 0; j < nev; ++j) {
      auto jptr = get_ptr(vecs) + j * block_size;
      CHECK(std::abs(dot(iptr, jptr) - double(i == j)) < 1e-10);
    }
  }
}
// Check orthogonality of basis vectors (generalized eigenproblem)
template<operator_kind MKind, typename M>
void check_basis_vectors(mpi::arpack_solver<MKind, raw_storage> const& ar,
                         M&& b) {
  auto vecs = get_basis_vectors(ar);

  constexpr bool const ComplexBasisVecs = MKind == ezarpack::Complex;
  mpi_mat_vec<ComplexBasisVecs> prod(ar.dim(), MPI_COMM_WORLD);
  mpi_dot<ComplexBasisVecs> dot(ar.dim(), MPI_COMM_WORLD);

  auto block_size = ar.local_block_size();
  int nev = ar.nconv();

  auto bx = make_buffer<scalar_t<MKind>>(block_size);
  for(int i = 0; i < nev; ++i) {
    auto iptr = get_ptr(vecs) + i * block_size;
    for(int j = 0; j < nev; ++j) {
      auto jptr = get_ptr(vecs) + j * block_size;
      prod(get_ptr(b), jptr, bx.get());
      CHECK(std::abs(dot(iptr, bx.get()) - double(i == j)) < 1e-10);
    }
  }
}

} // namespace mpi
} // namespace ezarpack
