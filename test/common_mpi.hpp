/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2024 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/
#pragma once

#include <vector>

#include <catch2/catch.hpp>

#include "ezarpack/mpi/mpi_util.hpp"
#include "ezarpack/mpi/solver_base.hpp"

// A common base for utility classes working with vectors distributed among
// MPI ranks.
class mpi_dist_op_base {
protected:
  MPI_Comm comm;       // MPI communicator
  int const comm_size; // MPI communicator size
  int const comm_rank; // Rank within the MPI communicator

  std::vector<int> const block_sizes; // Sizes of vector blocks

  int const local_block_start; // Start of the local block
  int const local_block_size;  // Size of the local block

public:
  mpi_dist_op_base(int N, MPI_Comm const& comm)
      : comm(comm),
        comm_size(ezarpack::mpi::size(comm)),
        comm_rank(ezarpack::mpi::rank(comm)),
        block_sizes(compute_block_sizes(N, comm_size)),
        local_block_start(
            ezarpack::mpi::compute_local_block_start(N, comm_size, comm_rank)),
        local_block_size(block_sizes[comm_rank]) {}

  static std::vector<int> compute_block_sizes(int N, int comm_size) {
    std::vector<int> sizes(comm_size);
    for(int rank = 0; rank < comm_size; ++rank)
      sizes[rank] = ezarpack::mpi::compute_local_block_size(N, comm_size, rank);
    return sizes;
  }
};

// Set initial residual vector
template<ezarpack::operator_kind OpKind, typename Backend>
void set_init_residual_vector(
    ezarpack::mpi::arpack_solver<OpKind, Backend>& ar) {
  int const N = ar.dim();
  int const block_start = ar.local_block_start();
  int const block_size = ar.local_block_size();
  for(int i = 0; i < block_size; ++i)
    ar.residual_vector()[i] = double(i + block_start) / N;
}

// Test constructor of mpi::arpack_solver
template<typename SolverType> void test_mpi_arpack_solver_ctor() {
  std::vector<std::vector<unsigned int>> block_sizes = {
      {100}, {50, 50}, {34, 33, 33}, {25, 25, 25, 25}};
  std::vector<std::vector<unsigned int>> block_starts = {
      {0}, {0, 50}, {0, 34, 67}, {0, 25, 50, 75}};

  auto comm_size = ezarpack::mpi::size(MPI_COMM_WORLD);
  auto comm_rank = ezarpack::mpi::rank(MPI_COMM_WORLD);

  if(comm_size <= 4) {
    SolverType ar1(100, MPI_COMM_WORLD);
    CHECK(ar1.local_block_size() == block_sizes[comm_size - 1][comm_rank]);
    CHECK(ar1.local_block_start() == block_starts[comm_size - 1][comm_rank]);

    SolverType ar2(block_sizes[comm_size - 1], MPI_COMM_WORLD);
    CHECK(ar2.local_block_size() == block_sizes[comm_size - 1][comm_rank]);
    CHECK(ar2.local_block_start() == block_starts[comm_size - 1][comm_rank]);
  }
}
