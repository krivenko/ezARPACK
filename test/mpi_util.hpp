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

#include <vector>

#include "ezarpack/mpi/mpi_util.hpp"

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
