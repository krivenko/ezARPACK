/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2025 Igor Krivenko
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/
/// @file ezarpack/mpi/mpi_util.hpp
/// @brief MPI utility functions.
#pragma once

#include <mpi.h>

namespace ezarpack {
namespace mpi {

/// Get size of an MPI communicator.
/// @param comm MPI communicator.
inline int size(MPI_Comm const& comm) {
  int comm_size;
  MPI_Comm_size(comm, &comm_size);
  return comm_size;
}

/// Get rank of the calling process within an MPI communicator.
/// @param comm MPI communicator.
inline int rank(MPI_Comm const& comm) {
  int comm_rank;
  MPI_Comm_rank(comm, &comm_rank);
  return comm_rank;
}

/// Compute the size of a vector block owned by an MPI process for the most
/// even distribution of the vector among all processes.
/// @param N size of the full vector.
/// @param comm_size Size of the MPI communicator.
/// @param comm_rank Rank within the MPI communicator.
inline int compute_local_block_size(int N, int comm_size, int comm_rank) {
  int small_block_size = N / comm_size;
  int big_block_size = small_block_size + 1;
  int n_big_blocks = N - small_block_size * comm_size;
  // The first n_big_blocks ranks own blocks of size big_block_size,
  // the rest of the ranks own blocks of size small_block_size.
  return (comm_rank < n_big_blocks) ? big_block_size : small_block_size;
}

/// Compute the starting index of a vector block owned by an MPI process for
/// the most even distribution of the vector among all processes.
/// @param N size of the full vector.
/// @param comm_size Size of the MPI communicator.
/// @param comm_rank Rank within the MPI communicator.
inline int compute_local_block_start(int N, int comm_size, int comm_rank) {
  int small_block_size = N / comm_size;
  int big_block_size = small_block_size + 1;
  int n_big_blocks = N - small_block_size * comm_size;
  // The first n_big_blocks ranks own blocks of size big_block_size,
  // the rest of the ranks own blocks of size small_block_size.
  return (comm_rank < n_big_blocks)
             ? (big_block_size * comm_rank)
             : (big_block_size * n_big_blocks +
                small_block_size * (comm_rank - n_big_blocks));
}

} // namespace mpi
} // namespace ezarpack
