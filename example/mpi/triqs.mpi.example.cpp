/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2026 Igor Krivenko
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/

#include <cmath>
#include <iostream>
#include <vector>

// This example shows how to use an MPI-parallelized solver of ezARPACK and
// the TRIQS storage backend to partially diagonalize a large sparse symmetric
// matrix and find a number of its low-lying eigenvalues.

#include <ezarpack/mpi/arpack_solver.hpp>
#include <ezarpack/storages/triqs.hpp>
#include <ezarpack/version.hpp>

using namespace ezarpack;
using namespace triqs::arrays;

// Size of the matrix
const int N = 10000;

// We are going to use a band matrix with this bandwidth
const int bandwidth = 5;

// The number of low-lying eigenvalues we want to compute
const int N_ev = 10;

int main(int argc, char* argv[]) {

  // Initialize MPI environment
  MPI_Init(&argc, &argv);

  // Call utility functions from namespace 'ezarpack::mpi' to find out
  // the world communicator size and the rank of the calling process.
  const int comm_size = ezarpack::mpi::size(MPI_COMM_WORLD);
  const int comm_rank = ezarpack::mpi::rank(MPI_COMM_WORLD);

  // Print ezARPACK version
  if(comm_rank == 0)
    std::cout << "Using ezARPACK version " << EZARPACK_VERSION << std::endl;

  // Construct an MPI-parallelized solver object for the symmetric case.
  // For the TRIQS storage backend, other options would be
  // * `mpi::arpack_solver<Asymmetric, triqs_storage>' for general
  //   real matrices;
  // * `mpi::arpack_solver<Complex, triqs_storage>' for general
  //   complex matrices.
  using solver_t = ezarpack::mpi::arpack_solver<Symmetric, triqs_storage>;
  solver_t solver(N, MPI_COMM_WORLD);

  // Specify parameters for the solver
  using params_t = solver_t::params_t;
  params_t params(N_ev,               // Number of low-lying eigenvalues
                  params_t::Smallest, // We want the smallest eigenvalues
                  true);              // Yes, we want the eigenvectors
                                      // (Ritz vectors) as well

  // Vectors from the N-dimensional space of the problem are partitioned
  // into contiguous blocks. These blocks are distributed among all
  // MPI processes in the communicator used to construct 'solver'.
  int block_start = solver.local_block_start();
  int block_size = solver.local_block_size();
  // Block owned by the calling process covers the index range
  // [block_start; block_start + block_size] within a full vector.

  // Compute and collect sizes of all rank-local blocks for later use.
  std::vector<int> block_sizes(comm_size);
  for(int rank = 0; rank < comm_size; ++rank)
    block_sizes[rank] =
        ezarpack::mpi::compute_local_block_size(N, comm_size, rank);

  // Temporary vector used in distributed matrix-vector multiplication
  vector<double> local_op_in(N);

  // Linear operator representing multiplication of a given vector by our matrix
  // The matrix to be diagonalized is defined as
  // A_{ij} = |i-j| / (1 + i + j), if |i-j| <= bandwidth, zero otherwise
  auto matrix_op = [&](auto in, auto out) {
    // 'in' and 'out' are views of the locally stored blocks of their respective
    // distributed N-dimensional vectors. Therefore, matrix-vector
    // multiplication has to be performed in two steps.

    // 1. Local multiplication of A's columns
    // [block_start; block_start + block_size] by 'in'. The result is an
    // N-dimensional vector stored in 'local_op_in'.
    local_op_in() = 0;
    for(int i = 0; i < N; ++i) {
      int j_min = std::max(block_start, i - bandwidth);
      int j_max = std::min(block_start + block_size - 1, i + bandwidth);
      for(int j = j_min; j <= j_max; ++j) {
        int j_local = j - block_start;
        local_op_in(i) += double(std::abs(i - j)) / (1 + i + j) * in(j_local);
      }
    }

    // 2. Sum up (MPI reduce) results from step 1 and scatter the sum over
    // 'out' blocks stored on different MPI ranks.
    MPI_Reduce_scatter(local_op_in.data_start(), out.data_start(),
                       block_sizes.data(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  };

  // Run diagonalization!
  solver(matrix_op, params);

  if(comm_rank == 0) {
    // Number of converged eigenvalues
    std::cout << solver.nconv() << " out of " << params.n_eigenvalues
              << " eigenvalues have converged" << std::endl;

    // Print found eigenvalues
    std::cout << "Eigenvalues (Ritz values):" << std::endl;
    std::cout << solver.eigenvalues() << std::endl;
  }

  // Check A*v = \lambda*v
  auto const& lambda = solver.eigenvalues();
  auto const& v = solver.eigenvectors();
  vector<double> lhs(block_size), rhs(block_size);

  for(int i = 0; i < N_ev; ++i) { // For each eigenpair ...
    auto const eigenvec = v(range(), i);
    matrix_op(eigenvec, lhs()); // calculate the local block of A*v
    rhs = lambda(i) * eigenvec; // and the local block of \lambda*v

    std::cout << i << ", block [" << block_start << ", "
              << (block_start + block_size - 1)
              << "]: deviation = " << norm2(rhs - lhs) / block_size
              << std::endl;
  }

  // Print some computation statistics
  if(comm_rank == 0) {
    auto stats = solver.stats();

    std::cout << "Number of Lanczos update iterations: " << stats.n_iter
              << std::endl;
    std::cout << "Total number of OP*x operations: " << stats.n_op_x_operations
              << std::endl;
    std::cout << "Total number of steps of re-orthogonalization: "
              << stats.n_reorth_steps << std::endl;
  }

  // Terminate MPI execution environment
  MPI_Finalize();

  return 0;
}
