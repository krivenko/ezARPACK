/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2020 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/

#include <algorithm>
#include <iostream>

// This example shows how to use ezARPACK and the Eigen3 storage backend
// to partially diagonalize a large sparse symmetric matrix
// and find a number of its low-lying eigenvalues.

#include <ezarpack/arpack_worker.hpp>
#include <ezarpack/storages/eigen.hpp>
#include <ezarpack/version.hpp>

using namespace ezarpack;
using namespace Eigen;

// Size of the matrix
const int N = 10000;

// We are going to use a band matrix with this bandwidth
const int bandwidth = 5;

// The number of low-lying eigenvalues we want to compute
const int N_ev = 10;

// Linear operator representing multiplication of a given vector by our matrix.
// The operator must act on the 'in' vector and store results in 'out'.
//
// NB: With C++14 one could use a generic lambda function instead,
//
//  auto matrix_op = [](auto in, auto out) {
//                   ...
//  };
//
struct {
  template<typename TIn, typename TOut> void operator()(TIn in, TOut out) {
    out.fill(0); // Clear result

    // out_i = \sum_j A_{ij} in_j
    // A_{ij} = |i-j| / (1 + i + j), if |i-j| <= bandwidth, zero otherwise
    for(int i = 0; i < N; ++i) {
      int j_min = std::max(0, i - bandwidth);
      int j_max = std::min(N - 1, i + bandwidth);
      for(int j = j_min; j <= j_max; ++j) {
        out(i) += double(std::abs(i - j)) / (1 + i + j) * in(j);
      }
    }
  };
} matrix_op;

int main() {

  // Print ezARPACK version
  std::cout << "Using ezARPACK version " << EZARPACK_VERSION << std::endl;

  // Construct a worker object for the symmetric case.
  // For the Eigen3 storage backend, other options would be
  // * `arpack_worker<ezarpack::Asymmetric, eigen_storage>' for general
  //   real matrices;
  // * `arpack_worker<ezarpack::Complex, eigen_storage>' for general
  //   complex matrices.
  arpack_worker<ezarpack::Symmetric, eigen_storage> worker(N);

  // Specify parameters for the worker
  using params_t = arpack_worker<ezarpack::Symmetric, eigen_storage>::params_t;
  params_t params(N_ev,               // Number of low-lying eigenvalues
                  params_t::Smallest, // We want the smallest eigenvalues
                  true);              // Yes, we want the eigenvectors
                                      // (Ritz vectors) as well

  // Run diagonalization!
  worker(matrix_op, params);

  // Print found eigenvalues
  std::cout << "Eigenvalues (Ritz values):" << std::endl;
  std::cout << worker.eigenvalues().transpose() << std::endl;

  // Check A*v = \lambda*v
  auto const& lambda = worker.eigenvalues();
  auto const& v = worker.eigenvectors();
  VectorXd lhs(N), rhs(N);

  for(int i = 0; i < N_ev; ++i) {     // For each eigenpair ...
    matrix_op(v.col(i), lhs.head(N)); // calculate A*v
    rhs = lambda(i) * v.col(i);       // and \lambda*v

    std::cout << i << ": deviation = " << (rhs - lhs).squaredNorm() / (N * N)
              << std::endl;
  }

  // Print some computation statistics
  auto stats = worker.stats();

  std::cout << "Number of Arnoldi update iterations: " << stats.n_iter
            << std::endl;
  std::cout << "Number of 'converged' Ritz values: " << stats.n_converged
            << std::endl;
  std::cout << "Total number of OP*x operations: " << stats.n_op_x_operations
            << std::endl;
  std::cout << "Total number of steps of re-orthogonalization: "
            << stats.n_reorth_steps << std::endl;

  return 0;
}
