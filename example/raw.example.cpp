/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2019 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/

#include <iostream>
#include <algorithm>

// This example shows how to use ezARPACK and the raw memory storage backend
// to partially diagonalize a large sparse symmetric matrix
// and find a number of its low-lying eigenvalues.

#include <ezarpack/storages/raw.hpp>
#include <ezarpack/arpack_worker.hpp>
#include <ezarpack/version.hpp>

using namespace ezarpack;

// Size of the matrix
const int N = 10000;

// We are going to use a band matrix with this bandwidth
const int bandwidth = 5;

// The number of low-lying eigenvalues we want to compute
const int N_ev = 10;

int main(int argc, char* argv[]) {

  // Print ezARPACK version
  std::cout << "Using ezARPACK version " << EZARPACK_VERSION << std::endl;

  // Construct a worker object for the symmetric case.
  // For the raw memory storage backend, other options would be
  // * `arpack_worker<Asymmetric, raw_storage>' for general real matrices;
  // * `arpack_worker<Complex, raw_storage>' for general complex matrices.
  arpack_worker<Symmetric, raw_storage> worker(N);

  using vector_view_t = arpack_worker<Symmetric, raw_storage>::vector_view_t;
  using vector_const_view_t = arpack_worker<Symmetric, raw_storage>::vector_const_view_t;

  // Linear operator representing multiplication of a given vector by our matrix.
  // The operator must act on the 'from' vector and store results in 'to'.
  auto matrix_op = [](vector_const_view_t from, vector_view_t to) {
    std::fill(to, to + N, 0); // Clear result

    // to_i = \sum_j A_{ij} from_j
    // A_{ij} = |i-j| / (1 + i + j), if |i-j| <= bandwidth, zero otherwise
    for(int i = 0; i < N; ++i) {
      int j_min = std::max(0, i - bandwidth);
      int j_max = std::min(N, i + bandwidth);
      for(int j = j_min; j <= j_max; ++j) {
        to[i] += double(std::abs(i - j)) / (1 + i + j) * from[j];
      }
    }
  };

  // Specify parameters for the worker
  using params_t = arpack_worker<Symmetric, raw_storage>::params_t;
  params_t params(N_ev,               // Number of low-lying eigenvalues
                  params_t::Smallest, // We want the smallest eigenvalues
                  true                // Yes, we want the eigenvectors (Ritz vectors) as well
                  );

  // Run diagonalization!
  worker(matrix_op, params);

  // Print found eigenvalues
  auto const& lambda = worker.eigenvalues();
  std::cout << "Eigenvalues (Ritz values):\n[";
  for(int i = 0; i < N_ev - 1; ++i) {
    std::cout << lambda[i] << ",";
  }
  std::cout << lambda[N_ev - 1] << "]" << std::endl;

  // Check A*v = \lambda*v
  // NB: Eigenvectors are stored in the column major order
  auto const& v = worker.eigenvectors();

  double * lhs = new double[N];
  double * rhs = new double[N];

  for(int i = 0; i < N_ev; ++i) {                      // For each eigenpair ...
    matrix_op(v + N*i, lhs);                            // calculate A*v
    std::transform(v + N*i, v + N*(i+1), rhs,           // and \lambda*v
                  [&](double x) { return lambda[i]*x;}
                  );

    double deviation = 0;
    for(int j = 0; j < N; ++j) {
      double d = rhs[j] - lhs[j];
      deviation += d * d;
    }
    deviation /= N*N;
    std::cout << i << ": deviation = " << deviation << std::endl;
  }

  delete[] lhs;
  delete[] rhs;

  // Print some computation statistics
  auto stats = worker.stats();

  std::cout << "Number of Arnoldi update iterations: " << stats.n_iter << std::endl;
  std::cout << "Number of 'converged' Ritz values: " << stats.n_converged << std::endl;
  std::cout << "Total number of OP*x operations: " << stats.n_op_x_operations << std::endl;
  std::cout << "Total number of steps of re-orthogonalization: " << stats.n_reorth_steps << std::endl;

  return 0;
}
