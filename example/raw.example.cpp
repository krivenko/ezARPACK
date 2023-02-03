/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2023 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/

#include <algorithm>
#include <cmath>
#include <iostream>

// This example shows how to use ezARPACK and the raw memory storage backend
// to partially diagonalize a large sparse symmetric matrix
// and find a number of its low-lying eigenvalues.

#include <ezarpack/arpack_solver.hpp>
#include <ezarpack/storages/raw.hpp>
#include <ezarpack/version.hpp>

using namespace ezarpack;

// Size of the matrix
const int N = 10000;

// We are going to use a band matrix with this bandwidth
const int bandwidth = 5;

// The number of low-lying eigenvalues we want to compute
const int N_ev = 10;

int main() {

  // Print ezARPACK version
  std::cout << "Using ezARPACK version " << EZARPACK_VERSION << std::endl;

  // Construct a solver object for the symmetric case.
  // For the raw memory storage backend, other options would be
  // * `arpack_solver<Asymmetric, raw_storage>' for general real matrices;
  // * `arpack_solver<Complex, raw_storage>' for general complex matrices.
  using solver_t = arpack_solver<Symmetric, raw_storage>;
  solver_t solver(N);

  // Specify parameters for the solver
  using params_t = solver_t::params_t;
  params_t params(N_ev,               // Number of low-lying eigenvalues
                  params_t::Smallest, // We want the smallest eigenvalues
                  true);              // Yes, we want the eigenvectors
                                      // (Ritz vectors) as well

  // Linear operator representing multiplication of a given vector by our matrix
  // The operator must act on the 'in' vector and store results in 'out'.
  auto matrix_op = [](solver_t::vector_const_view_t in,
                      solver_t::vector_view_t out) {
    std::fill(out, out + N, 0); // Clear result

    // out_i = \sum_j A_{ij} in_j
    // A_{ij} = |i-j| / (1 + i + j), if |i-j| <= bandwidth, zero otherwise
    for(int i = 0; i < N; ++i) {
      int j_min = std::max(0, i - bandwidth);
      int j_max = std::min(N - 1, i + bandwidth);
      for(int j = j_min; j <= j_max; ++j) {
        out[i] += double(std::abs(i - j)) / (1 + i + j) * in[j];
      }
    }
  };

  // Run diagonalization!
  solver(matrix_op, params);

  // Number of converged eigenvalues
  std::cout << solver.nconv() << " out of " << params.n_eigenvalues
            << " eigenvalues have converged" << std::endl;

  // Print found eigenvalues
  auto const& lambda = solver.eigenvalues();
  std::cout << "Eigenvalues (Ritz values):\n[";
  for(int i = 0; i < N_ev - 1; ++i) {
    std::cout << lambda[i] << ",";
  }
  std::cout << lambda[N_ev - 1] << "]" << std::endl;

  // Check A*v = \lambda*v
  // NB: Eigenvectors are stored in the column-major order
  auto const& v = solver.eigenvectors();

  double* lhs = new double[N];
  double* rhs = new double[N];

  for(int i = 0; i < N_ev; ++i) { // For each eigenpair ...
    auto const eigenvec = v + N * i;
    matrix_op(eigenvec, lhs);                   // calculate A*v
    std::transform(eigenvec, eigenvec + N, rhs, // and \lambda*v
                   [&](double x) { return lambda[i] * x; });

    double deviation = 0;
    for(int j = 0; j < N; ++j) {
      double d = rhs[j] - lhs[j];
      deviation += d * d;
    }
    deviation = std::sqrt(deviation) / N;
    std::cout << i << ": deviation = " << deviation << std::endl;
  }

  delete[] lhs;
  delete[] rhs;

  // Print some computation statistics
  auto stats = solver.stats();

  std::cout << "Number of Lanczos update iterations: " << stats.n_iter
            << std::endl;
  std::cout << "Total number of OP*x operations: " << stats.n_op_x_operations
            << std::endl;
  std::cout << "Total number of steps of re-orthogonalization: "
            << stats.n_reorth_steps << std::endl;

  return 0;
}
