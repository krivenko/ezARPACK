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

#include <algorithm>
#include <iostream>

// This example shows how to use ezARPACK and the Armadillo storage backend
// to partially diagonalize a large sparse symmetric matrix
// and find a number of its low-lying eigenvalues.

#include <ezarpack/arpack_solver.hpp>
#include <ezarpack/storages/armadillo.hpp>
#include <ezarpack/version.hpp>

using namespace ezarpack;
using namespace arma;

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
  // For the Armadillo storage backend, other options would be
  // * `arpack_solver<Asymmetric, armadillo_storage>' for general real matrices;
  // * `arpack_solver<Complex, armadillo_storage>' for general complex matrices.
  arpack_solver<Symmetric, armadillo_storage> solver(N);

  using vector_view_t =
      arpack_solver<Symmetric, armadillo_storage>::vector_view_t;
  using vector_const_view_t =
      arpack_solver<Symmetric, armadillo_storage>::vector_const_view_t;

  // Linear operator representing multiplication of a given vector by our matrix
  // The operator must act on the 'in' vector and store results in 'out'.
  auto matrix_op = [](vector_const_view_t in, vector_view_t out) {
    out.zeros(); // Clear result

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

  // Specify parameters for the solver
  using params_t = arpack_solver<Symmetric, armadillo_storage>::params_t;
  params_t params(N_ev,               // Number of low-lying eigenvalues
                  params_t::Smallest, // We want the smallest eigenvalues
                  true);              // Yes, we want the eigenvectors
                                      // (Ritz vectors) as well

  // Run diagonalization!
  solver(matrix_op, params);

  // Number of converged eigenvalues
  std::cout << solver.nconv() << " out of " << params.n_eigenvalues
            << " eigenvalues have converged" << std::endl;

  // Print found eigenvalues
  std::cout << "Eigenvalues (Ritz values):" << std::endl;
  std::cout << solver.eigenvalues().t() << std::endl;

  // Check A*v = \lambda*v
  auto const& lambda = solver.eigenvalues();
  auto const& v = solver.eigenvectors();
  vec lhs(N), rhs(N);

  for(int i = 0; i < N_ev; ++i) {        // For each eigenpair ...
    matrix_op(v.col(i), lhs(span::all)); // calculate A*v
    rhs = lambda[i] * v.col(i);          // and \lambda*v

    std::cout << i << ": deviation = " << std::pow(norm(rhs - lhs), 2) / (N * N)
              << std::endl;
  }

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
