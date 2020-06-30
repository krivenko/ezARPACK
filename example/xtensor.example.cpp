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

// This example shows how to use ezARPACK and the xtensor storage backend
// to partially diagonalize a large sparse symmetric matrix
// and find a number of its low-lying eigenvalues.

#include <ezarpack/arpack_worker.hpp>
#include <ezarpack/storages/xtensor.hpp>
#include <ezarpack/version.hpp>

#include <xtensor/xio.hpp>
#include <xtensor/xnorm.hpp>

using namespace ezarpack;
using namespace xt;

// Size of the matrix
const int N = 10000;

// We are going to use a band matrix with this bandwidth
const int bandwidth = 5;

// The number of low-lying eigenvalues we want to compute
const int N_ev = 10;

int main() {

  // Print ezARPACK version
  std::cout << "Using ezARPACK version " << EZARPACK_VERSION << std::endl;

  // Construct a worker object for the symmetric case.
  // For the xtensor storage backend, other options would be
  // * `arpack_worker<ezarpack::Asymmetric, xtensor_storage>' for general
  //   real matrices;
  // * `arpack_worker<ezarpack::Complex, xtensor_storage>' for general
  //   complex matrices.
  arpack_worker<ezarpack::Symmetric, xtensor_storage> worker(N);

  // Linear operator representing multiplication of a given vector by our matrix
  // The operator must act on the 'in' vector and store results in 'out'.
  auto matrix_op = [](auto in, auto out) {
    out.fill(.0); // Clear result

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

  // Specify parameters for the worker
  using params_t =
      arpack_worker<ezarpack::Symmetric, xtensor_storage>::params_t;
  params_t params(N_ev,               // Number of low-lying eigenvalues
                  params_t::Smallest, // We want the smallest eigenvalues
                  true);              // Yes, we want the eigenvectors
                                      // (Ritz vectors) as well

  // Run diagonalization!
  worker(matrix_op, params);

  // Print found eigenvalues
  std::cout << "Eigenvalues (Ritz values):" << std::endl;
  std::cout << worker.eigenvalues() << std::endl;

  // Check A*v = \lambda*v
  auto const& lambda = worker.eigenvalues();
  auto const& v = worker.eigenvectors();

  auto lhs = xt::xtensor<double, 1>::from_shape({N});
  auto rhs = xt::xtensor<double, 1>::from_shape({N});

  for(int i = 0; i < N_ev; ++i) { // For each eigenpair ...
    auto col = xt::view(v, xt::all(), i);
    matrix_op(col, xt::view(lhs, xt::all())); // calculate A*v
    rhs = lambda(i) * col;                    // and \lambda*v

    std::cout << i << ": deviation = " << xt::norm_sq(rhs - lhs) / (N * N)
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
