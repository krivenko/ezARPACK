ezARPACK
========

[![Build Status](https://travis-ci.org/krivenko/ezARPACK.svg?branch=master)](
https://travis-ci.org/krivenko/ezARPACK)
[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-red)](
https://krivenko.github.io/ezARPACK)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3930203.svg)](
https://doi.org/10.5281/zenodo.3930203)

ezARPACK is a C++11 wrapper around ARPACK-NG [1] that can be used in conjunction
with a number of C++ vector/matrix algebra libraries. It allows for solving
large scale eigenproblems for real symmetric, non-symmetric and complex matrices
with a minimal amount of boilerplate code.

When used directly, ARPACK-NG does not force the user to stick to a predefined
storage format of the matrix being diagonalized. Instead, on each iteration of
the Arnoldi/Lanczos algorithm the user code is expected to apply the
corresponding linear operator to the vector (memory buffer) passed to it and
store the result in another buffer. ezARPACK retains this idea allowing to use
any callable C++ object as the linear operator.

Another important feature of ezARPACK is its extensibility with respect to
compatible matrix algebra libraries. Currently, it supports the following
libraries (storage backends):

* Eigen 3 [2];
* Blaze >= 3 [3];
* Armadillo [4];
* Boost uBLAS >= 1.58 [5];
* TRIQS arrays >= 2.0 [6];
* xtensor >= 0.20.0 [7];
* Raw memory buffers *(for unit testing, not recommended for general use)*.

One can easily add support for their favorite vector/matrix framework by
defining a new specialization of the `storage_traits` structure.

Copyright (C) 2016-2020 Igor Krivenko <igor.s.krivenko @ gmail.com>

---
**ATTENTION**

Version 0.9 comes with a complete set of documentation and a few breaking
changes.

* A more common term 'solver' is now used instead of 'worker' everywhere in the
  code. In particular, `arpack_worker` has been renamed to `arpack_solver`.
* Method `arpack_solver::from_vector_n()` has been renamed to
  `arpack_solver::in_vector_n()`.
* Method `arpack_solver::to_vector_n()` has been renamed to
  `arpack_solver::out_vector_n()`.
* Computational mode `Invert` has been renamed to `Inverse`.
* Field `n_converged` has been removed from the `stats_t` structures.
  The number of converged Ritz values is now returned by method
  `arpack_solver::nconv()`.
---

Dependencies
------------

ezARPACK is a header-only library that has no external dependencies.

However, one will need a working installation of ARPACK-NG 3.6.0 or newer [1]
in order to compile examples and unit tests. Futhermore, specific examples and
tests will only be built if the respective matrix algebra library is detected by
CMake (does not apply to the raw memory storage backend).

Installation
------------

ezARPACK is usable without installation, just add
`-I/<path_to_ezARPACK_sources>/include` to the compiler command line and
`-L/<ARPACK-NG_installation_prefix>/lib -larpack` to the linker command line.

You will need CMake version 3.1.0 or newer [8] to build examples/unit tests and
to install ezARPACK so that it can be used from other CMake projects.

Assuming that ezARPACK is to be installed in `<ezARPACK_installation_prefix>`,
the installation normally proceeds in a few simple steps.

```
$ git clone https://github.com/krivenko/ezARPACK.git ezARPACK.git
$ mkdir ezARPACK.build && cd ezARPACK.build
$ cmake ../ezARPACK.git                                   \
$ -DCMAKE_INSTALL_PREFIX=<ezARPACK_installation_prefix>   \
  -DARPACK_NG_ROOT=<ARPACK-NG_installation_prefix>        \
  -DEigen3_ROOT=<Eigen3_installation_prefix>              \
  -Dblaze_ROOT=<Blaze_installation_prefix>                \
  -DArmadillo_ROOT=<Armadillo_installation_prefix>        \
  -DBOOST_ROOT=<Boost_installation_prefix>                \
  -DTRIQS_ROOT=<TRIQS_installation_prefix>                \
  -Dxtensor_ROOT=<xtensor_installation_prefix>            \
  -Dxtensor-blas_ROOT=<xtensor-blas_installation_prefix>  \
  -DExamples=ON                                           \
  -DTests=ON
$ make
$ make test
$ make install
```

Compilation of the tests can be disabled with CMake flag `-DTests=OFF`
*(not recommended)*.

Examples are compiled by default, disable them with `-DExamples=OFF`.

CMake options specific to individual storage backends (`Eigen3_ROOT`,
`blaze_ROOT`, `Armadillo_ROOT`, `BOOST_ROOT`, `TRIQS_ROOT`,
`xtensor_ROOT`/`xtensor-blas_ROOT`) can be omitted if the respective libraries
are installed in the standard system locations. If some of the libraries are not
found, CMake will skip the corresponding examples and unit tests.

Usage
-----

Once ezARPACK is installed, you can use it in your CMake project. Here is
a minimal example of an application `CMakeLists.txt` file.

```cmake
cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

project(myproject LANGUAGES CXX)

# ARPACK_NG_ROOT and EZARPACK_ROOT are installation prefixes of
# ARPACK-NG and ezARPACK respectively
set(arpack-ng_DIR ${ARPACK_NG_ROOT}/lib/cmake)
set(ezARPACK_DIR ${EZARPACK_ROOT}/lib/cmake)

# Import ARPACK-NG targets
find_package(arpack-ng 3.6.0 REQUIRED)
# Import ezARPACK targets
find_package(ezARPACK 0.9 CONFIG REQUIRED)
# Import Eigen (Blaze, Armadillo, etc) targets
find_package(Eigen3 CONFIG)

# Build an executable called 'test'
add_executable(test test.cpp)

# Make ezARPACK and Eigen headers visible to the compiler
# and link to ARPACK-NG libraries.
target_link_libraries(test PRIVATE
                      ezarpack Eigen3::Eigen ${arpack_ng_LIBRARIES})
```

Here is how `test.cpp` could look like.
```c++
#include <iostream>

// This example shows how to use ezARPACK and the Eigen3 storage backend
// to partially diagonalize a large sparse symmetric matrix
// and find a number of its low-lying eigenvalues.

#include <ezarpack/storages/eigen.hpp>
#include <ezarpack/arpack_solver.hpp>
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

  // Construct a solver object for the symmetric case.
  // For the Eigen3 storage backend, other options would be
  // * `arpack_solver<ezarpack::Asymmetric, eigen_storage>' for general
  //   real matrices;
  // * `arpack_solver<ezarpack::Complex, eigen_storage>' for general
  //   complex matrices.
  arpack_solver<ezarpack::Symmetric, eigen_storage> solver(N);

  // Specify parameters for the solver
  using params_t = arpack_solver<ezarpack::Symmetric, eigen_storage>::params_t;
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
  std::cout << solver.eigenvalues().transpose() << std::endl;

  // Check A*v = \lambda*v
  auto const& lambda = solver.eigenvalues();
  auto const& v = solver.eigenvectors();
  VectorXd lhs(N), rhs(N);

  for(int i = 0; i < N_ev; ++i) {     // For each eigenpair ...
    matrix_op(v.col(i), lhs.head(N)); // calculate A*v
    rhs = lambda(i) * v.col(i);       // and \lambda*v

    std::cout << i << ": deviation = "
              << (rhs - lhs).squaredNorm() / (N*N) << std::endl;
  }

  // Print some computation statistics
  auto stats = solver.stats();

  std::cout << "Number of Arnoldi update iterations: "
            << stats.n_iter << std::endl;
  std::cout << "Total number of OP*x operations: "
            << stats.n_op_x_operations << std::endl;
  std::cout << "Total number of steps of re-orthogonalization: "
            << stats.n_reorth_steps << std::endl;

  return 0;
}
```

Citing
------

If you find this library useful for your research, you can help me by citing it
using the following BibTeX entry.

```
@software{igor_krivenko_2020_3930203,
  author       = {Igor Krivenko},
  title        = {{ezARPACK - a C++ ARPACK-NG wrapper compatible with
                   multiple matrix/vector algebra libraries: Release
                   0.9}},
  month        = jul,
  year         = 2020,
  publisher    = {Zenodo},
  version      = {0.9},
  doi          = {10.5281/zenodo.3930203},
  url          = {https://doi.org/10.5281/zenodo.3930203}
}
```

Known issues
------------

* ezARPACK is still in beta, use with caution!
* Parallel ARPACK routines (PARPACK) are not supported.

License
-------

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

[1]: https://github.com/opencollab/arpack-ng
[2]: http://eigen.tuxfamily.org
[3]: https://bitbucket.org/blaze-lib/blaze
[4]: http://arma.sourceforge.net
[5]: https://www.boost.org/doc/libs/1_58_0/libs/numeric/ublas/doc
[6]: https://triqs.github.io/triqs/latest
[7]: https://github.com/xtensor-stack/xtensor
[8]: https://cmake.org/download
