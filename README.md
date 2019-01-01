ezARPACK
========

ezARPACK is a C++11 wrapper around ARPACK-NG [1] that can be used in conjunction with
a number of C++ vector/matrix algebra libraries. It allows for solving large scale
eigenproblems for symmetric, asymmetric and complex matrices with a minimal amount
of boilerplate code.

When used directly, ARPACK-NG does not force the user to stick to a predefined storage
format of the matrix being diagonalized. Instead, on each iteration of the
Arnoldi/Lanczos algorithm the user code is expected to apply the corresponding linear
operator to the vector (memory buffer) passed to it and store the result in another buffer.
ezARPACK retains this idea allowing to use any callable C++ object as the linear operator.

Another important feature of ezARPACK is its extensibility with respect to compatible
matrix algebra libraries. Currently, it supports the following libraries (storage backends):

* Raw memory buffers *(not recommended for general use)*;
* Eigen 3 [2];
* Blaze >= 3 [3];
* Armadillo [4];
* Boost uBLAS >= 1.58 [5];
* TRIQS arrays >= 2.0 [6].

One can easily add support for her favorite vector/matrix framework by defining
a new instance of the `storage_traits` structure (see, for example, `include/storages/eigen.hpp`).

Copyright (C) 2016-2019 by I. Krivenko

Dependencies
------------

ezARPACK is a header-only library that has no external dependencies.

However, one will need a working installation of ARPACK-NG 3.6.0 or newer [1] in order to compile
examples and unit tests. Futhermore, specific examples and tests will only be built if the
respective matrix algebra library is detected by CMake (does not apply to the raw memory storage
backend).

Installation
------------

ezARPACK is usable without installation, just add `-I/<path_to_ezARPACK_sources>/include`
to the compiler command line and `-L/<ARPACK-NG_installation_prefix>/lib -larpack` to
the linker command line.

You will need CMake version 3.0.2 or newer [7] to build examples/unit tests and to install ezARPACK
such that it can be used from other CMake projects.

Assuming that ezARPACK is to be installed in `<ezARPACK_installation_prefix>`, the installation
normally proceeds in a few simple steps.

```
$ git clone https://github.com/krivenko/ezARPACK.git ezARPACK.git
$ mkdir ezARPACK.build && cd ezARPACK.build
$ cmake ../ezARPACK.git                                 \
$ -DCMAKE_INSTALL_PREFIX=<ezARPACK_installation_prefix> \
  -DARPACK_NG_ROOT=<ARPACK-NG_installation_prefix>      \
  -DEigen3_ROOT=<Eigen3_installation_prefix>            \
  -Dblaze_ROOT=<Blaze_installation_prefix>              \
  -DArmadillo_ROOT=<Armadillo_installation_prefix>      \
  -DBOOST_ROOT=<Boost_installation_prefix>              \
  -DTRIQS_ROOT=<TRIQS_installation_prefix>              \
  -DExamples=ON                                         \
  -DTests=ON
$ make
$ make test
$ make install
```

Compilation of the tests can be disabled with CMake flag `-DTests=OFF` *(not recommended)*.

Examples are compiled by default, disable them with `-DExamples=OFF`.

CMake options specific to individual storage backends (`Eigen3_ROOT`, `blaze_ROOT`, `Armadillo_ROOT`,
`BOOST_ROOT`, `TRIQS_ROOT`) can be omitted if the respective libraries are installed in the
standard system locations. If some of the libraries are not found, CMake will skip the corresponding
examples and unit tests.

Usage
-----

Once ezARPACK is installed, you can use it in your CMake project. Here is a minimal
example of an application `CMakeLists.txt` file.

```cmake
cmake_minimum_required(VERSION 3.0.2 FATAL_ERROR)

project(myproject LANGUAGES CXX)

# ARPACK_NG_ROOT and EZARPACK_ROOT are installation prefixes of
# ARPACK-NG and ezARPACK respectively
set(arpack-ng_DIR ${ARPACK_NG_ROOT}/lib/cmake)
set(ezARPACK_DIR ${EZARPACK_ROOT}/lib/cmake)

# Import ARPACK-NG targets
find_package(arpack-ng 3.6.0 REQUIRED)
# Import ezARPACK targets
find_package(ezARPACK 0.6 CONFIG REQUIRED)

# Build an executable called 'test'
add_executable(test test.cpp)

# Make ezARPACK headers visible to the compiler
# and link to ARPACK-NG libraries
target_link_libraries(test ezarpack arpack)
```

Here is how `test.cpp` could look like.
```c++
#include <iostream>

// This example shows how to use ezARPACK and the Eigen3 storage backend
// to partially diagonalize a large sparse symmetric matrix
// and find a number of its low-lying eigenvalues.

#include <ezarpack/storages/eigen.hpp>
#include <ezarpack/arpack_worker.hpp>
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
// The operator must act on the 'from' vector and store results in 'to'.
//
// NB: With C++14 one could use a generic lambda function instead,
//
//  auto matrix_op = [](auto from, auto to) {
//                   ...
//  };
//
struct {
  template<typename TFrom, typename TTo> void operator()(TFrom from, TTo to) {
    to.fill(0); // Clear result

    // to_i = \sum_j A_{ij} from_j
    // A_{ij} = |i-j| / (1 + i + j), if |i-j| <= bandwidth, zero otherwise
    for(int i = 0; i < N; ++i) {
      int j_min = std::max(0, i - bandwidth);
      int j_max = std::min(N, i + bandwidth);
      for(int j = j_min; j <= j_max; ++j) {
        to(i) += double(std::abs(i - j)) / (1 + i + j) * from(j);
      }
    }
  };
} matrix_op;

int main(int argc, char* argv[]) {

  // Print ezARPACK version
  std::cout << "Using ezARPACK version " << EZARPACK_VERSION << std::endl;

  // Construct a worker object for the symmetric case.
  // For the Eigen3 storage backend, other options would be
  // * `arpack_worker<ezarpack::Asymmetric, eigen_storage>' for general real matrices;
  // * `arpack_worker<ezarpack::Complex, eigen_storage>' for general complex matrices.
  arpack_worker<ezarpack::Symmetric, eigen_storage> worker(N);

  // Specify parameters for the worker
  using params_t = arpack_worker<ezarpack::Symmetric, eigen_storage>::params_t;
  params_t params(N_ev,               // Number of low-lying eigenvalues
                  params_t::Smallest, // We want the smallest eigenvalues
                  true                // Yes, we want the eigenvectors (Ritz vectors) as well
                  );

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

    std::cout << i << ": deviation = " << (rhs - lhs).squaredNorm() / (N*N) << std::endl;
  }

  // Print some computation statistics
  auto stats = worker.stats();

  std::cout << "Number of Arnoldi update iterations: " << stats.n_iter << std::endl;
  std::cout << "Number of 'converged' Ritz values: " << stats.n_converged << std::endl;
  std::cout << "Total number of OP*x operations: " << stats.n_op_x_operations << std::endl;
  std::cout << "Total number of steps of re-orthogonalization: " << stats.n_reorth_steps << std::endl;

  return 0;
}
```

Documentation
-------------

For now, I only provide a few examples in the `example` directory.

Known issues
------------

* ezARPACK is still beta, use with caution!
* Parallel ARPACK routines (PARPACK) are not supported.
* `arpack_worker<Asymmetric, ...>` will refuse to run in `ShiftAndInvertReal` and
  `ShiftAndInvertImag` modes (`dnaupd` modes 3 and 4). This is a temporary
  workaround for [a `dneupd` issue](http://forge.scilab.org/index.php/p/arpack-ng/issues/1315/).

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
[6]: https://triqs.github.io/triqs/master
[7]: https://cmake.org/download
