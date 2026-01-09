ezARPACK
========

[![Build status](https://github.com/krivenko/ezARPACK/actions/workflows/CI.yml/badge.svg)](
https://github.com/krivenko/ezARPACK/actions/workflows/CI.yml)
[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-red)](
https://krivenko.github.io/ezARPACK)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7047515.svg)](
https://doi.org/10.5281/zenodo.7047515)

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
* TRIQS arrays (part of the TRIQS libraries 2.0--3.0.x) [6];
* TRIQS/nda library >= 1.1.0 [7];
* xtensor >= 0.20.0 [8];
* Raw memory buffers *(for unit testing, not recommended for general use)*.

One can easily add support for their favorite vector/matrix framework by
defining a new specialization of the `storage_traits` structure.

Since version 1.0, ezARPACK wraps Parallel ARPACK (PARPACK) routines with MPI
message passing layer in addition to the serial routines.

Copyright (C) 2016-2026 Igor Krivenko <iskrivenko [at] proton [dot] me>

Dependencies
------------

ezARPACK is a header-only library that has no external dependencies.

However, one will need a working installation of ARPACK-NG 3.6.0 or newer [1]
in order to compile examples and unit tests. Furthermore, specific examples and
tests will only be built if the respective matrix algebra library is detected by
CMake (does not apply to the raw memory storage backend).

Installation
------------

ezARPACK is usable without installation, just add
`-I/<path_to_ezARPACK_sources>/include` to the compiler command line and
`-L/<ARPACK-NG_installation_prefix>/lib -larpack` to the linker command line.

You will need CMake version 3.18.0 or newer [9] to build examples/unit tests and
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
  -Dnda_ROOT=<nda_installation_prefix>                    \
  -Dxtensor_ROOT=<xtensor_installation_prefix>            \
  -Dxtensor-blas_ROOT=<xtensor-blas_installation_prefix>  \
  -DEnableMPI=ON                                          \
  -DExamples=ON                                           \
  -DTests=ON
$ make
$ make test
$ make install
```

Compilation of the tests can be disabled with CMake flag `-DTests=OFF`
*(not recommended)*.

Examples are compiled by default, disable them with `-DExamples=OFF`.

Detection of an MPI implementation and compilation of the MPI-enabled unit tests
and examples can be skipped by setting `-DEnableMPI=OFF`.

CMake options specific to individual storage backends (`Eigen3_ROOT`,
`blaze_ROOT`, `Armadillo_ROOT`, `BOOST_ROOT`, `TRIQS_ROOT`, `nda_ROOT`
`xtensor_ROOT`/`xtensor-blas_ROOT`) can be omitted if the respective libraries
are installed in the standard system locations. If some of the libraries are not
found, CMake will skip the corresponding examples and unit tests.

Another way to install the most recent released version of ezARPACK is via
[`conda`](https://www.anaconda.com/) package manager.

```
$ conda install -c krivenko ezarpack
```

Usage
-----

Once ezARPACK is installed, you can use it in your CMake project. Here is
a minimal example of an application `CMakeLists.txt` file.

```cmake
cmake_minimum_required(VERSION 3.18.0 FATAL_ERROR)

project(myproject LANGUAGES CXX)

# ezARPACK_ROOT is the installation prefix of ezARPACK.
set(ezARPACK_DIR ${ezARPACK_ROOT}/lib/cmake)

# Import ezARPACK target.
find_package(ezARPACK 1.0 CONFIG REQUIRED)

# Import Eigen (Blaze, Armadillo, etc) targets.
find_package(Eigen3 CONFIG REQUIRED)

# Build an executable called `myprog`.
add_executable(myprog myprog.cpp)
target_link_libraries(myprog ezarpack Eigen3::Eigen)

# Find a usable version of ARPACK-NG.
# Macro find_arpackng() can be instructed to use a specific ARPACK-NG
# installation by setting the CMake variable `ARPACK_NG_ROOT`.
find_arpackng(3.6.0 REQUIRED)

# Link the executable to the ARPACK library.
target_link_libraries(myprog ${ARPACK_LIBRARIES})
```

Here is how `myprog.cpp` could look like.
```c++
#include <cmath>
#include <iostream>

// This example shows how to use ezARPACK and the Eigen3 storage backend
// to partially diagonalize a large sparse symmetric matrix
// and find a number of its low-lying eigenvalues.

#include <ezarpack/arpack_solver.hpp>
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

int main() {

  // Print ezARPACK version
  std::cout << "Using ezARPACK version " << EZARPACK_VERSION << std::endl;

  // Construct a solver object for the symmetric case.
  // For the Eigen3 storage backend, other options would be
  // * `arpack_solver<ezarpack::Asymmetric, eigen_storage>' for general
  //   real matrices;
  // * `arpack_solver<ezarpack::Complex, eigen_storage>' for general
  //   complex matrices.
  using solver_t = arpack_solver<ezarpack::Symmetric, eigen_storage>;
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

  for(int i = 0; i < N_ev; ++i) { // For each eigenpair ...
    const VectorXd eigenvec = v.col(i);
    matrix_op(eigenvec.head(N), lhs.head(N)); // calculate A*v
    rhs = lambda(i) * eigenvec;               // and \lambda*v

    std::cout << i << ": deviation = " << (rhs - lhs).norm() / N
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
```

The same eigenproblem can be solved using an MPI-parallelized solver that wraps
PARPACK routines. In this case one needs to additionally link the executable
to MPI libraries.

```cmake
# Parallel ARPACK (MPI)

# Build another executable `myprog_mpi`.
add_executable(myprog_mpi myprog_mpi.cpp)
target_link_libraries(myprog_mpi ezarpack Eigen3::Eigen)

# Detect an MPI-3.0 implementation.
find_package(MPI 3.0 REQUIRED)

# Link the executable to the Parallel ARPACK library and to the MPI.
target_include_directories(myprog_mpi PRIVATE ${MPI_CXX_INCLUDE_PATH})
target_link_libraries(myprog_mpi ${PARPACK_LIBRARIES} ${MPI_CXX_LIBRARIES})
```

```c++
#include <cmath>
#include <iostream>
#include <vector>

#include <ezarpack/mpi/arpack_solver.hpp>
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

int main(int argc, char* argv[]) {

  // Initialize MPI environment
  MPI_Init(&argc, &argv);

  // Call utility functions from namespace 'ezarpack::mpi' to find out
  // the world communicator size and the rank of the calling process.
  const int comm_size = mpi::size(MPI_COMM_WORLD);
  const int comm_rank = mpi::rank(MPI_COMM_WORLD);

  // Print ezARPACK version
  if(comm_rank == 0)
    std::cout << "Using ezARPACK version " << EZARPACK_VERSION << std::endl;

  // Construct an MPI-parallelized solver object for the symmetric case.
  // For the Eigen3 storage backend, other options would be
  // * `mpi::arpack_solver<ezarpack::Asymmetric, eigen_storage>' for general
  //   real matrices;
  // * `mpi::arpack_solver<ezarpack::Complex, eigen_storage>' for general
  //   complex matrices.
  using solver_t = mpi::arpack_solver<ezarpack::Symmetric, eigen_storage>;
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
    block_sizes[rank] = mpi::compute_local_block_size(N, comm_size, rank);

  // Temporary vector used in distributed matrix-vector multiplication
  VectorXd local_op_in = VectorXd(N);

  // Linear operator representing multiplication of a given vector by our matrix
  auto matrix_op = [&](solver_t::vector_const_view_t in,
                       solver_t::vector_view_t out) {
    // 'in' and 'out' are views of the locally stored blocks of their respective
    // distributed N-dimensional vectors. Therefore, matrix-vector
    // multiplication has to be performed in two steps.

    // 1. Local multiplication of A's columns
    // [block_start; block_start + block_size] by 'in'. The result is an
    // N-dimensional vector stored in 'local_op_in'.
    local_op_in.fill(0);
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
    MPI_Reduce_scatter(local_op_in.data(), out.data(), block_sizes.data(),
                       MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  };

  // Run diagonalization!
  solver(matrix_op, params);

  if(comm_rank == 0) {
    // Number of converged eigenvalues
    std::cout << solver.nconv() << " out of " << params.n_eigenvalues
              << " eigenvalues have converged" << std::endl;

    // Print found eigenvalues
    std::cout << "Eigenvalues (Ritz values):" << std::endl;
    std::cout << solver.eigenvalues().transpose() << std::endl;
  }

  // Check A*v = \lambda*v
  auto const& lambda = solver.eigenvalues();
  auto const& v = solver.eigenvectors();
  VectorXd lhs(block_size), rhs(block_size);

  for(int i = 0; i < N_ev; ++i) { // For each eigenpair ...
    const VectorXd eigenvec_block = v.col(i);
    matrix_op(eigenvec_block.head(block_size),
              lhs.head(block_size));  // calculate the local block of A*v
    rhs = lambda(i) * eigenvec_block; // and the local block of \lambda*v

    std::cout << i << ", block [" << block_start << ", "
              << (block_start + block_size - 1)
              << "]: deviation = " << (rhs - lhs).norm() / block_size
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
```

Citing
------

If you find this library useful for your research, you can help me by citing it
using the following BibTeX entry.

```
@software{igor_krivenko_2022_7047515,
  author       = {Igor Krivenko},
  title        = {{ezARPACK - a C++ ARPACK-NG wrapper compatible with
                   multiple matrix/vector algebra libraries: Release
                   1.0}},
  month        = sep,
  year         = 2022,
  publisher    = {Zenodo},
  version      = {1.0},
  doi          = {10.5281/zenodo.7047515},
  url          = {https://doi.org/10.5281/zenodo.7047515}
}
```

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
[7]: https://github.com/triqs/nda
[8]: https://github.com/xtensor-stack/xtensor
[9]: https://cmake.org/download
