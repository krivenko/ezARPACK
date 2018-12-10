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
* TRIQS arrays [2].

One can easily add support for her favorite vector/matrix framework by defining
a new instance of the `storage_traits` structure (see, for example, `include/storages/triqs.hpp`).

Copyright (C) 2016-2018 by I. Krivenko

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

You will need CMake version 3.0.2 or newer [3] to build examples/unit tests and to install ezARPACK
such that it can be used from other CMake projects.

Assuming that ezARPACK is to be installed in `<ezARPACK_installation_prefix>`, the installation
normally proceeds in a few simple steps.

```
$ git clone https://github.com/krivenko/ezARPACK.git ezARPACK.git
$ mkdir ezARPACK.build && cd ezARPACK.build
$ cmake ../ezARPACK.git                                 \
$ -DCMAKE_INSTALL_PREFIX=<ezARPACK_installation_prefix> \
  -DARPACK_NG_ROOT=<ARPACK-NG_installation_prefix>      \
  -DTRIQS_ROOT=<TRIQS_installation_prefix>              \
  -DExamples=ON                                         \
  -DTests=ON
$ make
$ make test
$ make install
```

Compilation of the tests can be disabled with CMake flag `-DTests=OFF` *(not recommended)*.

Examples are compiled by default, disable them with `-DExamples=OFF`.

CMake options specific to individual storage backends (`TRIQS_ROOT`) can be omitted if
the respective libraries are installed in the standard system locations. If some of the
libraries are not found, CMake will skip the corresponding examples and unit tests.

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
find_package(ezARPACK 0.2 CONFIG REQUIRED)

# Build an executable called 'test'
add_executable(test test.cpp)

# Make ezARPACK headers visible to the compiler
# and link to ARPACK-NG libraries
target_link_libraries(test ezarpack arpack)
```

Documentation
-------------

For now, I only provide a few examples in the `example` directory.

License
-------

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

[1]: https://github.com/opencollab/arpack-ng
[2]: https://triqs.github.io/triqs/master
[3]: https://cmake.org/download
