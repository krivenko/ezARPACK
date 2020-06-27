.. _usage:

How to use ezARPACK in your project
===================================

Depending on what build system your C++ project is based on, there are
two ways to use ezARPACK. Either way, you would need

* A reasonably recent compiler supporting C++11 or newer;
* a working `ARPACK-NG installation
  <https://github.com/opencollab/arpack-ng>`_ (version 3.6.0 or newer).
* One of supported :ref:`linear algebra libraries <backends>` (unless
  you are going to :ref:`add support for a new one <new_backend>`).

Makefiles/no build system
-------------------------

Assume that `<EZARPACK_ROOT>` is either a directory with unpacked
ezARPACK sources or its :ref:`installation <installation>`
directory. A compiler command line for your program can be as simple as

.. code:: shell

  g++ -O3 -I<EZARPACK_ROOT>/include -L<ARPACK_NG_ROOT>/lib -larpack \
      -o myprog myprog.cpp

(similar for `clang++` and other compilers). More `-I` flags might be needed
if the linear algebra framework you choose is not visible to the compiler by
default.

.. note::

  When linking to the static version of ARPACK-NG library
  (`libarpack.a`), you might also need to explicitly link your code to
  a BLAS/LAPACK implementation. A common symptom of this problem is
  undefined references to `dgemv_` and other LAPACK subroutines.

CMake
-----

Here is a minimal example of a root `CMakeLists.txt` script for your
project.

.. code:: cmake

  cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

  project(myproject LANGUAGES CXX)

  # EZARPACK_ROOT is the installation prefix of ezARPACK
  set(ezARPACK_DIR ${EZARPACK_ROOT}/lib/cmake)

  # Import ezARPACK target
  find_package(ezARPACK 0.9 CONFIG REQUIRED)

  # Import Eigen (Blaze, Armadillo, etc) targets
  find_package(Eigen3 CONFIG)

  # Build an executable called 'myprog'
  add_executable(myprog myprog.cpp)
  target_link_libraries(myprog PRIVATE ezarpack Eigen3::Eigen)

If no usable ARPACK-NG lib has been detected by ezARPACK during
installation, you will have to link `myprog` to a ARPACK-NG library
explicitly.

.. code:: cmake

  # ARPACK_NG_ROOT is the installation prefix of ARPACK-NG
  set(arpack-ng_DIR ${ARPACK_NG_ROOT}/lib/cmake)

  # Import ARPACK-NG targets
  find_package(arpack-ng 3.6.0 REQUIRED)

  # Link to ARPACK-NG
  target_link_libraries(myprog ${arpack_ng_LIBRARIES})

