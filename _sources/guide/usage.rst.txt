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
* A working implementation of MPI-3.0 or newer, such as
  `OpenMPI <https://www.open-mpi.org/>`_ or `MPICH <https://www.mpich.org/>`_,
  if you plan to use the :ref:`Parallel ARPACK solvers <mpi>`.

.. warning::
  The Parallel ARPACK eigensolvers for general complex matrices are
  `known to be unstable <https://github.com/opencollab/arpack-ng/pull/245>`_
  and can yield wrong results in ARPACK-NG versions prior to 3.8.0.
  Be sure to link your code to ARPACK-NG 3.8.0 or newer if you use the
  :ref:`complex <complex>` variant of :cpp:type:`ezarpack::mpi::arpack_solver`.

Makefiles/no build system
-------------------------

Assume that ``<EZARPACK_ROOT>`` is either a directory with unpacked
ezARPACK sources or its :ref:`installation <installation>`
directory. A compiler command line for your program can be as simple as

.. code:: shell

  g++ -O3 -I<EZARPACK_ROOT>/include -L<ARPACK_NG_ROOT>/lib -larpack \
      -o myprog myprog.cpp

(similar for ``clang++`` and other compilers). More ``-I`` flags might be needed
if the linear algebra framework you choose is not visible to the compiler by
default. Similarly, adding ``-L<MPI_ROOT>/lib -lmpi`` may be necessary if you
are using an MPI implementation installed in a non-system location.

.. note::

  When linking to the static version of ARPACK-NG library
  (``libarpack.a``), you might also need to explicitly link your executable to
  a BLAS/LAPACK implementation. A common symptom of this problem is
  undefined references to ``dgemv_`` and other LAPACK subroutines.

CMake
-----

Here is a minimal example of a root ``CMakeLists.txt`` script for your
project.

.. code:: cmake

  cmake_minimum_required(VERSION 3.13.0 FATAL_ERROR)

  project(myproject LANGUAGES CXX)
  set(CMAKE_CXX_STANDARD 11)

  # EZARPACK_ROOT is the installation prefix of ezARPACK
  set(ezARPACK_DIR ${EZARPACK_ROOT}/lib/cmake)

  # Import ezARPACK target
  find_package(ezARPACK 1.0 CONFIG REQUIRED)

  # Import Eigen (Blaze, Armadillo, etc) targets
  find_package(Eigen3 CONFIG)

  # Build an executable called 'myprog'
  add_executable(myprog myprog.cpp)
  target_link_libraries(myprog PRIVATE ezarpack Eigen3::Eigen)

  # ARPACK_NG_ROOT is the installation prefix of ARPACK-NG
  set(arpack-ng_DIR ${ARPACK_NG_ROOT}/lib/cmake)

  # Find ARPACK-NG
  find_package(arpack-ng 3.6.0 REQUIRED)

  # Link the executable to ARPACK-NG
  if(TARGET ARPACK::ARPACK) # ARPACK-NG 3.8.0 and newer
    target_link_libraries(myprog ARPACK::ARPACK)
  else() # ARPACK-NG prior to version 3.8.0
    target_link_libraries(myprog ${arpack_ng_LIBRARIES})
  endif()

In order to make the :ref:`Parallel ARPACK solvers <mpi>` work, ``myprog``
must be additionally linked to MPI and PARPACK (only for ARPACK-NG>=3.8.0) .

.. code:: cmake

  # Find MPI-3.0 or newer
  find_package(MPI 3.0)

  # Add MPI include directories and link to MPI libraries
  target_include_directories(myprog ${MPI_CXX_INCLUDE_PATH})
  target_link_libraries(myprog ${MPI_CXX_LIBRARIES})

  if(TARGET PARPACK::PARPACK) # ARPACK-NG 3.8.0 and newer
    target_link_libraries(myprog PARPACK::PARPACK)
  endif()
