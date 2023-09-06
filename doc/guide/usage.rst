.. _usage:

How to use ezARPACK in your project
===================================

Depending on what build system your C++ project is based on, there are
two ways to use ezARPACK. Either way, you would need

* A compiler supporting C++11 or newer;
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

  g++ -O3 -o myprog myprog.cpp \
      -I<EZARPACK_ROOT>/include -L<ARPACK_NG_ROOT>/lib -larpack

  # For executables using the Parallel ARPACK (MPI) solvers
  mpic++ -O3 -o myprog_mpi myprog_mpi.cpp -I<EZARPACK_ROOT>/include \
         -L<ARPACK_NG_ROOT>/lib -larpack -lparpack

(similar for ``clang++`` and other compilers). More ``-I`` flags may be needed
if the linear algebra framework you choose is not visible to the compiler by
default. Similarly, adding ``-I<MPI_ROOT>/include -L<MPI_ROOT>/lib`` may be
necessary if you are using an MPI implementation installed in a non-system
location.

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
  find_arpackng(3.6.0 REQUIRED)

  # Link the executable to the ARPACK library.
  target_link_libraries(myprog ${ARPACK_LIBRARIES})

Linking your targets to ARPACK-NG libraries can be a bit of a hassle, as CMake
interface of ARPACK-NG changed a few times since version 3.6.0. In particular,
CMake scripts of the versions prior to 3.8.0 do not export any targets. Instead,
they expose library information via a variable ``arpack_ng_LIBRARIES``.
In order to make linking more user-friendly, ezARPACK exports a macro called
``find_arpackng()``. It finds an ARPACK-NG installation while taking care of
said CMake interface differences. Upon successful detection of ARPACK-NG,
``find_arpackng()`` sets two variables that can be later passed to
``target_link_libraries()``,

  - ``ARPACK_LIBRARIES`` - list of ARPACK libraries and linker flags
  - ``PARPACK_LIBRARIES`` - list of Parallel ARPACK libraries and linker flags

Setting CMake variable ``ARPACK_NG_ROOT`` will instruct ``find_arpackng()``
to look for ARPACK-NG at a specific installation prefix before proceeding
to system locations.

Making the :ref:`Parallel ARPACK solvers <mpi>` work requires additional
linking to MPI and PARPACK.

.. code:: cmake

  # Parallel ARPACK (MPI)

  # Build another executable `myprog_mpi`.
  add_executable(myprog_mpi myprog_mpi.cpp)
  target_link_libraries(myprog_mpi ezarpack Eigen3::Eigen)

  # Detect an MPI-3.0 implementation.
  find_package(MPI 3.0 REQUIRED)

  # Link the executable to the Parallel ARPACK library and to the MPI.
  target_include_directories(myprog_mpi PRIVATE ${MPI_CXX_INCLUDE_PATH})
  target_link_libraries(myprog_mpi ${PARPACK_LIBRARIES} ${MPI_CXX_LIBRARIES})
