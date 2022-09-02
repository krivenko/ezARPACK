.. _mpi:

Parallel ARPACK versions of eigensolvers
========================================

As size :math:`N` of an eigenproblem grows, it becomes more and more beneficial
to employ
`MPI parallelism <https://en.wikipedia.org/wiki/Message_Passing_Interface>`_
offered by the
`Parallel ARPACK (PARPACK) <https://doi.org/10.1007/3-540-62095-8_51>`_ wrappers
:cpp:type:`ezarpack::mpi::arpack_solver`. Despite being very similar to the
serial versions API-wise, they have a major conceptual distinction.

All :math:`N`-dimensional vectors :cpp:type:`ezarpack::mpi::arpack_solver`
deals with are partitioned into contiguous blocks and stored in a distributed
manner: Each MPI process physically stores only one block of each such vector.
The distributed storage scheme applies to ARPACK's workspace vectors, as well
as to computed Ritz and Schur vectors. Below we list API differences of
:cpp:type:`ezarpack::mpi::arpack_solver` as compared to
:cpp:type:`ezarpack::arpack_solver`.

.. default-domain:: cpp

* Specializations of :cpp:type:`ezarpack::mpi::arpack_solver` are defined
  in header files under the ``ezarpack/mpi/`` subdirectory.

.. code-block:: cpp

    // Symmetric real eigenproblems
    #include <ezarpack/mpi/solver_symmetric.hpp>
    // General real eigenproblems
    #include <ezarpack/mpi/solver_asymmetric.hpp>
    // General complex eigenproblems
    #include <ezarpack/mpi/solver_complex.hpp>

    // Include all specializations at once
    #include <ezarpack/mpi/arpack_solver.hpp>

* Constructors take an MPI communicator.

  .. function:: arpack_solver::arpack_solver(unsigned int N, \
                                             MPI_Comm const& comm)
  .. function:: arpack_solver::arpack_solver(\
                std::vector<unsigned int> const& block_sizes, \
                MPI_Comm const& comm)

  The first constructor partitions :math:`N`-dimensional vectors and distributes
  their blocks among all processes in ``comm`` in the most even way.
  In particular, if size of the communicator :math:`P` divides :math:`N`,
  then all blocks will have the same size :math:`N/P`.

  The second constructor takes a list of :math:`P` positive numbers
  -- block sizes -- with one element per each process in ``comm``.
  It computes :math:`N` as a total of the list.

* Access to the MPI communicator used to construct a solver.

  .. function:: MPI_Comm arpack_solver::mpi_comm() const

* Access to the index range of the block stored by the calling process.

  .. function:: int arpack_solver::local_block_start() const
  .. function:: int arpack_solver::local_block_size() const

  Locally stored elements of a vector span the
  ``[local_block_start(); local_block_start() + local_block_size()[`` index
  range within the :math:`N`-dimensional space.

* Overloads of :cpp:func:`arpack_solver::operator()()` take callable objects
  that implement action of various sparse matrices on workspace vectors.
  In the parallel version of the solvers, vector views ``in`` and ``out``
  passed to those callables expose only **the locally stored blocks** of the
  full :math:`N`-dimensional vectors. Nonetheless, **the callables are still
  expected to act in the full space.** If there are non-zero matrix elements
  connecting blocks stored by different processes, it is user's responsibility
  to organize data exchange between MPI processes so that ``out`` is properly
  updated on all ranks.

* Methods :cpp:func:`arpack_solver::workspace_vector()`,
  :cpp:func:`arpack_solver::residual_vector()` and
  :cpp:func:`arpack_solver::Bx_vector()` return views of the corresponding
  locally stored blocks. Similarly, each column of the matrix views returned
  by :cpp:func:`arpack_solver::eigenvectors()` and
  :cpp:func:`arpack_solver::schur_vectors()` has length ``local_block_size()``
  and represents the locally stored block of a Ritz/Schur vector.

  Note, however, that :cpp:func:`arpack_solver::eigenvalues()` returns
  the full list of converged eigenvalues on all MPI ranks.
