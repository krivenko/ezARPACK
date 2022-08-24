.. _examples_mpi:

Usage examples: MPI-parallelized solvers
========================================

These examples show how to find a few smallest eigenvalues of a real symmetric
matrix and their respective eigenvectors using the MPI-parallelized (PARPACK)
versions of ezARPACK's eigensolvers. The matrix is a banded one with
the bandwidth equal to 5. Note that the matrix is never stored in memory and
is defined only by the rule of how it acts on a vector.

.. _example_eigen_mpi:

Eigen 3
-------

.. literalinclude:: ../../example/mpi/eigen.mpi.example.cpp
  :language: cpp
  :lines: 14-
  :linenos:

.. _example_blaze_mpi:

Blaze
-----

.. literalinclude:: ../../example/mpi/blaze.mpi.example.cpp
  :language: cpp
  :lines: 14-
  :linenos:

.. _example_armadillo_mpi:

Armadillo
---------

.. literalinclude:: ../../example/mpi/armadillo.mpi.example.cpp
  :language: cpp
  :lines: 14-
  :linenos:

.. _example_ublas_mpi:

Boost uBLAS
-----------

.. literalinclude:: ../../example/mpi/ublas.mpi.example.cpp
  :language: cpp
  :lines: 14-
  :linenos:

.. _example_triqs_mpi:

TRIQS arrays
------------

.. literalinclude:: ../../example/mpi/triqs.mpi.example.cpp
  :language: cpp
  :lines: 14-
  :linenos:

.. _example_nda_mpi:

TRIQS/nda
---------

.. literalinclude:: ../../example/mpi/nda.mpi.example.cpp
  :language: cpp
  :lines: 14-
  :linenos:

.. _example_xtensor_mpi:

xtensor
-------

.. literalinclude:: ../../example/mpi/xtensor.mpi.example.cpp
  :language: cpp
  :lines: 14-
  :linenos:
