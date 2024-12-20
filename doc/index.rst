ezARPACK
========

ezARPACK is a C++11 wrapper around `ARPACK-NG
<https://github.com/opencollab/arpack-ng>`_, a FORTRAN 77 library designed to
solve large scale sparse eigenproblems. It allows for solving
eigenproblems for symmetric, asymmetric and complex double precision
matrices with a minimal amount of boilerplate code. Both standard and
generalized problems are supported, as well as all advanced spectral
transformation modes from the original FORTRAN implementation.

The main goal of this header-only library is providing easy access to the
powerful ARPACK and Parallel ARPACK (PARPACK) algorithms from modern C++ code.
ezARPACK pursues this goal by implementing two crucial features.

- Following the ideology of ARPACK, ezARPACK does not force the user to
  stick with a predefined matrix storage format. Instead, user code is expected
  to apply a linear operator to a vector (memory buffer) passed to it and to
  store the result in another buffer at each iteration of
  the Arnoldi/Lanczos algorithm. Any callable C++ object with a suitable
  signature can be used to implement action of the operator on a vector.

- Input, output and temporary data arrays used by ARPACK subroutines are
  allocated and stored as matrix/vector objects from a supported C++ matrix
  algebra library. Currently, there are a few supported *storage backends*.

  - `Eigen 3 <http://eigen.tuxfamily.org>`_
  - `Blaze >= 3 <https://bitbucket.org/blaze-lib/blaze>`_
  - `Armadillo <http://arma.sourceforge.net>`_
  - `Boost uBLAS >= 1.58 \
    <https://www.boost.org/doc/libs/1_58_0/libs/numeric/ublas/doc>`_
  - `TRIQS arrays (part of the TRIQS libraries 2.0--3.0.x) \
    <https://triqs.github.io/triqs/3.0.x/>`_
  - `TRIQS/nda library >= 1.1.0 <https://github.com/triqs/nda/>`_
  - `xtensor >= 0.20.0 <https://github.com/xtensor-stack/xtensor>`_

  Upon choosing the right backend via a template parameter of the solver class
  :cpp:type:`ezarpack::arpack_solver`
  (or :cpp:type:`ezarpack::mpi::arpack_solver`), programmer can use callable
  objects to act on the vectors from their library of choice. The output
  eigenvalues and eigenvectors are accessible to the calling code in
  a compatible format too.

  Besides covering a range of popular matrix frameworks, ezARPACK
  allows for easy :ref:`addition of new backends <new_backend>` via
  specialization of a traits structure.

You can learn about library's capabilities and the most common usage patterns
from the :ref:`guide`. For an in-depth description of ezARPACK's API
consult the :ref:`reference`. The source code repository and the issue tracker
are hosted on `GitHub <https://github.com/krivenko/ezARPACK>`_.

Being a header-only library, ezARPACK does not require installation. However,
you may still choose to follow the :ref:`installation` to build and run unit
tests, build examples, and install the library so that it is discoverable by
other CMake-based projects.

ezARPACK is distributed under the terms of the *Mozilla Public License, v. 2.0*.
You can obtain a copy of the MPL at http://mozilla.org/MPL/2.0/.

Contents
========

.. toctree::
    :name: mastertoc
    :maxdepth: 3

    installation
    guide/index
    reference/index
    changelog
    genindex
    search
