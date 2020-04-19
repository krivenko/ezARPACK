.. _reference:

API reference
=============

* Specializations of the `arpack_worker` eigensolver class

    * :ref:`refworkersymmetric` - symmetric real eigenproblems
    * :ref:`refworkerasymmetric` - general real eigenproblems
    * :ref:`refworkercomplex` - general complex eigenproblems

* Data storage backends

    * :ref:`refeigen` - `Eigen 3 <http://eigen.tuxfamily.org>`_
    * :ref:`refblaze` - `Blaze >= 3 <https://bitbucket.org/blaze-lib/blaze>`_
    * :ref:`refarmadillo` - `Armadillo <http://arma.sourceforge.net>`_
    * :ref:`refublas` - `Boost uBLAS >= 1.58
      <https://www.boost.org/doc/libs/1_58_0/libs/numeric/ublas/doc>`_
    * :ref:`reftriqs` - `TRIQS arrays <https://triqs.github.io/triqs/latest/>`_
    * :ref:`refxtensor` - `xtensor <https://github.com/xtensor-stack/xtensor>`_

* Auxiliary headers

    * :ref:`refworker` - includes all `arpack_worker` variants at once
    * :ref:`refarpack` - Low level C++ interface to FORTRAN routines of ARPACK
    * :ref:`refversion` - ezARPACK version information
