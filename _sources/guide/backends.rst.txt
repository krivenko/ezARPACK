.. _backends:

List of supported matrix/vector algebra libraries
=================================================

One can choose which of supported linear algebra libraries to use by passing
a special tag type as the second template parameter (``Backend``) to
the :ref:`ezarpack::arpack_solver <refsolverbase>` class template.

The following table lists all supported backends, their respective tag types
and headers containing backends' definitions. The last columns contains links
to complete usage examples.

.. list-table::
  :header-rows: 1
  :align: left
  :widths: auto

  * - Storage backend
    - Tag type
    - Header
    - Examples

  * - `Eigen 3 <http://eigen.tuxfamily.org>`_
    - ``ezarpack::eigen_storage``
    - ``<ezarpack/storages/eigen.hpp>``
    - :ref:`example_eigen`

  * - `Blaze >= 3 <https://bitbucket.org/blaze-lib/blaze>`_
    - ``ezarpack::blaze_storage``
    - ``<ezarpack/storages/blaze.hpp>``
    - :ref:`example_blaze`

  * - `Armadillo <http://arma.sourceforge.net>`_
    - ``ezarpack::armadillo_storage``
    - ``<ezarpack/storages/armadillo.hpp>``
    - :ref:`example_armadillo`

  * - `Boost uBLAS >= 1.58 \
      <https://www.boost.org/doc/libs/1_58_0/libs/numeric/ublas/doc>`_
    - ``ezarpack::ublas_storage``
    - ``<ezarpack/storages/ublas.hpp>``
    - :ref:`example_ublas`

  * - `TRIQS arrays >= 2.0 <https://triqs.github.io/triqs/latest/>`_
    - ``ezarpack::triqs_storage``
    - ``<ezarpack/storages/triqs.hpp>``
    - :ref:`example_triqs`

  * -  `xtensor >= 0.20.0 <https://github.com/xtensor-stack/xtensor>`_
    - ``ezarpack::xtensor_storage``
    - ``<ezarpack/storages/xtensor.hpp>``
    - :ref:`example_xtensor`
