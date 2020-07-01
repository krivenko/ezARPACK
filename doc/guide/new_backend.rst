.. _new_backend:

Advanced: Adding a new storage backend
======================================

ezARPACK offers support for a few numerical linear algebra frameworks
:ref:`out of the box <backends>`. Thanks to the decoupled design of ezARPACK's
solver classes, it is fairly easy to extend the support to other matrix/vector
algebra libraries. In order to be conceptually compatible with ezARPACK's
architecture, a C++ library has to implement the following abstractions.

* A class that wraps a one-dimensional contiguous array of ``double`` elements
  (*real vector type*).
* A class that wraps a one-dimensional contiguous array of
  ``std::complex<double>`` elements (*complex vector type*).
* A class that wraps a one-dimensional contiguous array of ``int`` elements
  (*integer vector type*).
* A class that wraps a two-dimensional
  `column-major <https://en.wikipedia.org/wiki/Row-_and_column-major_order>`_
  contiguous array of ``double`` elements (*real matrix type*).
* A class that wraps a two-dimensional
  `column-major <https://en.wikipedia.org/wiki/Row-_and_column-major_order>`_
  contiguous array of ``std::complex<double>`` elements (*complex matrix type*).
* A class that implements a partial contiguous view of a real/complex vector,
  i.e. a subvector.
* A class that implements a rectangular contiguous view of a real/complex
  matrix, i.e. a matrix block. It is sufficient to support only the blocks
  including all matrix rows and a number of the leftmost columns.
* All vector and matrix types are subject to two additional requirements,

  - Their underlying data arrays must be resizable at run time;
  - There must be a way to acquire a pointer to the underlying data array.

Let us say one has a fictitious library ``mylib`` that meets all listed
requirements. One lets ezARPACK know about the new library by implementing a
*storage backend* in a new header file. In the following sections,
we give step-by-step instructions on how to write such a header.

For a complete example of a storage backend implementation, see
``<ezarpack/storages/eigen.hpp>`` (:ref:`Eigen 3 <refeigen>` backend).

Basic header structure
~~~~~~~~~~~~~~~~~~~~~~

Our storage backend header must include ``<complex>``,
``<ezarpack/storages/base.hpp>`` and all relevant headers of ``mylib``.
Normally, it is a good idea to also use some include guard macros of
``#pragma once``.

.. code:: cpp

  // mylib_storage.hpp - Storage backend for mylib.

  #pragma once

  #include <complex>
  #include <ezarpack/storages/base.hpp>

  // Include all relevant parts of 'mylib'.
  #include <mylib.hpp>

Every storage backend must define a tag type associated with it. The actual type
does not really matter, so it can be an empty structure.

.. code:: cpp

  // mylib storage backend tag.
  struct mylib_storage {};

Passing the defined tag type as the second template parameter to
``arpack_solver`` (e.g. :ref:`arpack_solver\<Symmetric,mylib_storage\>
<refsolversymmetric>`) will enable use of the ``mylib`` backend in the solver.

The most crucial and biggest part of the header is a specialization of the
``storage_traits`` structure. This traits structure is going to be a 'glue'
layer between :ref:`arpack_solver <refsolverbase>` and the new library.
``arpack_solver`` extracts ``mylib``-specific type
information from typedef members of ``storage_traits`` and calls its static
member functions to handle ``mylib``'s vector/matrix/view objects.

.. code:: cpp

  template<> struct storage_traits<mylib_storage> {
    // Member typedefs.
    ...
    // Static members functions.
    ...
  };

The rest of this HOWTO gives a detailed description of mandatory and optional
members of the ``storage_traits`` specialization.

Member type definitions of the traits structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The member type definitions of ``storage_traits`` form two groups. The first one
is for the container (vector/matrix) types.

.. code:: cpp

  // One-dimensional wrapper around a contiguous array of 'double'
  using real_vector_type = mylib::vector<double>;

  // One-dimensional wrapper around a contiguous array of
  // 'std::complex<double>'.
  using complex_vector_type = mylib::vector<std::complex<double>>;

  // One-dimensional wrapper around a contiguous array of 'int'.
  using int_vector_type = mylib::vector<int>;

  // Two-dimensional wrapper around a contiguous array of 'double'.
  // The storage order must be column-major.
  using real_matrix_type = mylib::matrix<double>;

  // Two-dimensional wrapper around a contiguous array of
  // 'std::complex<double>'. The storage order must be column-major.
  using complex_matrix_type = mylib::matrix<std::complex<double>>;

The second group includes all *view* type declarations. ezARPACK makes a
distinction between constant views and regular (read/write) views. The constant
views are returned/passed to the user code whenever a data array is meant to be
read and should be protected against external modifications. Although not
recommended, it is still possible to use ``mylib``'s read/write views as a
substitute for the constant views. This will result in functional albeit more
error-prone user code.

.. code:: cpp

  // Contiguous partial view of a real vector (subvector).
  using real_vector_view_type = mylib::vector_view<double>;

  // Contiguous partial constant view of a real vector (subvector).
  using real_vector_const_view_type = mylib::vector_cview<double>;

  // Contiguous partial view of a complex vector (subvector).
  using complex_vector_view_type = mylib::vector_view<std::complex<double>>;

  // Contiguous partial constant view of a complex vector (subvector).
  using complex_vector_const_view_type =
    mylib::vector_cview<std::complex<double>>;

  // Contiguous partial constant view of a real matrix (matrix block) that
  // includes a number of the leftmost columns.
  using real_matrix_const_view_type = mylib::matrix_cview<double>;

  // Contiguous partial constant view of a complex matrix (matrix block) that
  // includes a number of the leftmost columns.
  using complex_matrix_const_view_type =
    mylib::matrix_cview<std::complex<double>>;

Static member functions of the traits structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following member functions are mandatory for any specialization of
``storage_traits``.

* Vector object factories.

  .. code:: cpp

      // Make a real vector of a given size.
      static real_vector_type make_real_vector(int size) {
        // Call real_vector_type's constructor and return the result.
      }
      // Make a complex vector of a given size.
      static complex_vector_type make_complex_vector(int size) {
        // Call complex_vector_type's constructor and return the result.
      }
      // Make an integer vector of a given size.
      static int_vector_type make_int_vector(int size) {
        // Call int_vector_type's constructor and return the result.
      }

* Matrix object factories.

  .. code:: cpp

      // Make a real matrix with given dimensions.
      static real_matrix_type make_real_matrix(int rows, int cols) {
        // Call real_matrix_type's constructor and return the result.
      }
      // Make a complex matrix with given dimensions.
      static complex_matrix_type make_complex_matrix(int rows, int cols) {
        // Call complex_matrix_type's constructor and return the result.
      }

* Vector/matrix destructors.

  .. code:: cpp

    // Free memory occupied by vector 'v'
    static void destroy(real_vector_type & v) { ... }
    static void destroy(complex_vector_type & v) { ... }
    static void destroy(int_vector_type & v) { ... }

    // Free memory occupied by matrix 'm'
    static void destroy(real_matrix_type & m) { ... }
    static void destroy(complex_matrix_type & m) { ... }

  .. warning::

    The ``destroy()`` functions should free memory occupied by
    ``v`` and ``m`` **if and only if** it is not already done by
    ``v``'s and ``m``'s destructors. Since most libraries manage the memory on
    their own, bodies of ``destroy()`` should normally be left blank.

* Resize functions.

  .. code:: cpp

    // Resize vector 'v'.
    static void resize(real_vector_type & v, int new_size) { ... }
    static void resize(complex_vector_type & v, int new_size) { ... }
    static void resize(int_vector_type & v, int new_size) { ... }

    // Resize matrix 'm'.
    static void resize(real_matrix_type & m, int new_rows, int new_cols) {
      ...
    }
    static void resize(complex_matrix_type & m, int new_rows, int new_cols) {
      ...
    }

* Raw memory pointer accessors.

  .. code:: cpp

    // Return a pointer to the underlying data array owned by vector 'v'.
    static double* get_data_ptr(real_vector_type & v) { ... }
    static std::complex<double>* get_data_ptr(complex_vector_type &v) { ... }
    static int* get_data_ptr(int_vector_type & v) { ... }

    // Return a pointer to the underlying data array owned by matrix 'm'.
    static double* get_data_ptr(real_matrix_type & m) { ... }
    static std::complex<double>* get_data_ptr(complex_matrix_type & m) { ... }

* Vector view factories.

  .. code:: cpp

    // Make a complete view of vector 'v'.
    static real_vector_view_type make_vector_view(real_vector_type & v) {
      // Call real_vector_view_type's constructor and return the result.
    }
    static complex_vector_view_type make_vector_view(complex_vector_type & v) {
      // Call complex_vector_view_type's constructor and return the result.
    }

    // Make a partial view of vector 'v' starting at position 'start' and
    // including 'size' elements.
    static real_vector_view_type
    make_vector_view(real_vector_type & v, int start, int size) {
      // Call real_vector_view_type's constructor and return the result.
    }
    static complex_vector_view_type
    make_vector_view(complex_vector_type & v, int start, int size) {
      // Call complex_vector_view_type's constructor and return the result.
    }

    // Make a constant partial view of vector 'v' starting at position 'start'
    // and including 'size' elements.
    static real_vector_const_view_type
    make_vector_const_view(real_vector_type const& v, int start, int size) {
      // Call real_vector_const_view_type's constructor and return the result.
    }
    static complex_vector_const_view_type
    make_vector_const_view(complex_vector_type const& v, int start, int size) {
      // Call complex_vector_const_view_type's constructor and return the result
    }

* Matrix constant view factories.

  .. code:: cpp

    // Make a complete constant view of matrix 'm'.
    static real_matrix_const_view_type
    make_matrix_const_view(real_matrix_type const& m) {
      // Call real_matrix_const_view_type's constructor and return the result.
    }
    static complex_matrix_const_view_type
    make_matrix_const_view(complex_matrix_type const& m) {
      // Call complex_matrix_const_view_type's constructor and return the result
    }

    // Make a partial constant view of matrix 'm' including 'cols'
    // leftmost columns.
    static real_matrix_const_view_type
    make_matrix_const_view(real_matrix_type const& m, int rows, int cols) {
      // Call real_matrix_const_view_type's constructor and return the result.
    }
    static complex_matrix_const_view_type
    make_matrix_const_view(complex_matrix_type const& m, int rows, int cols) {
      // Call complex_matrix_const_view_type's constructor and return the result
    }

Some of the functions, such as ``destroy()`` and ``resize()``, do not have to be
defined separately for each argument type. It is acceptable to use function
templates instead.

With these functions implemented, one can already instantiate and use
:ref:`arpack_solver\<Symmetric,mylib_storage\>
<refsolversymmetric>` and
:ref:`arpack_solver\<Complex,mylib_storage\>
<refsolvercomplex>`. The asymmetric case, however, requires more work, as
described in the next section.

Optional: Eigenvalue/eigenvector post-processing functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Because of specifics of the internal data storage format and numerical
algorithm, extracting eigenvalues and eigenvectors after a completed
:ref:`arpack_solver\<Asymmetric,mylib_storage\> <refsolverasymmetric>` run needs
some post-processing that is not done by ARPACK-NG itself.
The storage traits structure may optionally implement three static
member functions, which will be called by the asymmetric solver to extract a
computed eigensystem from memory buffers and return it to the user in
a convenient form.

.. code:: cpp

  static complex_vector_type
  make_asymm_eigenvalues(real_vector_type const& dr,
                         real_vector_type const& di,
                         int nconv) {
    // Compute and return dr + i*di
  }

``make_asymm_eigenvalues()`` is the simplest of the three functions. It is
called to combine two real vectors -- lists of real (``dr``) and
imaginary (``di``) parts of computed eigenvalues -- into one complex vector.
``nconv`` is the total number of the computed eigenvalues. Exactly ``nconv``
first elements of ``dr`` and ``di`` must be used to form the complex vector
(``dr`` and ``di`` can, in general, be longer or not providing size information
at all).

.. code:: cpp

    inline static complex_matrix_type
    make_asymm_eigenvectors(real_vector_type const& z,
                            real_vector_type const& di,
                            int N,
                            int nconv) {
      // Extract and return eigenvectors as columns of a complex matrix.
    }

This function extracts eigenvectors from a real vector ``z`` according to
special rules. ``z`` holds components of the eigenvectors as a sequence of
``nconv`` length-``N`` chunks, where ``N`` is the dimension of the eigenproblem
and ``nconv`` has the same meaning as before. Meaning of each chunk depends on
the corresponding component of ``di``. If ``di[i]`` is zero, then the ``i``-th
chunk of ``z`` contains a real eigenvector. Otherwise,
``di[i] = -di[i+1] != 0``, in which case the ``i``-th and ``(i+1)``-th chunks
of ``z`` are real and imaginary parts of a complex eigenvector respectively.
Every such pair corresponds to a complex conjugate pair of eigenvectors,
so that the total amount of vectors stored in ``z`` is exactly ``nconv``.
The extracted eigenvectors must be returned as columns of a complex
``N`` x ``nconv`` matrix.

.. code:: cpp

    template<typename A>
    inline static complex_vector_type
    make_asymm_eigenvalues(real_vector_type const& z,
                          real_vector_type const& di,
                          A&& a,
                          int N,
                          int nconv) {
      // Compute eigenvalues as Rayleigh quotients and return them in
      // a complex vector.
    }

In the ``ShiftAndInvertReal`` and ``ShiftAndInvertImag`` spectral transformation
modes, ARPACK-NG computes eigenvalues of an auxiliary real matrix. Those
eigenvalues are implicitly related to the ones of the original eigenproblem.
One way to extract the original eigenvalues is via solution of a quadratic
equation. Unfortunately, this approach is not perfect, because solutions
of quadratic equations are not unique, and it can be difficult to match the
correct solution with a given eigenvector :math:`\mathbf{x}`. A robust
alternative approach is to compute the eigenvalue :math:`\lambda` of
:math:`\hat A\mathbf{x} = \lambda\hat M\mathbf{x}` as the Rayleigh quotient
:math:`\lambda = \frac{\mathbf{x}^\dagger \hat A \mathbf{x}}
{\mathbf{x}^\dagger\hat M \mathbf{x}}`, which is the purpose of the last of the
three functions. ``z``, ``di``, ``N`` and ``nconv`` have the same meaning as
before, and callable object ``a`` represents the linear operator :math:`\hat A`.
This overload of ``make_asymm_eigenvalues()`` should extract the eigenvectors
from ``z`` one by one and compute :math:`\lambda` for each of them. It is
beneficial to treat the real vectors differently from the complex ones,
as the Rayleigh quotient can be computed at lower memory and CPU costs
if :math:`\mathbf{x}^\dagger = \mathbf{x}^T`.

.. note:: Despite the name, the quotients amount to just the numerators.
          ARPACK-NG guarantees that :math:`\mathbf{x}^\dagger\hat M
          \mathbf{x} = 1`, so there is no need to consider matrix
          :math:`\hat M` at all.
