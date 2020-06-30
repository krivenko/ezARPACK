.. _restarting:

Advanced: Customization of Lanczos/Arnoldi implicit restarting
==============================================================

ARPACK-NG implements a version of Lanczos/Arnoldi process called the
Implicitly Restarted Lanczos (Arnoldi) Iteration. At each step of this scheme,
a complete Arnoldi factorization of size :math:`m = k + p` is compressed to
size :math:`k` with :math:`k` being the requested number of eigenvalues.
The algorithm analyzes spectrum (*Ritz values*)
of current tridiagonal/upper Hessenberg representation of the matrix in
question, as well as *Ritz estimates*, which show how far the Ritz values are
from matrix's true eigenvalues. By default, the Ritz values are sorted into
'wanted' and 'unwanted' subsets. The unwanted values are then used as 'shifts'
in a basis compression process, which helps filter unwanted information and
enhance the components of the basis vectors belonging to the relevant
subspace. A complete description of this restarting scheme can be found in
`Section 4.4 of ARPACK Users' Guide \
<https://www.caam.rice.edu/software/ARPACK/UG/node50.html>`_.

In some (rather rare) cases, one may want to deviate from the default
"Exact Shift Strategy" and change the way those shifts are selected/computed.
ezARPACK allows for customization of the implicit restarting via an optional
extra argument of ``operator()`` in specializations of
:ref:`ezarpack::arpack_worker <refworkerbase>`. The argument is expected to be
a callable object ``shifts_f`` that supplies a list of shifts. It must have one
of the following signatures.

* For :ref:`ezarpack::arpack_worker\<Symmetric, Backend\> <refworkersymmetric>`:

.. code:: cpp

  using worker_t = ezarpack::arpack_worker<Symmetric, Backend>;

  // ritz_values - View of a real vector with current m Ritz values.
  // ritz_bounds - View of a real vector with current m Ritz estimates.
  // shifts - Real vector view to receive the computed shifts.
  shifts_f(worker_t::real_vector_const_view_t ritz_values,
           worker_t::real_vector_const_view_t ritz_bounds,
           worker_t::real_vector_view_t shifts);

* For
  :ref:`ezarpack::arpack_worker\<Asymmetric, Backend\> <refworkerasymmetric>`:

.. code:: cpp

  using worker_t = ezarpack::arpack_worker<Asymmetric, Backend>;

  // ritz_values_re - View of a vector with real parts of current m Ritz values.
  // ritz_values_im - View of a vector with imaginary parts of
  // current m Ritz values.
  // ritz_bounds - View of a real vector with current m Ritz estimates.
  // shifts_re - Real vector view to receive real parts of the computed shifts.
  // shifts_im - Real vector view to receive imaginary parts of the computed
  // shifts.
  shifts_f(worker_t::real_vector_const_view_t ritz_values_re,
           worker_t::real_vector_const_view_t ritz_values_im,
           worker_t::real_vector_const_view_t ritz_bounds,
           worker_t::real_vector_view_t shifts_re,
           worker_t::real_vector_view_t shifts_im);

* For :ref:`ezarpack::arpack_worker\<Complex, Backend\> <refworkercomplex>`:

.. code:: cpp

  using worker_t = ezarpack::arpack_worker<Complex, Backend>;

  // ritz_values - View of a complex vector with current m Ritz values.
  // ritz_bounds - View of a complex vector with current m Ritz estimates.
  // shifts - Complex vector view to receive the computed shifts.
  shifts_f(worker_t::complex_vector_const_view_t ritz_values,
           worker_t::complex_vector_const_view_t ritz_bounds,
           worker_t::complex_vector_view_t shifts);

