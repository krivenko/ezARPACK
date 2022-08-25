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
<http://li.mit.edu/Archive/Activities/Archive/CourseWork/Ju_Li/MITCourses/18.335/Doc/ARPACK/Lehoucq97.pdf#page=64>`_.

In some (rather rare) cases, one may want to deviate from the default
"Exact Shift Strategy" and change the way those shifts are selected/computed.
ezARPACK allows for customization of the implicit restarting via an optional
extra argument of ``operator()`` in specializations of
:ref:`ezarpack::arpack_solver <refsolverbase>`. The argument is expected to be
a callable object ``shifts_f`` that supplies a list of shifts. It must have one
of the following signatures.

* For :ref:`ezarpack::arpack_solver\<Symmetric, Backend\> <refsolversymmetric>`:

.. code:: cpp

  using solver_t = ezarpack::arpack_solver<Symmetric, Backend>;

  // ritz_values - View of a real vector with current m Ritz values.
  // ritz_bounds - View of a real vector with current m Ritz estimates.
  // shifts - Real vector view to receive the computed shifts.
  shifts_f(solver_t::real_vector_const_view_t ritz_values,
           solver_t::real_vector_const_view_t ritz_bounds,
           solver_t::real_vector_view_t shifts);

* For
  :ref:`ezarpack::arpack_solver\<Asymmetric, Backend\> <refsolverasymmetric>`:

.. code:: cpp

  using solver_t = ezarpack::arpack_solver<Asymmetric, Backend>;

  // ritz_values_re - View of a vector with real parts of current m Ritz values.
  // ritz_values_im - View of a vector with imaginary parts of
  // current m Ritz values.
  // ritz_bounds - View of a real vector with current m Ritz estimates.
  // shifts_re - Real vector view to receive real parts of the computed shifts.
  // shifts_im - Real vector view to receive imaginary parts of the computed
  // shifts.
  shifts_f(solver_t::real_vector_const_view_t ritz_values_re,
           solver_t::real_vector_const_view_t ritz_values_im,
           solver_t::real_vector_const_view_t ritz_bounds,
           solver_t::real_vector_view_t shifts_re,
           solver_t::real_vector_view_t shifts_im);

* For :ref:`ezarpack::arpack_solver\<Complex, Backend\> <refsolvercomplex>`:

.. code:: cpp

  using solver_t = ezarpack::arpack_solver<Complex, Backend>;

  // ritz_values - View of a complex vector with current m Ritz values.
  // ritz_bounds - View of a complex vector with current m Ritz estimates.
  // shifts - Complex vector view to receive the computed shifts.
  shifts_f(solver_t::complex_vector_const_view_t ritz_values,
           solver_t::complex_vector_const_view_t ritz_bounds,
           solver_t::complex_vector_view_t shifts);

