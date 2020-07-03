.. _complex:

General complex eigenproblems
=============================

This page is a walkthrough showing how to use
:ref:`ezarpack::arpack_solver\<Complex, \Backend> <refsolvercomplex>` in
your C++ code to compute a few eigenpairs :math:`(\lambda,\mathbf{x})` of

.. math::

  \hat A  \mathbf{x} = \lambda \hat M \mathbf{x}

with a complex matrix :math:`\hat A` and a complex Hermitian positive
semi-definite matrix :math:`\hat M`. The complex solver class supports a few
computational modes, where the original eigenproblem is recast into

.. math::

  \hat O \mathbf{x} = \mu \mathbf{x}.

Operator :math:`\hat O` acts in a vector space equipped with an inner product
defined by the matrix :math:`\hat B = \hat M`,

.. math::

  \langle \mathbf{x}, \mathbf{y} \rangle = \mathbf{x}^\dagger \hat B \mathbf{y}.

There are explicit relations between the original eigenvalues :math:`\lambda`
and their transformed counterparts :math:`\mu`.

.. note::

  If both matrices :math:`\hat A` and :math:`\hat M` are real,
  :ref:`ezarpack::arpack_solver\<Symmetric, \Backend> <refsolversymmetric>` or
  :ref:`ezarpack::arpack_solver\<Asymmetric, \Backend> <refsolverasymmetric>`
  should be used instead.

Typical steps needed to compute the eigenpairs are as follows.

1. Decide what :ref:`storage backend<backends>` you want to use or whether it is
   appropriate to :ref:`implement a new one<new_backend>`. In the following, we
   will assume that the `Eigen <http://eigen.tuxfamily.org>`_ backend has been
   selected.

2. Include ``<ezarpack/arpack_solver.hpp>`` and the relevant backend header.

  .. code::

    #include <ezarpack/arpack_solver.hpp>
    #include <ezarpack/storages/eigen.hpp>

  .. note::

    ``<ezarpack/arpack_solver.hpp>`` includes
    :ref:`all three specializations of ezarpack::arpack_solver<refsolver>`
    at once. If you want to speed up compilation a little bit, you can
    include ``<ezarpack/solver_base.hpp>`` and
    ``<ezarpack/solver_complex.hpp>`` instead.

3. Create a solver object.

  .. code:: cpp

    const int N = 1000; // Size of matrices A and M

    using namespace ezarpack;

    // Shorthand for solver's type.
    using solver_t = arpack_solver<Complex, eigen_storage>;
    // Solver object.
    solver_t solver(N);

4. Fill a ``params_t`` structure with calculation parameters.

  .. code:: cpp

    // params_t is a structure holding parameters of
    // the Implicitly Restarted Arnoldi iteration.
    using params_t = solver_t::params_t;

    // Requested number of eigenvalues to compute.
    const int nev = 10;

    // Compute the eigenvalues with the largest magnitudes.
    auto eigenvalues_select = params_t::LargestMagnitude;

    // Compute Ritz vectors (eigenvectors).
    params_t::compute_vectors_t compute_vectors = params_t::Ritz;

    // params_t's constructor takes three arguments -- mandatory parameters
    // that need be set explicitly.
    params_t params(nev, eigenvalues_select, compute_vectors);

  The following table contains an annotated list of all supported parameters.

  .. _params:

  .. list-table::
    :header-rows: 1
    :align: left
    :widths: auto

    * - Parameter name
      - Type
      - Default value
      - Description

        .. _complex_n_eigenvalues:
    * - ``n_eigenvalues``
      - ``unsigned int``
      - n/a
      - Number of eigenvalues to compute.

        .. _complex_eigenvalues_select:
    * - ``eigenvalues_select``
      - ``params_t::eigenvalues_select_t`` (enumeration)
      - n/a
      - Part of the spectrum to target. Acceptable values are
        ``LargestMagnitude`` (largest eigenvalues in magnitude),
        ``SmallestMagnitude`` (smallest eigenvalues in magnitude),
        ``LargestReal`` (eigenvalues of largest real part),
        ``SmallestReal`` (eigenvalues of smallest real part),
        ``LargestImag`` (eigenvalues of largest imaginary part) and
        ``SmallestImag`` (eigenvalues of smallest imaginary part).

        .. _complex_ncv:
    * - ``ncv``
      - ``int``
      - min(2 * ``n_eigenvalues`` + 2, ``N``)
      - How many Arnoldi vectors to generate at each iteration.

        .. _complex_compute_vectors:
    * - ``compute_vectors``
      - ``compute_vectors_t`` (enumeration)
      - n/a
      - ``Schur`` -- compute only Schur vectors (orthogonal basis vectors of
        the ``n_eigenvalues``-dimensional subspace), ``Ritz`` -- compute Ritz
        vectors (eigenvectors) in addition to the Schur vectors, ``None`` --
        compute neither Schur nor Ritz vectors.

        .. _complex_random_residual_vector:
    * - ``random_residual_vector``
      - ``bool``
      - ``true``
      - Use a randomly generated initial residual vector?

        .. _complex_sigma:
    * - ``sigma``
      - ``std::complex<double>``
      - `0`
      - Complex eigenvalue shift :math:`\sigma` for spectral transformation
        modes.

        .. _complex_tolerance:
    * - ``tolerance``
      - ``double``
      - Machine precision
      - Relative tolerance for Ritz value (eigenvalue) convergence.

        .. _complex_max_iter:
    * - ``max_iter``
      - ``unsigned int``
      - ``INT_MAX``
      - Maximum number of Arnoldi update iterations allowed.

  .. note::

    In the Shift-and-Invert mode, values of
    :ref:`eigenvalues_select <complex_eigenvalues_select>` refer
    to the spectrum of the **transformed** problem, not the original one. For
    instance, ``LargestMagnitude`` with a complex shift :math:`\sigma`
    will pick eigenvalues :math:`\lambda` closest to
    :math:`\sigma`, because they correspond to the eigenvalues
    :math:`\mu = 1/(\lambda - \sigma)` that have the largest magnitude.

5. Optionally set the initial vector for Arnoldi iteration if a better choice
   than a random vector is known.
   :ref:`random_residual_vector <complex_random_residual_vector>` parameter must
   be set to ``false`` for the changes made to the initial vector
   to take effect.

   A view of the residual vector is accessible via method
   ``residual_vector()`` of the solver.

   .. code:: cpp

     // Set all components of the initial vector to 1.
     auto rv = solver.residual_vector();
     for(int i = 0; i < N; ++i) rv[i] = 1.0;

   One may also call ``residual_vector()`` later, after a diagonalization run
   has started, to retrieve the current residual vector.

6. Choose one of supported computational modes and perform diagonalization.
   In this part, user is supposed to call the ``solver`` object and pass the
   parameter structure as well as callable objects (*e.g.* lambda-functions)
   that represent action of operators :math:`\hat O` and :math:`\hat B` on
   a given vector. The supplied objects will be called to generate Arnoldi
   vectors. Syntax and semantics of the C++ code vary between
   the computational modes and will be explained individually for each of
   them.

     .. _complex_standard:

   - **Standard mode** (for standard eigenproblems, :math:`\hat M = \hat I`).

     .. code:: cpp

       using vector_view_t = solver_t::vector_view_t;
       using vector_const_view_t = solver_t::vector_const_view_t;

       auto Aop = [](vector_const_view_t in, vector_view_t out) {
         // Code implementing action of matrix A on vector 'in':
         // out = A * in
       };

       solver(Aop, params);

     .. _complex_Inverse:

   - **Regular inverse mode** (for positive-definite :math:`\hat M`).

     In this mode, the transformed eigenproblem is defined by
     :math:`\hat O = \hat M^{-1} \hat A`, :math:`\hat B = \hat M` and
     :math:`\lambda = \mu`.

     .. code:: cpp

       using vector_view_t = solver_t::vector_view_t;
       using vector_const_view_t = solver_t::vector_const_view_t;

       auto op = [](vector_const_view_t in, vector_view_t out) {
         // Code implementing action of matrix M^{-1} A on vector 'in':
         // out = M^{-1} * A * in
       };
       auto Bop = [](vector_const_view_t in, vector_view_t out) {
         // Code implementing action of matrix M on vector 'in':
         // out = M * in
       };

       solver(op, Bop, solver_t::Inverse, params);

     Inverting a sparse matrix :math:`\hat M` will likely make it dense, which
     is usually undesirable from the storage standpoint. A more practical
     solution is to compute the sparse LU or Cholesky factorization of
     :math:`\hat M` once (outside of the lambda-function's body), and write
     the lambda-function so that it computes ``out`` as the solution of
     the linear system ``M * out = A * in`` using the precomputed factorization.

     .. _complex_ShiftAndInvert:

   - **Shift-and-Invert mode**.

     In this mode, the transformed eigenproblem is defined by
     :math:`\hat O = (\hat A -\sigma \hat M)^{-1} \hat M`,
     :math:`\hat B = \hat M` and :math:`\lambda = 1/\mu + \sigma`.
     The complex spectral shift :math:`\sigma` must be set in the parameters
     structure, see the :ref:`list of parameters <complex_sigma>`.

     .. code:: cpp

       using vector_view_t = solver_t::vector_view_t;
       using vector_const_view_t = solver_t::vector_const_view_t;

       auto op = [](vector_view_t in, vector_view_t out) {
         // Code implementing action of matrix (A - sigma*M)^{-1} * M on 'in':
         // out = (A - sigma*M)^{-1} * M * in
       };
       auto Bop = [](vector_const_view_t in, vector_view_t out) {
         // Code implementing action of matrix M on vector 'in':
         // out = M * in
       };

       solver(op, Bop, solver_t::ShiftAndInvert, params);

     Inverting a sparse matrix :math:`\hat A - \sigma\hat M` will likely make it
     dense, which is usually undesirable from the storage standpoint. A more
     practical solution is to compute the sparse LU or Cholesky factorization of
     :math:`\hat A - \sigma\hat M` once (outside of the lambda-function's body),
     and write the lambda-function so that it (1) computes ``M * in`` and
     (2) computes ``out`` as the solution of the linear system
     ``(A - \sigma M) * out = M * in`` using the precomputed factorization.

   .. note::

     In most computational modes above, it is seemingly necessary to apply
     operator :math:`\hat B` to the same vector twice per generated Arnoldi
     vector, once in functor ``op`` and once in ``Bop``. It is actually possible
     to spare one of the applications. Calling ``solver.Bx_available()`` inside
     ``op`` will tell whether ``Bop`` has already been called at the current
     iteration, and ``solver.Bx_vector()`` will return a constant view of the
     application result :math:`\hat B \mathbf{x}`.

   The ``in`` and ``out`` views passed to the callable objects always expose one
   of three length-:math:`N` vectors stored inside the solver object. There is
   another, indirect way to access them.

   .. code:: cpp

     // Get index (0-2) of the current 'in' vector and request a view of it
     auto in_view = solver.workspace_vector(solver.in_vector_n());
     // Similar for the 'out' vector
     auto out_view = solver.workspace_vector(solver.out_vector_n());

   In advanced usage scenarios, the implicit restarting procedure can be
   customized via an extra argument of ``solver``'s call operator.
   See ':ref:`restarting`' for more details.

   .. code:: cpp

     auto shifts_f = [](solver_t::complex_vector_const_view_t ritz_values,
                        solver_t::complex_vector_const_view_t ritz_bounds,
                        solver_t::complex_vector_view_t shifts) {
                          // Compute shifts for the implicit restarting
                        };

     // Standard mode
     solver(op, params, shifts_f);
     // Other modes, e.g. Inverse
     solver(op, Bop, solver_t::Inverse, params, shifts_f);

   ``solver_t::operator()`` can throw two special exception types.

   - :ref:`ezarpack::maxiter_reached <maxiter_reached>` - Maximum number of
     implicitly restarted Arnoldi iterations has been reached.
   - :ref:`ezarpack::ncv_insufficient <ncv_insufficient>` - No shifts could be
     applied during a cycle of the Implicitly restarted Arnoldi iteration.
     Consider increasing the number of Arnoldi vectors generated at each
     iteration (:ref:`ncv <complex_ncv>` parameter).

   The rest of possible problems reported by ARPACK-NG result in generic
   `std::runtime_error <https://en.cppreference.com/w/cpp/error/runtime_error>`_
   exceptions.

   7. Request computed eigenvalues and eigenvectors. For the eigenvectors, the
   :ref:`compute_vectors <complex_compute_vectors>` parameter must be set to
   ``params_t::Ritz``.

   .. code:: cpp

     auto lambda = solver.eigenvalues();
     auto vecs = solver.eigenvectors();

   The eigenvectors are columns of the complex matrix view ``vecs``.
   If the diagonalization run has ended prematurely (for example, when
   the maximum number of iterations has been reached), then it is still
   possible to extract ``solver.nconv()`` converged eigenpairs.

8. Optionally request the Schur vectors, i.e. :math:`\hat B`-orthogonal basis
   vectors of the relevant vector subspace
   (:ref:`compute_vectors <complex_compute_vectors>` must be either
   ``params_t::Schur`` or ``params_t::Ritz``).

   .. code:: cpp

     auto basis = solver.schur_vectors();

   The basis vectors are ``solver.nconv()`` columns of the complex matrix
   view ``basis``.

9. Optionally request statistics about the completed run.

   .. code:: cpp

     // Print some computation statistics
     auto stats = solver.stats();

     std::cout << "Number of Arnoldi update iterations: " << stats.n_iter
               << std::endl;
     std::cout << "Total number of O*x operations: " << stats.n_op_x_operations
               << std::endl;
     std::cout << "Total number of B*x operations: " << stats.n_b_x_operations
               << std::endl;
     std::cout << "Total number of steps of re-orthogonalization: "
               << stats.n_reorth_steps << std::endl;
