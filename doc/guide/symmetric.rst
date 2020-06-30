.. _symmetric:

Symmetric real eigenproblems
============================

This page is a walkthrough showing how to use
:ref:`ezarpack::arpack_worker\<Symmetric, \Backend> <refworkersymmetric>` in
your C++ code to compute a few eigenpairs :math:`(\lambda,\mathbf{x})` of

.. math::

  \hat A  \mathbf{x} = \lambda \hat M \mathbf{x}

with real symmetric matrices :math:`\hat A` and :math:`\hat M`. The symmetric
worker class supports a few computational modes where the original eigenproblem
is recast into

.. math::

  \hat O \mathbf{x} = \mu \mathbf{x}.

The new matrix :math:`\hat O` is symmetric w.r.t. to an inner product
defined by a symmetric positive semi-definite matrix :math:`\hat B`,

.. math::

  \langle \mathbf{x}, \hat O \mathbf{y} \rangle =
  \langle \hat O \mathbf{x}, \mathbf{y} \rangle, \quad
  \langle \mathbf{x}, \mathbf{y} \rangle = \mathbf{x}^T \hat B \mathbf{y}

or, equivalently,

.. math::

  \hat B \hat O = \hat O^T \hat B.

There are explicit relations between the original eigenvalues :math:`\lambda`
and their transformed counterparts :math:`\mu`. :math:`\lambda` are guaranteed
to be real, and computed eigenvectors :math:`\mathbf{x}_i` form an
orthonormal system w.r.t. the :math:`\hat B`-weighted inner product,
:math:`\langle \mathbf{x}_i, \mathbf{x}_j \rangle = \delta_{i,j}`.

Typical steps needed to compute the eigenpairs are as follows.

1. Decide what :ref:`storage backend<backends>` you want to use or whether it is
   appropriate to :ref:`implement a new one<new_backend>`. In the following, we
   will assume that the `Eigen <http://eigen.tuxfamily.org>`_ backend has been
   selected.

2. Include `arpack_worker.hpp` and the relevant backend header.

  .. code::

    #include <ezarpack/arpack_worker.hpp>
    #include <ezarpack/storages/eigen.hpp>

  .. note::

    `arpack_worker.hpp` includes
    :ref:`all three specializations of ezarpack::arpack_worker<refworker>`
    at once. If you want to speed up compilation a little bit, you can
    include `worker_base.hpp` and `worker_symmetric.hpp` instead.

3. Create a worker object.

  .. code:: cpp

    const int N = 1000; // Size of matrices A and M

    using namespace ezarpack;

    // Shorthand for worker's type.
    using worker_t = arpack_worker<Symmetric, eigen_storage>;
    // Worker object.
    worker_t worker(N);

4. Fill a `params_t` structure with calculation parameters.

  .. code:: cpp

    // params_t is a structure holding parameters of
    // the Implicitly Restarted Lanczos iteration.
    using params_t = worker_t::params_t;

    // Requested number of eigenvalues to compute.
    const int nev = 10;

    // Compute the smallest eigenvalues.
    auto eigenvalues_select = params_t::Smallest;

    // Compute eigenvectors too?
    bool compute_eigenvectors = true;

    // params_t's constructor takes three arguments -- mandatory parameters
    // that need be set explicitly.
    params_t params(10, eigenvalues_select, compute_eigenvectors);

  The following table contains an annotated list of all supported parameters.

  .. list-table::
    :header-rows: 1
    :align: left
    :widths: auto

    * - Parameter name
      - Type
      - Default value
      - Description

    * - `n_eigenvalues`
      - `unsigned int`
      - n/a
      - Number of eigenvalues to compute.

    * - `eigenvalues_select`
      - `params_t::eigenvalues_select_t` (enumeration)
      - n/a
      - Part of the spectrum to target. Acceptable values are
        `Largest` (algebraically largest eigenvalues),
        `Smallest` (algebraically smallest eigenvalues),
        `LargestMagnitude` (largest eigenvalues in magnitude),
        `SmallestMagnitude` (smallest eigenvalues in magnitude) and
        `BothEnds` (eigenvalues at both ends of the spectrum;
        If `n_eigenvalues` is odd, compute one more from the high end
        than from the low end).

    * - `ncv`
      - `int`
      - min(2*n_eigenvalues + 2, N)
      - How many Lanczos vectors to generate at each iteration.

    * - `compute_eigenvectors`
      - `bool`
      - n/a
      - Request computation of eigenvectors in addition to the eigenvalues.

    * - `random_residual_vector`
      - `bool`
      - `true`
      - Use a randomly generated initial residual vector?

    * - `sigma`
      - `double`
      - `0`
      - Real eigenvalue shift :math:`\sigma` for spectral transformation modes.

    * - `tolerance`
      - `double`
      - Machine precision
      - Relative tolerance for Ritz value (eigenvalue) convergence.

    * - `max_iter`
      - `unsigned int`
      - `INT_MAX`
      - Maximum number of Lanczos update iterations allowed.

  .. note::

    In the spectral transformation modes, values of `eigenvalues_select` refer
    to the spectrum of the **transformed** problem, not the original one. For
    instance, `LargestMagnitude` used in the shift-and-invert mode will pick
    eigenvalues :math:`\lambda` closest to the shift :math:`\sigma`, because
    they correspond to the eigenvalues :math:`\mu = 1/(\lambda - \sigma)`
    that have the largest magnitude.

5. Optionally set the initial vector for Lanczos iteration if a better choice
   than a random vector is known. `random_residual_vector` parameter must
   be set to `false` for the changes made to the initial vector to take effect.

   A view of the residual vector is accessible via the method
   `residual_vector()` of the worker.

   .. code:: cpp

     // Set all components of the initial vector to 1.
     auto rv = worker.residual_vector();
     for(int i = 0; i < N; ++i) rv[i] = 1.0;

   One may also call `residual_vector()` later, after a diagonalization run
   has started, to retrieve the current residual vector.

6. Choose one of supported computational modes and perform diagonalization.
   In this part, user is supposed to call the `worker` object and pass the
   parameter structure as well as callable objects (*e.g.* lambda-functions)
   that represent action of operators :math:`\hat O` and :math:`\hat B` on
   a given vector. The supplied objects will be called to generate Lanczos
   vectors. Syntax and semantics of the C++ code vary between
   the computational modes and will be explained individually for each of
   them.

   - **Standard mode** (for standard eigenproblems, :math:`\hat M = \hat I`).

     .. code:: cpp

       using vector_view_t = worker_t::vector_view_t;
       using vector_const_view_t = worker_t::vector_const_view_t;

       auto Aop = [](vector_const_view_t from, vector_view_t to) {
         // Code implementing action of matrix A on vector 'from':
         // to = A * from
       };

       ar(Aop, params);

   - **Regular inverse mode** (for symmetric positive-definite :math:`\hat M`).

     In this mode, the transformed eigenproblem is defined by
     :math:`\hat O = \hat M^{-1} \hat A`, :math:`\hat B = \hat M` and
     :math:`\lambda = \mu`.

     .. code:: cpp

       using vector_view_t = worker_t::vector_view_t;
       using vector_const_view_t = worker_t::vector_const_view_t;

       auto op = [](vector_view_t from, vector_view_t to) {
         // Code implementing action of matrices M^{-1} and A according to
         // from = A * from;
         // to = invM * from;
         //
         // Note that unlike in the other computational modes, both 'from' and
         // 'to' must be updated!
       };
       auto Bop = [](vector_const_view_t from, vector_view_t to) {
         // Code implementing action of matrix M on vector 'from':
         // to = M * from
       };

       ar(op, Bop, worker_t::Inverse, params);

     Inverting a sparse matrix :math:`\hat M` will likely make it dense, which
     is usually undesirable from the storage standpoint. A more practical
     solution is to compute the sparse LU or Cholesky factorization of
     :math:`\hat M` once (outside of the lambda-function's body), and write
     the lambda-function so that it (1) sets `from = A * from` and (2) computes
     `to` as the solution of the linear system `M * to = from` using the
     precomputed factorization.

   - **Shift-and-Invert mode** (for symmetric positive semi-definite
     :math:`\hat M`).

     In this mode, the transformed eigenproblem is defined by
     :math:`\hat O = (\hat A -\sigma \hat M)^{-1} \hat M`,
     :math:`\hat B = \hat M` and :math:`\lambda = 1/\mu + \sigma`.
     The real spectral shift :math:`\sigma` must be set in the parameters
     structure, see table in step 4.

     .. code:: cpp

       using vector_view_t = worker_t::vector_view_t;
       using vector_const_view_t = worker_t::vector_const_view_t;

       auto op = [](vector_view_t from, vector_view_t to) {
         // Code implementing action of matrix (A - sigma*M)^{-1} * M on 'from'
         // to = (A - sigma*M)^{-1} * M * from;
       };
       auto Bop = [](vector_const_view_t from, vector_view_t to) {
         // Code implementing action of matrix M on vector 'from':
         // to = M * from
       };

       ar(op, Bop, worker_t::ShiftAndInvert, params);

     Inverting a sparse matrix :math:`\hat A - \sigma\hat M` will likely make it
     dense, which is usually undesirable from the storage standpoint. A more
     practical solution is to compute the sparse LU or Cholesky factorization of
     :math:`\hat A - \sigma\hat M` once (outside of the lambda-function's body),
     and write the lambda-function so that it (1) computes `M * from` and
     (2) computes `to` as the solution of the linear system
     `(A - \\sigma M) * to = M * from` using the precomputed factorization.

   - **Buckling mode** (for symmetric positive semi-definite
     :math:`\hat A` and symmetric indefinite :math:`\hat M`).

     In this mode, the transformed eigenproblem is defined by
     :math:`\hat O = (\hat A -\sigma \hat M)^{-1} \hat A`,
     :math:`\hat B = \hat A`, and :math:`\lambda = \sigma \frac{\mu}{\mu-1}`.
     The real spectral shift :math:`\sigma` must be set in the parameters
     structure, see table in step 4.

     .. code:: cpp

       using vector_view_t = worker_t::vector_view_t;
       using vector_const_view_t = worker_t::vector_const_view_t;

       auto op = [](vector_view_t from, vector_view_t to) {
         // Code implementing action of matrix
         // (A - sigma*M)^{-1} * A on 'from'
         // to = (A - sigma*M)^{-1} A * from;
       };
       auto Bop = [](vector_const_view_t from, vector_view_t to) {
         // Code implementing action of matrix A on vector 'from':
         // to = A * from
       };

       ar(op, Bop, worker_t::Buckling, params);

     Inverting a sparse matrix :math:`\hat A - \sigma\hat M` will likely make it
     dense, which is usually undesirable from the storage standpoint. A more
     practical solution is to compute the sparse LU or Cholesky factorization of
     :math:`\hat A - \sigma\hat M` once (outside of the lambda-function's body),
     and write the lambda-function so that it (1) computes
     `A * from` and (2) computes `to` as the solution of the linear
     system `(A - \\sigma M) * to = A * from` using the precomputed
     factorization.

   - **Cayley mode** (for symmetric positive semi-definite
     :math:`\hat M`).

     In this mode, the transformed eigenproblem is defined by
     :math:`\hat O = (\hat A -\sigma \hat M)^{-1} (\hat A + \sigma \hat M)`,
     :math:`\hat B = \hat M` and
     :math:`\lambda = \sigma\left(\frac{1+\mu}{1-\mu}\right)`.
     The real spectral shift :math:`\sigma` must be set in the parameters
     structure, see table in step 4.

     .. code:: cpp

       using vector_view_t = worker_t::vector_view_t;
       using vector_const_view_t = worker_t::vector_const_view_t;

       auto op = [](vector_view_t from, vector_view_t to) {
         // Code implementing action of matrix
         // (A - sigma*M)^{-1} * (A + sigma*M) on 'from'
         // to = (A - sigma*M)^{-1} * (A + sigma*M) * from;
       };
       auto Bop = [](vector_const_view_t from, vector_view_t to) {
         // Code implementing action of matrix M on vector 'from':
         // to = M * from
       };

       ar(op, Bop, worker_t::Cayley, params);

     Inverting a sparse matrix :math:`\hat A - \sigma\hat M` will likely make it
     dense, which is usually undesirable from the storage standpoint. A more
     practical solution is to compute the sparse LU or Cholesky factorization of
     :math:`\hat A - \sigma\hat M` once (outside of the lambda-function's body),
     and write the lambda-function so that it (1) computes
     `(A + \\sigma M) * from` and (2) computes `to` as the solution of the linear
     system `(A - \\sigma M) * to = (A + \\sigma M) * from` using the
     precomputed factorization.

   .. note::

     In most computational modes above, it is seemingly necessary to apply
     operator :math:`\hat B` to the same vector twice per generated Lanczos
     vector, once in functor `op` and once in `Bop`. It is actually possible to
     spare one of the applications. Calling `worker.Bx_available()` inside `op`
     will tell whether `Bop` has already been called at the current iteration,
     and `worker.Bx_vector()` will return a constant view of the application
     result :math:`\hat B \mathbf{x}`.

   The `from` and `to` views passed to the callable objects always expose one
   of three length-:math:`N` vectors stored inside the worker object. There is
   another, indirect way to access them.

   .. code:: cpp

     // Get index (0-2) of the current 'from' vector and request a view of it
     auto from_view = worker.workspace_vector(worker.from_vector_n());
     // Similar for the 'to' vector
     auto to_view = worker.workspace_vector(worker.to_vector_n());

   In advanced usage scenarios, the implicit restarting procedure can be
   customized via an extra argument of `worker`'s call operator.
   See :ref:`restarting` for more details.

   .. code:: cpp

     auto shifts_f = [](worker_t::real_vector_const_view_t ritz_values,
                        worker_t::real_vector_const_view_t ritz_bounds,
                        worker_t::real_vector_view_t shifts) {
                          // Compute shifts for the implicit restarting
                        };

     // Standard mode
     worker(op, params, shifts_f);
     // Other modes, e.g. Inverse
     worker(op, Bop, worker_t::Inverse, params, shifts_f);

   `worker_t::operator()` can throw two special exception types.

   - `maxiter_reached` - Maximum number of implicitly restarted Lanczos
     iterations has been reached.
   - `ncv_insufficient` - No shifts could be applied during a cycle of
     the Implicitly restarted Lanczos iteration. Consider increasing the number
     of Lanczos vectors generated at each iteration (`ncv` parameter).

   The rest of possible problems reported by ARPACK-NG result in generic
   `std::runtime_error` exceptions.

7. Request computed eigenvalues and eigenvectors (provided the
   `compute_eigenvectors` parameter has been enabled).

   .. code:: cpp

     auto lambda = worker.eigenvalues();
     auto vecs = worker.eigenvectors();

   The eigenvectors are columns of the real matrix view `vecs`.

8. Optionally request statistics about the completed run.

   .. code:: cpp

     // Print some computation statistics
     auto stats = worker.stats();

     std::cout << "Number of Arnoldi update iterations: " << stats.n_iter
               << std::endl;
     std::cout << "Number of 'converged' Ritz values: " << stats.n_converged
               << std::endl;
     std::cout << "Total number of O*x operations: " << stats.n_op_x_operations
               << std::endl;
     std::cout << "Total number of B*x operations: " << stats.n_b_x_operations
               << std::endl;
     std::cout << "Total number of steps of re-orthogonalization: "
               << stats.n_reorth_steps << std::endl;

   If a diagonalization run has ended prematurely (for example, when the maximum
   number of iterations has been reached), then it may still be possible to
   extract the first `stats.n_converged` eigenpairs.

