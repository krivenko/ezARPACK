.. _choice:

What algorithm variant is right for your problem?
=================================================

ARPACK-NG provides tools for solving a wide range of generalized eigenproblems

.. math::

  \hat A \mathbf{x} = \lambda \hat M \mathbf{x}

with large sparse matrices :math:`\hat A` and :math:`\hat M` (a standard
eigenproblem is one with :math:`\hat M = \hat I`). All supported classed of the
eigenproblems can be transformed to the canonical form

.. math::

  \hat O \mathbf{x} = \mu \hat B \mathbf{x},

where :math:`\hat B` is a Hermitian positive semi-definite matrix defining a
semi-inner product of vectors :math:`\mathbf{x}`.
Depending on what properties matrices :math:`\hat A` and :math:`\hat M` posses
-- being real or complex, symmetric, positive-definite, etc -- a different
algorithm and computational mode should be chosen for the optimal performance
and stability.

ezARPACK wraps ARPACK-NG's low-level subroutines in a form of three major
worker classes.

* :ref:`ezarpack::arpack_worker\<Symmetric, Backend\> <symmetric>`;
* :ref:`ezarpack::arpack_worker\<Asymmetric, Backend\> <asymmetric>`;
* :ref:`ezarpack::arpack_worker\<Complex, Backend\> <complex>`.

Here, they are listed in the order of generality of the eigenproblems they can
solve -- the symmetric variant is the most specialized and fastest one, whereas
the complex version is for the most general problems. Furthermore, each of the
variants supports a few computational modes. Picking the right combination of
worker and computational mode can be an overwhelming task for a non-expert.

In the table below, we give a classification of all supported eigenproblems with
recommendations on what worker variant/computational mode to choose.
In many cases, there are multiple acceptable choices, and the 'Notes' column
contains some more elaborate detail. Forms of the transformed matrices
:math:`\hat O` and :math:`\hat B` for each mode are also shown.

.. list-table::
  :header-rows: 1
  :align: left
  :widths: auto

  * - Type of :math:`\hat A`
    - Type of :math:`\hat M`
    - Worker variant
    - Computational mode
    - :math:`\hat O`
    - :math:`\hat B`
    - Notes

  * - Real symmetric
    - :math:`\hat I`
    - `Symmetric`
    - Standard
    - :math:`\hat A`
    - :math:`\hat I`
    - Optimal for finding extremal eigenvalues.

  * -
    -
    - `Symmetric`
    - `ShiftAndInvert`
    - :math:`(\hat A - \sigma)^{-1}`
    - :math:`\hat I`
    - Optimal for finding eigenvalues in the interior of the spectrum, clustered
      around the real shift :math:`\sigma`.

  * - Real symmetric
    - Real symmetric, positive-definite
    - `Symmetric`
    - `Invert`
    - :math:`\hat M^{-1} \hat A`
    - :math:`\hat M`
    - Optimal for finding extremal eigenvalues when :math:`\hat M` is
      well-conditioned.

  * -
    -
    - `Symmetric`
    - `ShiftAndInvert`
    - :math:`(\hat A-\sigma \hat M)^{-1} \hat M`
    - :math:`\hat M`
    - Optimal for finding eigenvalues in the interior of the spectrum, clustered
      around the real shift :math:`\sigma`.

  * - Real symmetric
    - Real symmetric, positive semi-definite
    - `Symmetric`
    - `ShiftAndInvert`
    - :math:`(\hat A - \sigma \hat M)^{-1} \hat M`
    - :math:`\hat M`
    - Optimal for finding eigenvalues clustered around the real shift
      :math:`\sigma`.

  * -
    -
    - `Symmetric`
    - `Cayley`
    - :math:`(\hat A - \sigma \hat M)^{-1} (\hat A + \sigma \hat M)`
    - :math:`\hat M`
    - Cayley-transformed eigenproblem. Another option for finding eigenvalues
      clustered around the real shift :math:`\sigma`. The transformation
      becomes ill-defined near :math:`\sigma = 0`.

  * - Real symmetric, positive semi-definite
    - Real symmetric
    - `Symmetric`
    - `Buckling`
    - :math:`(\hat A - \sigma \hat M)^{-1} \hat A`
    - :math:`\hat A`
    - Buckling-transformed eigenproblem. Optimal for finding eigenvalues
      clustered around the real shift :math:`\sigma`. The transformation
      becomes ill-defined near :math:`\sigma = 0`.

  * - General real
    - :math:`\hat I`
    - `Asymmetric`
    - Standard
    - :math:`\hat A`
    - :math:`\hat I`
    - Optimal for finding eigenvalues at the extreme points of the convex
      hull of the spectrum.

  * -
    -
    - `Asymmetric`
    - `ShiftAndInvertReal`
    - :math:`\Re [(\hat A - \sigma)^{-1}]`
    - :math:`\hat I`
    - Optimal for finding eigenvalues in the interior of the spectrum, clustered
      around the complex shift :math:`\sigma`. This mode must be chosen over
      `ShiftAndInvertImag` if :math:`\Im\sigma = 0`.

  * -
    -
    - `Asymmetric`
    - `ShiftAndInvertImag`
    - :math:`\Im [(\hat A - \sigma)^{-1}]`
    - :math:`\hat I`
    - Optimal for finding eigenvalues in the interior of the spectrum, clustered
      around the complex shift :math:`\sigma`. As :math:`\lambda` goes to
      infinity, the eigenvalues are damped more strongly in this mode than in
      `ShiftAndInvertReal`.

  * - General real
    - Real symmetric, positive-definite
    - `Asymmetric`
    - `Invert`
    - :math:`\hat M^{-1} \hat A`
    - :math:`\hat M`
    - Optimal for finding eigenvalues at the extreme points of the convex
      hull of the spectrum.

  * -
    -
    - `Asymmetric`
    - `ShiftAndInvertReal`
    - :math:`\Re [(\hat A - \sigma\hat M)^{-1} \hat M]`
    - :math:`\hat M`
    - Optimal for finding eigenvalues in the interior of the spectrum, clustered
      around the complex shift :math:`\sigma`. This mode must be chosen over
      `ShiftAndInvertImag` if :math:`\Im\sigma = 0`.

  * -
    -
    - `Asymmetric`
    - `ShiftAndInvertImag`
    - :math:`\Im [(\hat A - \sigma\hat M)^{-1} \hat M]`
    - :math:`\hat M`
    - Optimal for finding eigenvalues in the interior of the spectrum, clustered
      around the complex shift :math:`\sigma`. As :math:`\lambda` goes to
      infinity, the eigenvalues are damped more strongly in this mode than in
      `ShiftAndInvertImag`.

  * - General real
    - Real symmetric, positive semi-definite
    - `Asymmetric`
    - `ShiftAndInvertReal`
    - :math:`\Re [(\hat A - \sigma\hat M)^{-1} \hat M]`
    - :math:`\hat M`
    - Optimal for finding eigenvalues in the interior of the spectrum, clustered
      around the complex shift :math:`\sigma`. This mode must be chosen over
      `ShiftAndInvertImag` if :math:`\Im\sigma = 0`.

  * -
    -
    - `Asymmetric`
    - `ShiftAndInvertImag`
    - :math:`\Im [(\hat A - \sigma\hat M)^{-1} \hat M]`
    - :math:`\hat M`
    - Optimal for finding eigenvalues in the interior of the spectrum, clustered
      around the complex shift :math:`\sigma`. As :math:`\lambda` goes to
      infinity, the eigenvalues are damped more strongly in this mode than in
      `ShiftAndInvertImag`.

  * - General real
    - General real, invertible
    - `Asymmetric`
    - Standard
    - :math:`\hat M^{-1} \hat A`
    - :math:`\hat I`
    - **Not directly supported by ARPACK-NG**.

      One can manually form operator
      :math:`\hat O = \hat M^{-1} \hat A` and use the Asymmetric worker in the
      standard mode. Best used when :math:`\hat M` is well-conditioned and
      the eigenvalues of interest are at extreme points of the convex
      hull of the spectrum.

  * - General real
    - General real
    - `Asymmetric`
    - Standard
    - :math:`(\hat A - \sigma\hat M)^{-1} \hat M`
    - :math:`\hat I`
    - **Not directly supported by ARPACK-NG**.

      One can manually form operator
      :math:`\hat O = (\hat A - \sigma\hat M)^{-1}\hat M` and use the Asymmetric
      worker in the standard mode. Best used when :math:`\hat M` is nearly
      singular and/or for finding eigenvalues in the interior of the spectrum,
      clustered around the complex shift :math:`\sigma`. The eigenvalues
      :math:`\mu` computed by the worker must be manually back-transformed
      according to :math:`\lambda = \mu^{-1} + \sigma`.

  * - Complex
    - :math:`\hat I`
    - `Complex`
    - Standard
    - :math:`\hat A`
    - :math:`\hat I`
    - Optimal for finding eigenvalues at the extreme points of the convex
      hull of the spectrum.

  * -
    -
    - `Complex`
    - `ShiftAndInvert`
    - :math:`(\hat A - \sigma)^{-1}`
    - :math:`\hat I`
    - Optimal for finding eigenvalues in the interior of the spectrum, clustered
      around the complex shift :math:`\sigma`.

  * - Complex
    - Complex, Hermitian, positive-definite
    - `Complex`
    - `Invert`
    - :math:`\hat M^{-1} \hat A`
    - :math:`\hat M`
    - Optimal for finding eigenvalues at the extreme points of the convex
      hull of the spectrum when :math:`\hat M` is well-conditioned.

  * -
    -
    - `Complex`
    - `ShiftAndInvert`
    - :math:`(\hat A - \sigma \hat M)^{-1} \hat M`
    - :math:`\hat M`
    - Optimal for finding eigenvalues in the interior of the spectrum, clustered
      around the complex shift :math:`\sigma`.

  * - Complex
    - Complex, Hermitian, positive semi-definite
    - `Complex`
    - `ShiftAndInvert`
    - :math:`(\hat A - \sigma \hat M)^{-1} \hat M`
    - :math:`\hat M`
    - Optimal for finding eigenvalues in the interior of the spectrum, clustered
      around the complex shift :math:`\sigma`.

  * - Complex
    - Complex, invertible
    - `Complex`
    - Standard
    - :math:`\hat M^{-1} \hat A`
    - :math:`\hat I`
    - **Not directly supported by ARPACK-NG**.

      One can manually form operator
      :math:`\hat O = \hat M^{-1} \hat A` and use the Complex worker in the
      standard mode. Best used when :math:`\hat M` is well-conditioned and
      the eigenvalues of interest are at extreme points of the convex
      hull of the spectrum.

  * - Complex
    - General complex
    - `Complex`
    - Standard
    - :math:`(\hat A - \sigma\hat M)^{-1} \hat M`
    - :math:`\hat I`
    - **Not directly supported by ARPACK-NG**.

      One can manually form operator
      :math:`\hat O = (\hat A - \sigma\hat M)^{-1}\hat M` and use the Complex
      worker in the standard mode. Best used when :math:`\hat M` is nearly
      singular and/or for finding eigenvalues in the interior of the spectrum,
      clustered around the complex shift :math:`\sigma`. The eigenvalues
      :math:`\mu` computed by the worker must be manually back-transformed
      according to :math:`\lambda = \mu^{-1} + \sigma`.

Matrix :math:`\hat M` being well-conditioned means that it has a moderate
condition number :math:`||\hat M||_2\cdot||\hat M^{-1}||_2`.
The shift :math:`\sigma` used in various Shift-and-Invert modes has to be
provided by the user based on *a priori* knowledge about the spectrum. The
fastest convergence is achieved when it is close to the selected eigenvalues of
interest.

The table presented here is meant to give only some basic guidance. For a much
deeper overview of ARPACK-NG's capabilities you are referred to the definitive

  ARPACK Users' Guide: Solution of Large Scale Eigenvalue Problems
  with Implicitly Restarted Arnoldi Methods (R. B. Lehoucq, D. C. Sorensen,
  C. Yang, SIAM, 1998),
  https://www.caam.rice.edu/software/ARPACK/UG/

