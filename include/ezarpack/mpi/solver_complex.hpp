/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2022 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/
/// @file ezarpack/mpi/solver_complex.hpp
/// @brief Specialization of `mpi::arpack_solver` class for the case of general
/// complex eigenproblems.
#pragma once

#include <algorithm>
#include <numeric>
#include <utility>
#include <vector>

#include "mpi_util.hpp"

namespace ezarpack {
namespace mpi {

/// @brief Main solver class wrapping the Implicitly Restarted Arnoldi
/// Method (IRAM) for general complex eigenproblems (MPI version).
///
/// This specialization of `arpack_solver` calls ARPACK-NG routines `pznaupd()`
/// and `pzneupd()` to compute approximations to a few eigenpairs of a complex
/// linear operator @f$ \hat O @f$ with respect to a semi-inner product defined
/// by a Hermitian positive semi-definite matrix @f$ \hat B @f$.
///
/// @note If both @f$ \hat O @f$ and @f$ \hat B@f$ are real and symmetric,
/// then @ref arpack_solver<Symmetric, Backend> should be used instead.
///
/// @tparam Backend Tag type specifying what *storage backend* (matrix/vector
/// algebra library) must be used by `arpack_solver`. The storage backend
/// determines types of internally stored data arrays and input/output view
/// objects returned by methods of the class.
template<typename Backend> class arpack_solver<Complex, Backend> {

  using storage = storage_traits<Backend>;

public:
  /// @name Backend-specific array and view types

  /// @{

  /// One-dimensional data array (vector) of real numbers.
  using real_vector_t = typename storage::real_vector_type;
  /// One-dimensional data array (vector) of complex numbers.
  using complex_vector_t = typename storage::complex_vector_type;
  /// Two-dimensional data array (matrix) of complex numbers.
  using complex_matrix_t = typename storage::complex_matrix_type;
  /// One-dimensional data array (vector) of integers.
  using int_vector_t = typename storage::int_vector_type;

  /// Partial view (slice) of a complex vector.
  using complex_vector_view_t = typename storage::complex_vector_view_type;
  /// Partial constant view (slice) of a complex vector.
  using complex_vector_const_view_t =
      typename storage::complex_vector_const_view_type;
  /// Partial constant view (slice) of a complex matrix.
  using complex_matrix_const_view_t =
      typename storage::complex_matrix_const_view_type;

  /// Storage-specific view type to expose complex input vectors
  /// @f$ \mathbf{x} @f$. An argument of this type is passed as input to
  /// callable objects representing linear operators @f$ \hat O @f$ and
  /// @f$ \hat B @f$.
  using vector_const_view_t = complex_vector_const_view_t;

  /// Storage-specific view type to expose complex output vectors
  /// @f$ \mathbf{y} @f$. An argument of this type receives output from
  /// callable objects representing linear operators @f$ \hat O @f$ and
  /// @f$ \hat B @f$.
  using vector_view_t = complex_vector_view_t;

  /// @}

private:
  MPI_Comm comm; // MPI communicator
  int comm_size; // MPI communicator size
  int comm_rank; // Rank within the MPI communicator

  int N; // Matrix size

  int block_start; // Index of the first element in the local block
  int block_size;  // Size of the local block

  const char* which;          // WHICH parameter
  int nev = 0;                // Number of eigenvalues
  double tol;                 // Relative tolerance for Ritz value convergence
  complex_vector_t resid;     // Residual vector
  complex_vector_t workd;     // Working space
  int ncv = 0;                // Number of Lanczos vectors to be generated
  complex_matrix_t v;         // Matrix with Arnoldi basis vectors
  int ldv = 0;                // Leading dimension of v
  complex_matrix_t z;         // Matrix with Ritz vectors
  int ldz = 0;                // Leading dimension of z
  complex_vector_t d;         // Ritz values (real and imaginary parts)
  int iparam[11];             // Various input/output parameters
  int ipntr[14];              // Starting locations in workd and workl
  int info = 0;               // !=0 to use resid, 0 otherwise
  int rvec;                   // RVEC parameter of pzneupd
  char howmny;                // HOWMNY parameter of pzneupd
  int_vector_t select;        // SELECT parameter of pzneupd
  bool Bx_available_ = false; // Has B*x already been computed?

public:
  /// Input parameters of the Implicitly Restarted Arnoldi Method (IRAM).
  struct params_t {

    /// Number of eigenvalues (Ritz values) to compute.
    unsigned int n_eigenvalues;

    /// Categories of eigenvalues to compute.
    enum eigenvalues_select_t {
      LargestMagnitude,  /**< Largest eigenvalues in magnitude. */
      SmallestMagnitude, /**< Smallest eigenvalues in magnitude. */
      LargestReal,       /**< Eigenvalues of largest real part. */
      SmallestReal,      /**< Eigenvalues of smallest real part. */
      LargestImag,       /**< Eigenvalues of largest imaginary part. */
      SmallestImag       /**< Eigenvalues of smallest imaginary part. */
    };

    /// Which of the eigenvalues to compute?
    eigenvalues_select_t eigenvalues_select;

    /// Number of Arnoldi vectors to be generated.
    /// `-1` stands for the default value `min(2*n_eigenvalues + 2, N)`,
    /// where `N` is the dimension of the problem.
    int ncv = -1;

    /// Kinds of vectors to compute.
    enum compute_vectors_t {
      None,  /**< Do not compute neither Ritz nor Schur vectors. */
      Schur, /**< Compute Schur vectors (orthogonal basis vectors of
                 the @ref n_eigenvalues -dimensional subspace). */
      Ritz   /**< Compute Ritz vectors (eigenvectors) in addition
                   to the orthogonal basis vectors (Schur vectors). */
    };

    /// Compute Ritz or Schur vectors?
    compute_vectors_t compute_vectors;

    /// Use a randomly generated initial residual vector?
    bool random_residual_vector = true;

    /// Eigenvalue shift @f$ \sigma @f$ used if a spectral transformation is
    /// employed.
    dcomplex sigma = 0;

    /// Relative tolerance for Ritz value convergence. The default setting is
    /// machine precision.
    double tolerance = 0;

    /// Maximum number of IRAM iterations allowed.
    unsigned int max_iter = INT_MAX;

    /// Constructs an IRAM parameter object with given
    /// @ref n_eigenvalues, @ref eigenvalues_select and
    /// @ref compute_vectors.
    /// The rest of the parameters are set to their defaults.
    params_t(unsigned int n_eigenvalues,
             eigenvalues_select_t eigenvalues_select,
             compute_vectors_t compute_vectors)
        : n_eigenvalues(n_eigenvalues),
          eigenvalues_select(eigenvalues_select),
          compute_vectors(compute_vectors) {}
  };

  /// Constructs a solver object and allocates internal data buffers to be
  /// used by ARPACK-NG.
  ///
  /// This constructor partitions N-dimensional vectors into blocks and
  /// distributes those blocks among all MPI ranks in a communicator in the most
  /// even way.
  /// @param N Dimension of the eigenproblem.
  /// @param comm MPI communicator.
  arpack_solver(unsigned int N, MPI_Comm const& comm)
      : comm(comm),
        comm_size(size(comm)),
        comm_rank(rank(comm)),
        N(N),
        block_start(compute_local_block_start(N, comm_size, comm_rank)),
        block_size(compute_local_block_size(N, comm_size, comm_rank)),
        resid(storage::make_complex_vector(block_size)),
        workd(storage::make_complex_vector(3 * block_size)),
        v(storage::make_complex_matrix(block_size, 0)),
        z(storage::make_complex_matrix(0, 0)),
        d(storage::make_complex_vector(nev + 1)),
        select(storage::make_int_vector(0)) {
    if(comm_size > N)
      throw ARPACK_SOLVER_ERROR("MPI communicator size cannot exceed dimension "
                                "of the eigenproblem (got " +
                                std::to_string(comm_size) + " vs " +
                                std::to_string(N) + ")");
    iparam[3] = 1;
  }

  /// Constructs a solver object and allocates internal data buffers to be
  /// used by ARPACK-NG.
  ///
  /// This constructor accepts a list of sizes of MPI rank-local vector blocks.
  /// @param block_sizes Sizes of MPI rank-local vector blocks, one element per
  /// MPI rank.
  /// @param comm MPI communicator.
  arpack_solver(std::vector<unsigned int> const& block_sizes,
                MPI_Comm const& comm)
      : comm(comm),
        comm_size(size(comm)),
        comm_rank(rank(comm)),
        N(std::accumulate(block_sizes.begin(), block_sizes.end(), 0)),
        block_start(0),
        block_size(0),
        resid(storage::make_complex_vector(block_size)),
        workd(storage::make_complex_vector(3 * block_size)),
        v(storage::make_complex_matrix(block_size, 0)),
        z(storage::make_complex_matrix(0, 0)),
        d(storage::make_complex_vector(nev + 1)),
        select(storage::make_int_vector(0)) {
    if(block_sizes.size() != comm_size)
      throw ARPACK_SOLVER_ERROR("Size of 'block_sizes' must coincide with MPI "
                                "communicator size (got " +
                                std::to_string(block_sizes.size()) + " vs " +
                                std::to_string(comm_size) + ")");
    block_size = block_sizes[comm_rank];
    if(block_size == 0)
      throw ARPACK_SOLVER_ERROR("Vector block " + std::to_string(comm_rank) +
                                " has zero size");
    block_start = std::accumulate(block_sizes.begin(),
                                  block_sizes.begin() + comm_rank, 0);

    iparam[3] = 1;
  }

  ~arpack_solver() {
    storage::destroy(resid);
    storage::destroy(workd);
    storage::destroy(v);
    storage::destroy(z);
    storage::destroy(d);
    storage::destroy(select);
  }

  arpack_solver(arpack_solver const&) = delete;
  // clang-format off
  arpack_solver(arpack_solver&&) noexcept(
    noexcept(int_vector_t(std::declval<int_vector_t>())) &&
    noexcept(complex_vector_t(std::declval<complex_vector_t>())) &&
    noexcept(complex_matrix_t(std::declval<complex_matrix_t>()))) = default;
  // clang-format on

private:
  /// @internal Prepare values of input parameters and resize containers.
  void prepare(params_t const& params) {

    // Check n_eigenvalues
    nev = params.n_eigenvalues;
    int nev_min = 1;
    int nev_max = N - 2;

    if(nev < nev_min || nev > nev_max)
      throw ARPACK_SOLVER_ERROR("n_eigenvalues must be within [" +
                                std::to_string(nev_min) + ";" +
                                std::to_string(nev_max) + "]");

    // Character codes for eigenvalues_select
    static const std::array<const char*, 6> wh = {"LM", "SM", "LR",
                                                  "SR", "LI", "SI"};
    which = wh[int(params.eigenvalues_select)];

    // Check ncv
    ncv = params.ncv;
    if(ncv == -1)
      ncv = std::min(2 * int(params.n_eigenvalues) + 2, N);
    else if(ncv <= int(params.n_eigenvalues + 1) || ncv > N)
      throw ARPACK_SOLVER_ERROR("ncv must be within ]" +
                                std::to_string(params.n_eigenvalues + 1) + ";" +
                                std::to_string(N) + "]");

    storage::resize(v, block_size, ncv);
    ldv = storage::get_col_spacing(v) >= 0 ? storage::get_col_spacing(v)
                                           : block_size;

    // Eigenvectors
    rvec = (params.compute_vectors != params_t::None);
    howmny = params.compute_vectors == params_t::Schur ? 'P' : 'A';
    storage::resize(select, ncv);

    // According to pzneupd() docs, 'z' is not referenced if howmny == 'P'.
    // In fact, however, passing a zero-size 'z' in the Schur vector mode
    // results in a SEGFAULT.
    if(rvec) {
      storage::resize(z, block_size, nev + 1);
      ldz = storage::get_col_spacing(z) >= 0 ? storage::get_col_spacing(z)
                                             : block_size;
    } else {
      storage::resize(z, 1, 1);
      ldz = 1;
    }

    // Tolerance
    tol = std::max(.0, params.tolerance);

    // Use random initial residual vector?
    info = !params.random_residual_vector;

    iparam[2] = int(params.max_iter); // Max number of iterations
    if(iparam[2] <= 0)
      throw ARPACK_SOLVER_ERROR(
          "Maximum number of Arnoldi update iterations must be positive");
  }

public:
  /// If this functor is used to provide shifts for implicit restart,
  /// then the default ARPACK-NG's shift strategy (Exact Shift Strategy)
  /// will be employed.
  /// @sa Paragraph 4.4.1 of ARPACK Users' Guide:
  /// Solution of Large Scale Eigenvalue Problems
  /// with Implicitly Restarted Arnoldi Methods (R. B. Lehoucq, D. C. Sorensen,
  /// C. Yang, SIAM, 1998),
  /// @arpack_manual_url#page=64
  struct exact_shifts_f {
    /// Trivial call operator. The actual shifts will be internally computed by
    /// ARPACK-NG.
    ///
    /// @param[in] ritz_values View of a complex vector with current
    /// @ref params_t::ncv Ritz values.
    /// @param[in] ritz_bounds View of a complex vector with current estimated
    /// error bounds of the Ritz values.
    /// @param[out] shifts Complex vector view to receive the computed shifts.
    void operator()(complex_vector_const_view_t ritz_values,
                    complex_vector_const_view_t ritz_bounds,
                    complex_vector_view_t shifts) {}
  };

  /// Solve a standard eigenproblem @f$ \hat A\mathbf{x} = \lambda\mathbf{x}@f$.
  ///
  /// @param a A callable object representing the linear operator
  /// @f$ \hat A @f$. It must take two arguments,
  /// @code
  /// a(vector_const_view_t in, vector_view_t out)
  /// @endcode
  /// `a` is expected to act on the vector view `in` and write the result into
  /// the vector view `out`, `out = a*in`. Both `in` and `out` are MPI
  /// rank-local blocks of their respective N-dimensional vectors. If `a` has
  /// matrix elements connecting blocks stored on different MPI ranks, it is
  /// user's duty to pass data between those ranks and to appropriately update
  /// `out` on all ranks. Given an instance `as` of the
  /// arpack_solver< Complex, Backend > class, `in` is also indirectly
  /// accessible as
  /// @code
  /// as.workspace_vector(as.in_vector_n())
  /// @endcode
  /// and `out` is accessible as
  /// @code
  /// as.workspace_vector(as.out_vector_n())
  /// @endcode
  ///
  /// @param params Set of input parameters for the Implicitly Restarted
  /// Arnoldi Method.
  ///
  /// @param shifts_f Functor that implements a shift selection strategy for the
  /// implicit restarts. For the expected signature of the functor, see
  /// @ref exact_shifts_f::operator()(). When this argument is omitted, the
  /// default @ref exact_shifts_f "\"Exact Shift Strategy\"" is used, which is
  /// the right choice in most cases.
  ///
  /// @throws ezarpack::ncv_insufficient No shifts could be applied during
  /// a cycle of the IRA iteration.
  /// @throws ezarpack::maxiter_reached Maximum number of IRA iterations has
  /// been reached. All possible eigenvalues of @f$ \hat O @f$ have been found.
  /// @throws std::runtime_error Invalid input parameters and other errors
  /// reported by ARPACK-NG routines `pznaupd()` and `pzneupd()`.
  template<typename A, typename ShiftsF = exact_shifts_f>
  void operator()(A&& a, params_t const& params, ShiftsF shifts_f = {}) {

    prepare(params);

    iparam[0] = (std::is_same<ShiftsF, exact_shifts_f>::value ? 1 : 0);
    iparam[6] = 1; // Mode 1, standard eigenproblem

    const int workl_size = 3 * ncv * ncv + 5 * ncv;
    complex_vector_t workl = storage::make_complex_vector(workl_size);
    real_vector_t rwork = storage::make_real_vector(ncv);

    rci_flag ido = Init;
    Bx_available_ = false;
    do {
      f77::paupd(comm, ido, "I", block_size, which, nev, tol,
                 storage::get_data_ptr(resid), ncv, storage::get_data_ptr(v),
                 ldv, iparam, ipntr, storage::get_data_ptr(workd),
                 storage::get_data_ptr(workl), workl_size,
                 storage::get_data_ptr(rwork), info);
      switch(ido) {
        case ApplyOpInit:
        case ApplyOp: {
          int in_pos = in_vector_n() * block_size;
          int out_pos = out_vector_n() * block_size;
          a(storage::make_vector_const_view(workd, in_pos, block_size),
            storage::make_vector_view(workd, out_pos, block_size));
        } break;
        case Shifts: {
          shifts_f(storage::make_vector_const_view(workl, ipntr[5] - 1, ncv),
                   storage::make_vector_const_view(workl, ipntr[7] - 1, ncv),
                   storage::make_vector_view(workl, ipntr[13] - 1, iparam[7]));
        } break;
        case Done: break;
        default: {
          storage::destroy(rwork);
          storage::destroy(workl);
          throw ARPACK_SOLVER_ERROR("Reverse communication interface error");
        }
      }
    } while(ido != Done);

    handle_paupd_error_codes(info, rwork, workl);

    storage::resize(d, nev + 1);
    complex_vector_t workev = storage::make_complex_vector(2 * ncv);

    f77::peupd(comm, rvec, &howmny, storage::get_data_ptr(select),
               storage::get_data_ptr(d), storage::get_data_ptr(z), ldz,
               params.sigma, storage::get_data_ptr(workev), "I", block_size,
               which, nev, tol, storage::get_data_ptr(resid), ncv,
               storage::get_data_ptr(v), ldv, iparam, ipntr,
               storage::get_data_ptr(workd), storage::get_data_ptr(workl),
               workl_size, storage::get_data_ptr(rwork), info);

    storage::destroy(workev);
    storage::destroy(rwork);
    storage::destroy(workl);

    handle_peupd_error_codes(info);
  }

  // clang-format off
  /// Computational modes for generalized eigenproblems.
  enum Mode : int {
    Inverse = 2,
    /**< Regular inverse mode.

    Solve a generalized eigenproblem
    @f$ \hat A\mathbf{x} = \lambda \hat M\mathbf{x} @f$ by reduction to
    the canonical form with @f$ \hat O = \hat M^{-1} \hat A @f$ and
    @f$ \hat B = \hat M @f$, where @f$ \hat M @f$ is Hermitian positive
    definite.
    */
    ShiftAndInvert = 3
    /**< Shift-and-Invert mode.

    Solve a generalized eigenproblem
    @f$ \hat A\mathbf{x} = \lambda \hat M\mathbf{x} @f$ by reduction to
    the canonical form with @f$ \hat O = (\hat A - \sigma\hat M)^{-1} \hat M @f$
    and @f$ \hat B = \hat M @f$, where @f$ \hat M @f$ is Hermitian positive
    semi-definite.
    */
  };
  // clang-format on

  /// Solve a generalized eigenproblem @f$ \hat A\mathbf{x} =
  /// \lambda\hat M\mathbf{x}@f$.
  ///
  /// Operators @f$ \hat O @f$ and @f$ \hat B @f$
  /// mentioned below are related to @f$ \hat A @f$ and @f$ \hat M @f$ via
  /// a @ref Mode -dependent spectral transformation.
  ///
  /// @param op A callable object representing the linear operator
  /// @f$ \hat O @f$. It must take two arguments,
  /// @code
  /// op(vector_const_view_t in, vector_view_t out)
  /// @endcode
  /// `op` is expected to act on the vector view `in` and write the result
  /// into the vector view `out`, `out = op*in`.
  /// Both `in` and `out` are MPI rank-local blocks of their respective
  /// N-dimensional vectors. If `op` has matrix elements connecting blocks
  /// stored on different MPI ranks, it is user's duty to pass data between
  /// those ranks and to appropriately update `out` on all ranks.
  /// Given an instance `as` of the arpack_solver< Complex, Backend > class,
  /// `in` is also indirectly accessible as
  /// @code
  /// as.workspace_vector(as.in_vector_n())
  /// @endcode
  /// and `out` is accessible as
  /// @code
  /// as.workspace_vector(as.out_vector_n())as
  /// @endcode
  ///
  /// @param b A callable object representing the linear operator
  /// @f$ \hat B @f$. It must take two arguments,
  /// @code
  /// b(vector_const_view_t in, vector_view_t out)
  /// @endcode
  /// `b` is expected to act on the vector view `in` and write the result into
  /// the vector view `out`, `out = b*in`. Both `in` and `out` are MPI
  /// rank-local blocks of their respective N-dimensional vectors. If `b` has
  /// matrix elements connecting blocks stored on different MPI ranks, it is
  /// user's duty to pass data between those ranks and to appropriately update
  /// `out` on all ranks.
  ///
  /// @param mode @ref Mode "Computational mode" to be used.
  /// @param params Set of input parameters for the Implicitly Restarted
  /// Arnoldi Method.
  ///
  /// @param shifts_f Functor that implements a shift selection strategy for the
  /// implicit restarts. For the expected signature of the functor, see
  /// @ref exact_shifts_f::operator()(). When this argument is omitted, the
  /// default @ref exact_shifts_f "\"Exact Shift Strategy\"" is used, which is
  /// the right choice in most cases.
  ///
  /// @throws ezarpack::ncv_insufficient No shifts could be applied during
  /// a cycle of the IRA iteration.
  /// @throws ezarpack::maxiter_reached Maximum number of IRA iterations has
  /// been reached. All possible eigenvalues of @f$ \hat O @f$ have been found.
  /// @throws std::runtime_error Invalid input parameters and other errors
  /// reported by ARPACK-NG routines `pznaupd()` and `pzneupd()`.
  template<typename OP, typename B, typename ShiftsF = exact_shifts_f>
  void operator()(OP&& op,
                  B&& b,
                  Mode mode,
                  params_t const& params,
                  ShiftsF shifts_f = {}) {

    prepare(params);

    iparam[0] = (std::is_same<ShiftsF, exact_shifts_f>::value ? 1 : 0);
    iparam[6] = mode; // Modes 2-3, generalized eigenproblem

    const int workl_size = 3 * ncv * ncv + 5 * ncv;
    complex_vector_t workl = storage::make_complex_vector(workl_size);
    real_vector_t rwork = storage::make_real_vector(ncv);

    rci_flag ido = Init;
    Bx_available_ = false;
    do {
      f77::paupd(comm, ido, "G", block_size, which, nev, tol,
                 storage::get_data_ptr(resid), ncv, storage::get_data_ptr(v),
                 ldv, iparam, ipntr, storage::get_data_ptr(workd),
                 storage::get_data_ptr(workl), workl_size,
                 storage::get_data_ptr(rwork), info);
      switch(ido) {
        case ApplyOpInit: {
          int in_pos = in_vector_n() * block_size;
          int out_pos = out_vector_n() * block_size;
          Bx_available_ = false;
          op(storage::make_vector_const_view(workd, in_pos, block_size),
             storage::make_vector_view(workd, out_pos, block_size));
        } break;
        case ApplyOp: {
          int in_pos = in_vector_n() * block_size;
          int out_pos = out_vector_n() * block_size;
          // B*x is available via Bx_vector()
          Bx_available_ = true;
          op(storage::make_vector_const_view(workd, in_pos, block_size),
             storage::make_vector_view(workd, out_pos, block_size));
        } break;
        case ApplyB: {
          int in_pos = in_vector_n() * block_size;
          int out_pos = out_vector_n() * block_size;
          b(storage::make_vector_const_view(workd, in_pos, block_size),
            storage::make_vector_view(workd, out_pos, block_size));
        } break;
        case Shifts: {
          shifts_f(storage::make_vector_const_view(workl, ipntr[5] - 1, ncv),
                   storage::make_vector_const_view(workl, ipntr[7] - 1, ncv),
                   storage::make_vector_view(workl, ipntr[13] - 1, iparam[7]));
        } break;
        case Done: break;
        default: {
          storage::destroy(rwork);
          storage::destroy(workl);
          throw ARPACK_SOLVER_ERROR("Reverse communication interface error");
        }
      }
    } while(ido != Done);

    handle_paupd_error_codes(info, rwork, workl);

    storage::resize(d, nev + 1);
    complex_vector_t workev = storage::make_complex_vector(2 * ncv);

    f77::peupd(comm, rvec, &howmny, storage::get_data_ptr(select),
               storage::get_data_ptr(d), storage::get_data_ptr(z), ldz,
               params.sigma, storage::get_data_ptr(workev), "G", block_size,
               which, nev, tol, storage::get_data_ptr(resid), ncv,
               storage::get_data_ptr(v), ldv, iparam, ipntr,
               storage::get_data_ptr(workd), storage::get_data_ptr(workl),
               workl_size, storage::get_data_ptr(rwork), info);

    storage::destroy(workev);
    storage::destroy(rwork);
    storage::destroy(workl);

    handle_peupd_error_codes(info);
  }

  /// Returns dimension of the eigenproblem.
  inline int dim() const { return N; }

  /// Returns MPI communicator used to construct this solver.
  inline MPI_Comm mpi_comm() const { return comm; }

  /// Returns the index of the first vector element in the local block.
  inline int local_block_start() const { return block_start; }
  /// Returns the size of the local block.
  inline int local_block_size() const { return block_size; }

  /// Returns the index of the workspace vector, which is currently expected to
  /// be acted upon by linear operator @f$ \hat O @f$ or @f$ \hat B @f$.
  inline int in_vector_n() const { return (ipntr[0] - 1) / block_size; }

  /// Returns the index of the workspace vector, which is currently expected to
  /// receive result from application of linear operator @f$ \hat O @f$ or
  /// @f$ \hat B @f$.
  inline int out_vector_n() const { return (ipntr[1] - 1) / block_size; }

  /// Returns a view of an MPI rank-local block of a vector within ARPACK-NG's
  /// workspace array.
  ///
  /// @param n Index of the workspace vector. Valid values are 0, 1 and 2.
  /// @throws std::runtime_error Invalid index value.
  complex_vector_view_t workspace_vector(int n) {
    if(n < 0 || n > 2)
      throw ARPACK_SOLVER_ERROR(
          "Valid indices of workspace vectors are 0, 1 and 2 (got " +
          std::to_string(n) + ")");
    return storage::make_vector_view(workd, n * block_size, block_size);
  }

  /// Number of "converged" Ritz values.
  unsigned int nconv() const { return iparam[4]; }

  /// Returns a constant view of a list of @ref nconv() converged
  /// eigenvalues.
  ///
  /// \note In the generalized eigenproblem @ref Mode "modes", this
  /// method always returns eigenvalues of the **original** problem.
  complex_vector_const_view_t eigenvalues() const {
    return storage::make_vector_const_view(d, 0, nconv());
  }

  /// Returns a matrix, whose @ref nconv() columns are MPI rank-local blocks of
  /// converged Ritz vectors (eigenvectors).
  ///
  /// @throws std::runtime_error Ritz vectors have not been computed in the
  /// last IRAM run.
  complex_matrix_const_view_t eigenvectors() const {
    if((!rvec) || (howmny != 'A'))
      throw ARPACK_SOLVER_ERROR(
          "Invalid method call: Ritz vectors have not been computed");
    return storage::make_matrix_const_view(z, block_size, nconv());
  }

  /// Returns a view of a matrix, whose @ref nconv() columns are
  /// MPI rank-local blocks of Schur basis vectors.
  /// @throws std::runtime_error Schur vectors have not been computed in the
  /// last IRAM run.
  complex_matrix_const_view_t schur_vectors() const {
    if(!rvec)
      throw ARPACK_SOLVER_ERROR(
          "Invalid method call: Schur vectors have not been computed");
    return storage::make_matrix_const_view(v, block_size, nconv());
  }

  /// Returns a view of the MPI rank-local block of the current residual vector.
  ///
  /// When @ref params_t::random_residual_vector is set to `false`, the view
  /// returned by this accessor can be used to set the initial residual vector.
  complex_vector_view_t residual_vector() {
    return storage::make_vector_view(resid);
  }

  /// Has @f$ \hat B\mathbf{x} @f$ already been computed at the current
  /// IRAM iteration?
  bool Bx_available() const { return Bx_available_; }

  /// Returns a constant view of the MPI rank-local block of the most recently
  /// computed vector @f$ \hat B\mathbf{x} @f$.
  complex_vector_const_view_t Bx_vector() const {
    int n = ipntr[2] - 1;
    return storage::make_vector_const_view(workd, n, N);
  }

  /// Statistics regarding a completed IRAM run.
  struct stats_t {
    /// Number of Arnoldi update iterations taken.
    unsigned int n_iter;
    /// Total number of @f$ \hat O \mathbf{x} @f$ operations.
    unsigned int n_op_x_operations;
    /// Total number of @f$ \hat B \mathbf{x} @f$ operations.
    unsigned int n_b_x_operations;
    /// Total number of steps of re-orthogonalization.
    unsigned int n_reorth_steps;
  };

  /// Returns computation statistics from the last IRAM run.
  stats_t stats() const {
    stats_t s;
    s.n_iter = iparam[2];
    s.n_op_x_operations = iparam[8];
    s.n_b_x_operations = iparam[9];
    s.n_reorth_steps = iparam[10];
    return s;
  }

private:
  /// @internal Translate pznaupd's INFO codes into C++ exceptions.
  ///
  /// @param error_code pznaupd's INFO code.
  void handle_paupd_error_codes(int error_code,
                                real_vector_t& rwork,
                                complex_vector_t& workl) {
    if(error_code == 0) return;

    storage::destroy(rwork);
    storage::destroy(workl);
    switch(error_code) {
      case 1: throw(maxiter_reached(iparam[2]));
      case 3: throw(ncv_insufficient(ncv));
      case -8:
        throw ARPACK_SOLVER_ERROR("Error in LAPACK eigenvalue calculation");
      case -9: throw ARPACK_SOLVER_ERROR("Starting vector is zero");
      case -9999:
        throw ARPACK_SOLVER_ERROR(
            "Could not build an Arnoldi factorization. "
            "The size of the current Arnoldi factorization is " +
            std::to_string(nconv()));
      default:
        throw ARPACK_SOLVER_ERROR("pznaupd failed with error code " +
                                  std::to_string(error_code));
    }
  }

  /// @internal Translate pzneupd's INFO codes into C++ exceptions.
  ///
  /// @param error_code pzneupd's INFO code.
  void handle_peupd_error_codes(int error_code) {
    switch(error_code) {
      case 0: return;
      case 1:
        throw ARPACK_SOLVER_ERROR(
            "The Schur form computed by LAPACK routine csheqr "
            "could not be reordered by LAPACK routine ztrsen");
      case -8:
        throw ARPACK_SOLVER_ERROR("Error in LAPACK eigenvalue calculation");
      case -9:
        throw ARPACK_SOLVER_ERROR(
            "Error in LAPACK eigenvectors calculation (ztrevc)");
      case -14:
        throw ARPACK_SOLVER_ERROR(
            "pznaupd did not find any eigenvalues to sufficient accuracy");
      default:
        throw ARPACK_SOLVER_ERROR("pzneupd failed with error code " +
                                  std::to_string(error_code));
    }
  }
};

} // namespace mpi
} // namespace ezarpack
