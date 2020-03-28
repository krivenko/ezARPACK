/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2019 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/
#pragma once

namespace ezarpack {

/*******************************************************
 * ARPACK worker object: Case of real symmetric matrix A
 *******************************************************/
template<typename Backend> class arpack_worker<Symmetric, Backend> {

  using storage = storage_traits<Backend>;

  using real_vector_t = typename storage::real_vector_type;
  using real_matrix_t = typename storage::real_matrix_type;
  using int_vector_t = typename storage::int_vector_type;

  using real_vector_view_t = typename storage::real_vector_view_type;
  using real_vector_const_view_t =
      typename storage::real_vector_const_view_type;
  using real_matrix_const_view_t =
      typename storage::real_matrix_const_view_type;

  int N;                      // Matrix size
  const char* which;          // WHICH parameter
  int nev = 0;                // Number of eigenvalues
  double tol;                 // Relative tolerance for Ritz value convergence
  real_vector_t resid;        // Residual vector
  real_vector_t workd;        // Working space
  int ncv = 0;                // Number of Lanczos vectors to be generated
  real_matrix_t v;            // Matrix with Lanczos basis vectors
  real_vector_t d;            // Ritz values
  int iparam[11];             // Various input/output parameters
  int ipntr[11];              // Starting locations in workd and workl
  int info = 0;               // !=0 to use resid, 0 otherwise
  int rvec;                   // RVEC parameter of dseupd
  int_vector_t select;        // SELECT parameter of dseupd
  bool Bx_available_ = false; // Has B*x already been computed?

public:
  using vector_view_t = real_vector_view_t;
  using vector_const_view_t = real_vector_const_view_t;

  struct params_t {
    // Number of eigenvalues (Ritz values) to compute
    unsigned int n_eigenvalues;
    // Which of the Ritz values to compute
    enum {
      Largest,
      Smallest,
      LargestMagnitude,
      SmallestMagnitude,
      BothEnds
    } eigenvalues_select;

    // Expert option: number of Lanczos vectors to be generated
    // default: min(2*n_eigenvalues, N)
    int ncv = -1;

    // Compute Ritz vectors?
    bool compute_eigenvectors;

    // Use a randomly generated initial residual vector
    bool random_residual_vector = true;

    // Eigenvalue shift
    double sigma = 0;

    // Relative tolerance for Ritz value convergence
    double tolerance = 0;
    // Maximum number of Arnoldi update iterations
    unsigned int max_iter = INT_MAX;

    params_t(unsigned int n_eigenvalues,
             decltype(eigenvalues_select) ev_select,
             bool compute_eigenvectors)
        : n_eigenvalues(n_eigenvalues),
          eigenvalues_select(ev_select),
          compute_eigenvectors(compute_eigenvectors) {}
  };

  arpack_worker(unsigned int N)
      : N(N),
        resid(storage::make_real_vector(N)),
        workd(storage::make_real_vector(3 * N)),
        v(storage::make_real_matrix(N, 0)),
        d(storage::make_real_vector(nev)),
        select(storage::make_int_vector(0)) {
    iparam[3] = 1;
  }

  ~arpack_worker() {
    storage::destroy(resid);
    storage::destroy(workd);
    storage::destroy(v);
    storage::destroy(d);
    storage::destroy(select);
  }

  arpack_worker(arpack_worker const&) = delete;
  arpack_worker(arpack_worker&&) noexcept = delete;

  // Prepare values of input parameters
  void prepare(params_t const& params) {

    // Check n_eigenvalues
    nev = params.n_eigenvalues;
    int nev_min = params.eigenvalues_select == params_t::BothEnds ? 2 : 1;
    int nev_max = N - 1;

    if(nev < nev_min || nev > nev_max)
      throw ARPACK_WORKER_ERROR("n_eigenvalues must be within [" +
                                std::to_string(nev_min) + ";" +
                                std::to_string(nev_max) + "]");

    // Character codes for eigenvalues_select
    static const std::array<const char*, 5> wh = {"LA", "SA", "LM", "SM", "BE"};
    which = wh[int(params.eigenvalues_select)];

    // Check ncv
    ncv = params.ncv;
    if(ncv == -1)
      ncv = std::min(2 * int(params.n_eigenvalues) + 2, N);
    else if(ncv <= int(params.n_eigenvalues) || ncv > N)
      throw ARPACK_WORKER_ERROR("ncv must be within ]" +
                                std::to_string(params.n_eigenvalues) + ";" +
                                std::to_string(N) + "]");

    storage::resize(v, N, ncv);

    // Eigenvectors
    rvec = params.compute_eigenvectors;
    storage::resize(select, ncv);

    // Tolerance
    tol = std::max(.0, params.tolerance);

    // Use random initial residual vector?
    info = !params.random_residual_vector;

    iparam[2] = int(params.max_iter); // Max number of iterations
    if(iparam[2] <= 0)
      throw ARPACK_WORKER_ERROR(
          "Maximum number of Arnoldi update iterations must be positive");
  }

  struct trivial_shifts_f {
    void operator()(real_vector_view_t shifts) {}
  };

  /**********************************
   * Solve the standard eigenproblem
   * A*x = \lambda*x
   **********************************/

  // a: callable taking 2 arguments
  // a(vector_const_view_t from, vector_view_t to)
  // 'a' is expected to act on 'from' and write the result to 'to': to = A*from
  //
  // 'from' is also indirectly available as
  // this->workspace_vector(this->from_vector_n())
  // 'to' is also indirectly available as
  // this->workspace_vector(this->to_vector_n())
  //
  // shifts_f: callable taking one argument
  // shifts_f(vector_view_t shifts)
  // 'shifts_f' is expected to place the shifts for implicit restart
  // into 'shifts'
  template<typename A, typename ShiftsF = trivial_shifts_f>
  void operator()(A&& a, params_t const& params, ShiftsF shifts_f = {}) {

    prepare(params);

    iparam[0] = (std::is_same<ShiftsF, trivial_shifts_f>::value ? 1 : 0);
    iparam[6] = 1; // Mode 1, standard eigenproblem

    const int workl_size = ncv * ncv + 8 * ncv;
    real_vector_t workl = storage::make_real_vector(workl_size);

    rci_flag ido = Init;
    Bx_available_ = false;
    do {
      f77::aupd<true>(ido, "I", N, which, nev, tol,
                      storage::get_data_ptr(resid), ncv,
                      storage::get_data_ptr(v), N, iparam, ipntr,
                      storage::get_data_ptr(workd),
                      storage::get_data_ptr(workl), workl_size, info);
      switch(ido) {
        case ApplyOpInit:
        case ApplyOp: {
          int from_pos = from_vector_n() * N;
          int to_pos = to_vector_n() * N;
          a(storage::make_vector_const_view(workd, from_pos, N),
            storage::make_vector_view(workd, to_pos, N));
        } break;
        case Shifts:
          shifts_f(storage::make_vector_view(workl, ipntr[10] - 1, iparam[7]));
          break;
        case Done: break;
        default: {
          storage::destroy(workl);
          throw ARPACK_WORKER_ERROR("Reverse communication interface error");
        }
      }
    } while(ido != Done);

    handle_aupd_error_codes(info, workl);

    storage::resize(d, nev);

    f77::eupd(rvec, "A", storage::get_data_ptr(select),
              storage::get_data_ptr(d), storage::get_data_ptr(v), N,
              params.sigma, "I", N, which, nev, tol,
              storage::get_data_ptr(resid), ncv, storage::get_data_ptr(v), N,
              iparam, ipntr, storage::get_data_ptr(workd),
              storage::get_data_ptr(workl), workl_size, info);

    storage::destroy(workl);

    handle_eupd_error_codes(info);
  }

  /***********************************
   * Solve the generalized eigenproblem
   * A*x = \lambda*M*x
   ************************************/

  // clang-format off
  enum Mode : int {Invert = 2,           // OP = inv[M]*A and B = M
                   ShiftAndInvert = 3,   // OP = (inv[A - sigma*M])*M and B = M
                   Buckling = 4,         // OP = (inv[A - sigma*M])*A and B = A
                   Cayley = 5};          // OP = inv[A - sigma*M]*[A + sigma*M]
                                         // and B = M
  // clang-format on

  // op: callable taking 2 arguments
  // op(vector_view_t from, vector_view_t to)
  // In all modes except for 'Invert', 'op' is expected to act on 'from'
  // and write the result to 'to': to = OP*from.
  // In the 'Invert' mode, however, 'op' must do the following:
  //  from = A * from
  //  to = inv[M] * from
  //
  // b: callable taking 2 arguments
  // b(vector_const_view_t from, vector_view_t to)
  // 'b' is expected to act on 'from' and write the result to 'to': to = B*from
  //
  // 'from' is also indirectly available as
  // this->workspace_vector(this->vector_from_n())
  // 'to' is also indirectly available as
  // this->workspace_vector(this->vector_to_n())
  //
  // this->Bx_available() indicates whether B*x has already been computed
  // and available as this->Bx_vector()
  //
  // shifts_f: callable taking one argument
  // shifts_f(vector_view_t shifts)
  // 'shifts_f' is expected to place the shifts for implicit restart
  // into 'shifts'
  template<typename OP, typename B, typename ShiftsF = trivial_shifts_f>
  void operator()(OP&& op,
                  B&& b,
                  Mode mode,
                  params_t const& params,
                  ShiftsF shifts_f = {}) {

    prepare(params);

    iparam[0] = (std::is_same<ShiftsF, trivial_shifts_f>::value ? 1 : 0);
    iparam[6] = mode; // Modes 2-5, generalized eigenproblem

    const int workl_size = ncv * ncv + 8 * ncv;
    real_vector_t workl = storage::make_real_vector(workl_size);

    rci_flag ido = Init;
    Bx_available_ = false;
    do {
      f77::aupd<true>(ido, "G", N, which, nev, tol,
                      storage::get_data_ptr(resid), ncv,
                      storage::get_data_ptr(v), N, iparam, ipntr,
                      storage::get_data_ptr(workd),
                      storage::get_data_ptr(workl), workl_size, info);
      switch(ido) {
        case ApplyOpInit: {
          int from_pos = from_vector_n() * N;
          int to_pos = to_vector_n() * N;
          Bx_available_ = false;
          op(storage::make_vector_view(workd, from_pos, N),
             storage::make_vector_view(workd, to_pos, N));
        } break;
        case ApplyOp: {
          int from_pos = from_vector_n() * N;
          int to_pos = to_vector_n() * N;
          // B*x is available via Bx_vector()
          Bx_available_ = true;
          op(storage::make_vector_view(workd, from_pos, N),
             storage::make_vector_view(workd, to_pos, N));
        } break;
        case ApplyB: {
          int from_pos = from_vector_n() * N;
          int to_pos = to_vector_n() * N;
          b(storage::make_vector_const_view(workd, from_pos, N),
            storage::make_vector_view(workd, to_pos, N));
        } break;
        case Shifts:
          shifts_f(storage::make_vector_view(workl, ipntr[10] - 1, iparam[7]));
          break;
        case Done: break;
        default: {
          storage::destroy(workl);
          throw ARPACK_WORKER_ERROR("Reverse communication interface error");
        }
      }
    } while(ido != Done);

    handle_aupd_error_codes(info, workl);

    storage::resize(d, nev);
    double sigma = (mode != Invert) ? params.sigma : 0;

    f77::eupd(rvec, "A", storage::get_data_ptr(select),
              storage::get_data_ptr(d), storage::get_data_ptr(v), N, sigma, "G",
              N, which, nev, tol, storage::get_data_ptr(resid), ncv,
              storage::get_data_ptr(v), N, iparam, ipntr,
              storage::get_data_ptr(workd), storage::get_data_ptr(workl),
              workl_size, info);

    storage::destroy(workl);

    handle_eupd_error_codes(info);
  }

  // Index of workspace vector, which is expected to be acted on
  inline int from_vector_n() const { return (ipntr[0] - 1) / N; }

  // Index of workspace vector, which will receive result of the operator action
  inline int to_vector_n() const { return (ipntr[1] - 1) / N; }

  // Get view of a workspace vector
  real_vector_view_t workspace_vector(int n) const {
    if(n < 0 || n > 2)
      throw ARPACK_WORKER_ERROR(
          "Valid indices of workspace vectors are 0, 1 and 2 (got " +
          std::to_string(n) + ")");
    return storage::make_vector_view(workd, n * N, N);
  }

  // Access eigenvalues
  real_vector_const_view_t eigenvalues() const {
    return storage::make_vector_const_view(d, 0, iparam[4]);
  }

  // Access eigenvectors
  real_matrix_const_view_t eigenvectors() const {
    return storage::make_matrix_const_view(v, N, iparam[4]);
  }

  // Access residual vector
  real_vector_view_t residual_vector() {
    return storage::make_vector_view(resid);
  }

  // Has B*x already been computed?
  bool Bx_available() const { return Bx_available_; }

  // Previously computed vector B*x
  real_vector_const_view_t Bx_vector() const {
    unsigned int n = ipntr[2] - 1;
    return storage::make_vector_const_view(workd, n, N);
  }

  struct stats_t {
    // Number of Arnoldi update iterations
    unsigned int n_iter;
    // Number of "converged" Ritz values
    unsigned int n_converged;
    // Total number of OP*x operations
    unsigned int n_op_x_operations;
    // Total number of B*x operations
    unsigned int n_b_x_operations;
    // Total number of steps of re-orthogonalization
    unsigned int n_reorth_steps;
  };

  // Return computation statistics
  stats_t stats() const {
    stats_t s;
    s.n_iter = iparam[2];
    s.n_converged = iparam[4];
    s.n_op_x_operations = iparam[8];
    s.n_b_x_operations = iparam[9];
    s.n_reorth_steps = iparam[10];
    return s;
  }

private:
  // Error handling

  void handle_aupd_error_codes(int info, real_vector_t& workl) {
    if(info == 0) return;

    storage::destroy(workl);
    switch(info) {
      case 1: throw(maxiter_reached(iparam[2]));
      case 3: throw(ncv_insufficient(ncv));
      case -8:
        throw ARPACK_WORKER_ERROR(
            "Error in LAPACK tridiagonal eigenvalue calculation (dsteqr)");
      case -9: throw ARPACK_WORKER_ERROR("Starting vector is zero");
      case -13:
        throw ARPACK_WORKER_ERROR("n_eigenvalues = 1 is incompatible with "
                                  "eigenvalues_select = BothEnds");
      case -9999:
        throw ARPACK_WORKER_ERROR(
            "Could not build an Arnoldi factorization. "
            "The size of the current Arnoldi factorization is " +
            std::to_string(iparam[4]));
      default:
        throw ARPACK_WORKER_ERROR("dsaupd failed with error code " +
                                  std::to_string(info));
    }
  }
  void handle_eupd_error_codes(int info) {
    switch(info) {
      case 0: return;
      case -8:
        throw ARPACK_WORKER_ERROR(
            "Error in LAPACK tridiagonal eigenvalue calculation (dsteqr)");
      case -9: throw ARPACK_WORKER_ERROR("Starting vector is zero");
      case -12:
        throw ARPACK_WORKER_ERROR("n_eigenvalues = 1 is incompatible with "
                                  "eigenvalues_select = BothEnds");
      case -14:
        throw ARPACK_WORKER_ERROR(
            "dsaupd did not find any eigenvalues to sufficient accuracy");
      default:
        throw ARPACK_WORKER_ERROR("dseupd failed with error code " +
                                  std::to_string(info));
    }
  }
};

} // namespace ezarpack
