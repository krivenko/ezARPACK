/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2018 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/
#pragma once

namespace ezarpack {

/********************************************************
 * ARPACK worker object: Case of real asymmetric matrix A
 ********************************************************/
template<typename Backend> class arpack_worker<Asymmetric, Backend> {

  using storage = storage_traits<Backend>;

  using real_vector_t = typename storage::real_vector_type;
  using real_matrix_t = typename storage::real_matrix_type;
  using complex_vector_t = typename storage::complex_vector_type;
  using complex_matrix_t = typename storage::complex_matrix_type;
  using int_vector_t = typename storage::int_vector_type;

  using real_vector_view_t = typename storage::real_vector_view_type;
  using real_vector_const_view_t = typename storage::real_vector_const_view_type;

  int N;                       // Matrix size
  const char* which;           // WHICH parameter
  int nev = 0;                 // Number of eigenvalues
  int tol;                     // Relative tolerance for Ritz value convergence
  real_vector_t resid;         // Residual vector
  int ncv = 0;                 // Number of Lanczos vectors to be generated
  real_matrix_t v;             // Matrix with Arnoldi basis vectors
  real_matrix_t z;             // Matrix with Ritz vectors
  real_vector_t dr, di;        // Ritz values (real and imaginary parts)
  int iparam[11];              // Various input/output parameters
  int ipntr[14];               // Pointer to mark the starting locations in the workd and workl
  real_vector_t workd;         // Working space
  int info = 0;                // !=0 to use resid, 0 otherwise
  int rvec = false;            // RVEC parameter of dseupd
  char howmny;                 // HOWMNY parameter of dseupd
  int_vector_t select;         // SELECT parameter of dseupd
  double sigmar, sigmai;       // SIGMAR and SIGMAI parameters of dseupd
  bool Bx_available_ = false;  // Has B*x already been computed?

public:

  struct params_t {
    // Number of eigenvalues (Ritz values) to compute
    unsigned int n_eigenvalues;
    // Which of the Ritz values to compute
    enum {LargestMagnitude, SmallestMagnitude,
          LargestReal, SmallestReal, LargestImag, SmallestImag} eigenvalues_select;

    // Expert option: number of Lanczos vectors to be generated
    // default: min(2*n_eigenvalues+1, N)
    int ncv = -1;

    // Compute Ritz/Schur vectors?
    // None: do not compute anything;
    // Ritz: compute eigenvectors of A;
    // Schur: compute orthogonal basis vectors of the n_eigenvalues-dimensional subspace.
    enum {None, Ritz, Schur} compute_vectors;

    // Use a randomly generated initial residual vector
    bool random_residual_vector = true;

    // Eigenvalue shift
    dcomplex sigma = 0;

    // Relative tolerance for Ritz value convergence
    double tolerance = 0;
    // Maximum number of Arnoldi update iterations
    unsigned int max_iter = INT_MAX;

    params_t(unsigned int n_eigenvalues,
             decltype(eigenvalues_select) ev_select,
             decltype(compute_vectors) compute_evec)
      : n_eigenvalues(n_eigenvalues),
        eigenvalues_select(ev_select),
        compute_vectors(compute_evec)
    {}
  };

  arpack_worker(unsigned int N)
    : N(N),
      resid(storage::make_real_vector(N)),
      workd(storage::make_real_vector(3*N)),
      v(storage::make_real_matrix(N, 0)),
      z(storage::make_real_matrix(0, 0)),
      dr(storage::make_real_vector(nev+1)),
      di(storage::make_real_vector(nev+1)),
      select(storage::make_int_vector(0)) {
    iparam[3] = 1;
  }

  ~arpack_worker() {
    storage::destroy(resid, N);
    storage::destroy(workd, 3*N);
    storage::destroy(v, N, ncv);
    if(rvec)
      storage::destroy(z, N, nev+1);
    else
      storage::destroy(z, 0, 0);
    storage::destroy(dr, nev+1);
    storage::destroy(di, nev+1);
    storage::destroy(select, ncv);
  }

  arpack_worker(arpack_worker const&) = delete;
  arpack_worker(arpack_worker &&) noexcept = delete;

  // Prepare values of input parameters
  void prepare(params_t const& params) {

    // Check n_eigenvalues
    nev = params.n_eigenvalues;
    unsigned int nev_min = 1;
    unsigned int nev_max = N-2;

    if(nev < nev_min || nev > nev_max)
      throw std::runtime_error("arpack_worker: n_eigenvalues must be within ["
                               + std::to_string(nev_min) + ";"
                               + std::to_string(nev_max) + "]");

    // Character codes for eigenvalues_select
    static const std::array<const char*,6> wh = {"LM","SM","LR","SR","LI","SI"};
    which = wh[int(params.eigenvalues_select)];

    // Check ncv
    ncv = params.ncv;
    if(ncv == -1) ncv = std::min(2*int(params.n_eigenvalues)+2, N);
    else if(ncv <= params.n_eigenvalues+2 || ncv > N)
      throw std::runtime_error("arpack_worker: ncv must be within ]"
                               + std::to_string(params.n_eigenvalues+2) + ";"
                               + std::to_string(N) + "]");

    storage::resize(v, N, ncv);

    // Eigenvectors
    rvec = (params.compute_vectors != params_t::None);
    howmny = params.compute_vectors == params_t::Schur ? 'P' : 'A';
    storage::resize(select, ncv);
    if(rvec)
      storage::resize(z, N, nev+1);
    else
      storage::resize(z, 0, 0);

    // Tolerance
    tol = std::max(.0, params.tolerance);
    if(params.tolerance < 0) {
      std::cerr << "arpack_worker: negative tolerance " << params.tolerance
                << " is interpreted as machine epsilon." << std::endl;
    }

    // Use random initial residual vector?
    info = !params.random_residual_vector;

    iparam[2] = int(params.max_iter); // Max number of iterations
    if(iparam[2] <= 0)
      throw std::runtime_error("arpack_worker: maximum number of Arnoldi update iterations must be positive");
  }

  struct trivial_shifts_f {
    void operator()(real_vector_view_t shifts_re, real_vector_view_t shifts_im){}
  };

  /**********************************
  * Solve the standard eigenproblem
  * A*x = \lambda*x
  **********************************/

  // a: callable taking 2 arguments
  // a(real_vector_const_view_t from, real_vector_view_t to)
  // 'a' is expected to act on 'from' and write the result to 'to': to = A*from
  //
  // 'from' is also indirectly available as this->workspace_vector(this->vector_from_n())
  // 'to' is also indirectly available as this->workspace_vector(this->vector_to_n())
  //
  // shifts_f: callable taking two arguments
  // shifts_f(real_vector_view_t shifts_re, real_vector_view_t shifts_im)
  // 'shifts_f' is expected to place real and imaginary parts of the shifts for implicit restart
  // into 'shifts_re' and 'shifts_im' respectively.
  template<typename A, typename ShiftsF = trivial_shifts_f>
  void operator()(A a, params_t const& params, ShiftsF shifts_f = {}) {

    prepare(params);

    iparam[0] = (std::is_same<ShiftsF,trivial_shifts_f>::value ? 1 : 0);
    iparam[6] = 1; // Mode 1, standard eigenproblem

    const int workl_size = 3*ncv*ncv + 6*ncv;
    real_vector_t workl = storage::make_real_vector(workl_size);

    rci_flag ido = Init;
    Bx_available_ = false;
    do {
      f77::aupd<false>((int&)ido, "I", N, which,
                       nev, tol, storage::get_data_ptr(resid), ncv,
                       storage::get_data_ptr(v), N,
                       iparam, ipntr, storage::get_data_ptr(workd),
                       storage::get_data_ptr(workl), workl_size,
                       info);
      switch(ido) {
        case ApplyOpInit:
        case ApplyOp: {
          int from_pos = from_vector_n() * N;
          int to_pos = to_vector_n() * N;
          a(storage::make_vector_const_view(workd, from_pos, N),
            storage::make_vector_view(workd, to_pos, N)
          );
        }
        break;
        case Shifts: {
          int np = iparam[7];
          int shifts_re_n = ipntr[13]-1;
          int shifts_im_n = ipntr[13]-1+np;
          shifts_f(storage::make_vector_view(workl, shifts_re_n, np),
                   storage::make_vector_view(workl, shifts_im_n, np)
          );
        }
        break;
        case Done: break;
        default: {
          storage::destroy(workl, workl_size);
          throw std::runtime_error("arpack_worker: reverse communication interface error");
        }
      }
    } while(ido != Done);

    switch(info) {
      case 0: break;
      case 1: {
        storage::destroy(workl, workl_size);
        throw(maxiter_reached(iparam[2]));
      }
      case 3: {
        storage::destroy(workl, workl_size);
        throw(ncv_insufficient(ncv));
      }
      default: {
        storage::destroy(workl, workl_size);
        throw std::runtime_error("arpack_worker: dnaupd failed with error code " + std::to_string(info));
      }
    }

    storage::resize(dr, nev+1);
    storage::resize(di, nev+1);
    double sigmar, sigmai;
    real_vector_t workev = storage::make_real_vector(3*ncv);

    f77::eupd(rvec, &howmny, storage::get_data_ptr(select),
              storage::get_data_ptr(dr), storage::get_data_ptr(di),
              storage::get_data_ptr(z), (rvec ? N  : 0),
              sigmar, sigmai, storage::get_data_ptr(workev),
              "I", N, which, nev, tol, storage::get_data_ptr(resid), ncv,
              storage::get_data_ptr(v), N,
              iparam, ipntr, storage::get_data_ptr(workd),
              storage::get_data_ptr(workl), workl_size,
              info);

    storage::destroy(workl, workl_size);

    if(info) throw std::runtime_error("arpack_worker: dneupd failed with error code " + std::to_string(info));
  }

  /***********************************
  * Solve the generalized eigenproblem
  * A*x = \lambda*M*x
  ************************************/

  enum Mode : int {Invert = 2,               // OP = inv[M]*A and B = M
                   ShiftAndInvertReal = 3,   // OP = Real_Part{ inv[A - sigma*M]*M } and B = M
                   ShiftAndInvertImag = 4};  // OP = Imaginary_Part{ inv[A - sigma*M]*M } and B = M

  // op: callable taking 2 arguments
  // op(real_vector_const_view_t from, real_vector_view_t to)
  // 'op' is expected to act on 'from' and write the result to 'to': to = OP*from
  //
  // b: callable taking 2 arguments
  // b(real_vector_const_view_t from, real_vector_view_t to)
  // 'b' is expected to act on 'from' and write the result to 'to': to = B*from
  //
  // 'from' is also indirectly available as this->workspace_vector(this->from_vector_n())
  // 'to' is also indirectly available as this->workspace_vector(this->to_vector_n())
  //
  // this->Bx_available() indicates whether B*x has already been computed and available as this->Bx_vector()
  //
  // shifts_f: callable taking two arguments
  // shifts_f(real_vector_view_t shifts_re, real_vector_view_t shifts_im)
  // 'shifts_f' is expected to place real and imaginary parts of the shifts for implicit restart
  // into 'shifts_re' and 'shifts_im' respectively.
  template<typename OP, typename B, typename ShiftsF = trivial_shifts_f>
  void operator()(OP op, B b, Mode mode, params_t const& params, ShiftsF shifts_f = {}) {

    // https://github.com/opencollab/arpack-ng/issues/3
    // http://forge.scilab.org/index.php/p/arpack-ng/issues/1315/
    if(mode == ShiftAndInvertReal || mode == ShiftAndInvertImag)
      throw std::runtime_error("arpack_worker: ShiftAndInvertReal and ShiftAndInvertImag modes "
                               "are disabled due to potential problems with dneupd");

    prepare(params);

    iparam[0] = (std::is_same<ShiftsF,trivial_shifts_f>::value ? 1 : 0);
    iparam[6] = mode; // Modes 2-5, generalized eigenproblem

    const int workl_size = 3*ncv*ncv + 6*ncv;
    real_vector_t workl = storage::make_real_vector(workl_size);

    rci_flag ido = Init;
    Bx_available_ = false;
    do {
      f77::aupd<false>((int&)ido, "G", N, which,
                       nev, tol, storage::get_data_ptr(resid), ncv,
                       storage::get_data_ptr(v), N,
                       iparam, ipntr, storage::get_data_ptr(workd),
                       storage::get_data_ptr(workl), workl_size,
                       info);
      switch(ido) {
        case ApplyOpInit: {
          int from_pos = from_vector_n() * N;
          int to_pos = to_vector_n() * N;
          Bx_available_ = false;
          op(storage::make_vector_const_view(workd, from_pos, N),
             storage::make_vector_view(workd, to_pos, N)
          );
        }
        break;
        case ApplyOp: {
          int from_pos = from_vector_n() * N;
          int to_pos = to_vector_n() * N;
          // B*x is available via Bx_vector()
          Bx_available_ = true;
          op(storage::make_vector_const_view(workd, from_pos, N),
             storage::make_vector_view(workd, to_pos, N)
          );
        }
        break;
        case ApplyB: {
          int from_pos = from_vector_n() * N;
          int to_pos = to_vector_n() * N;
          b(storage::make_vector_const_view(workd, from_pos, N),
            storage::make_vector_view(workd, to_pos, N)
          );
        }
        break;
        case Shifts: {
          int np = iparam[7];
          int shifts_re_n = ipntr[13]-1;
          int shifts_im_n = ipntr[13]-1+np;
          shifts_f(storage::make_vector_view(workl, shifts_re_n, np),
                   storage::make_vector_view(workl, shifts_im_n, np)
          );
        }
        break;
        case Done: break;
        default: {
          storage::destroy(workl, workl_size);
          throw std::runtime_error("arpack_worker: reverse communication interface error");
        }
      }
    } while(ido != Done);

    switch(info) {
      case 0: break;
      case 1: {
        storage::destroy(workl, workl_size);
        throw(maxiter_reached(iparam[2]));
      }
      case 3: {
        storage::destroy(workl, workl_size);
        throw(ncv_insufficient(ncv));
      }
      default: {
        storage::destroy(workl, workl_size);
        throw std::runtime_error("arpack_worker: dnaupd failed with error code " + std::to_string(info));
      }
    }

    storage::resize(dr, nev+1);
    storage::resize(di, nev+1);
    if(mode != Invert) {
      sigmar = params.sigma.real();
      sigmai = params.sigma.imag();
    }
    real_vector_t workev = storage::make_real_vector(3*ncv);

    f77::eupd(rvec, &howmny, storage::get_data_ptr(select),
              storage::get_data_ptr(dr), storage::get_data_ptr(di),
              storage::get_data_ptr(z), (rvec ? N  : 0),
              sigmar, sigmai, storage::get_data_ptr(workev),
              "G", N, which, nev, tol, storage::get_data_ptr(resid), ncv,
              storage::get_data_ptr(v), N,
              iparam, ipntr, storage::get_data_ptr(workd),
              storage::get_data_ptr(workl), workl_size,
              info);

    storage::destroy(workev, 3*ncv);
    storage::destroy(workl, workl_size);

    if(info) throw std::runtime_error("arpack_worker: dneupd failed with error code " + std::to_string(info));
  }

  // Index of workspace vector, which is expected to be acted on
  inline int from_vector_n() const { return (ipntr[0]-1)/N; }

  // Index of workspace vector, which will receive result of the operator action
  inline int to_vector_n() const { return (ipntr[1]-1)/N; }

  // Get view of a workspace vector
  real_vector_view_t workspace_vector(int n) const {
    if(n < 0 || n > 2)
      throw std::runtime_error("arpack_worker: valid indices of workspace vectors are 0, 1 and 2"
                               " (got " + std::to_string(n) + ")");
    return storage::make_vector_view(workd, n*N, N);
  }

  // Access eigenvalues
  // Cannot be used in ShiftAndInvertReal and ShiftAndInvertImag modes if imag(sigma) != 0.
  complex_vector_t eigenvalues() const {
    assert(!((iparam[6] == ShiftAndInvertReal || iparam[6] == ShiftAndInvertImag) && sigmai != 0));
    return storage::combine_re_im(dr, di, nev);
  }

  // Access Ritz/Schur vectors
  complex_matrix_t eigenvectors() const {
    if(iparam[4] < nev)
      std::cerr << "arpack_worker: only " << iparam[4] << " out of " << nev
                << " Ritz vectors have converged.";
    return storage::extract_eigenvectors(z, di, N, nev);
  }

  // Access residual vector
  real_vector_view_t residual_vector() { return storage::make_vector_view(resid); }

  // Has B*x already been computed?
  bool Bx_available() const { return Bx_available_; }

  // Previously computed vector B*x
  real_vector_const_view_t Bx_vector() const {
    int n = ipntr[2]-1;
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
};

} // namespace ezarpack
