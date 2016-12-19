/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2016 I. Krivenko
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once

namespace triqs { namespace arrays { namespace arpack {

using namespace triqs::arrays;

/*************************
 * Real symmetric matrix A
 *************************/
template<> class arpack_worker<Symmetric> {

 int N;                 // Matrix size
 const char* which;     // WHICH parameter
 int nev;               // Number of eigenvalues
 int tol;               // Relative tolerance for Ritz value convergence
 vector<double> resid;  // Residual vector
 int ncv;               // Number of Lanczos vectors to be generated
 matrix<double> v;      // Matrix with Lanczos basis vectors
 vector<double> d;      // Ritz values
 int iparam[11];        // Various input/output parameters
 int ipntr[11];         // Pointer to mark the starting locations in the workd and workl
 vector<double> workd;  // Working space
 int info = 0;          // !=0 to use resid, 0 otherwise
 int rvec;              // RVEC parameter of dseupd
 vector<int> select;    // SELECT parameter of dseupd

public:

 struct params_t {
  // Number of eigenvalues (Ritz values) to compute
  unsigned int n_eigenvalues;
  // Which of the Ritz values to compute
  enum {Largest, Smallest, LargestMagnitude, SmallestMagnitude, BothEnds} eigenvalues_select;

  // Expert option: number of Lanczos vectors to be generated
  // default: min(2*n_eigenvalues, N)
  int ncv = -1;

  // Compute Ritz vectors?
  bool compute_eigenvectors;

  // Use a randomly generated initial residual vector
  bool random_residual_vector = true;
  // Initial residual vector
  vector<double> init_residual_vector = {};

  // Eigenvalue shift
  double sigma = 0;

  // relative tolerance for Ritz value convergence
  double tolerance = 0;
  // Maximum number of Arnoldi update iterations
  unsigned int max_iter = INT_MAX;

  params_t(unsigned int n_eigenvalues, decltype(eigenvalues_select) ev_select, bool compute_eigenvectors) :
   n_eigenvalues(n_eigenvalues),
   eigenvalues_select(ev_select),
   compute_eigenvectors(compute_eigenvectors)
  {}
 };

 arpack_worker(unsigned int N) : N(N), resid(N), workd(3*N), v(0,0,FORTRAN_LAYOUT) {
  iparam[3] = 1;
 }

 // Prepare values of input parameters
 void prepare(params_t const& params) {

  // Check n_eigenvalues
  nev = params.n_eigenvalues;
  unsigned int nev_min = params.eigenvalues_select == params_t::BothEnds ? 2 : 1;
  unsigned int nev_max = N-1;

  if(nev < nev_min || nev > nev_max)
   TRIQS_RUNTIME_ERROR << "arpack_worker: n_eigenvalues must be within [" << nev_min
                       << ";" << nev_max << "]";

  // Character codes for eigenvalues_select
  static const std::array<const char*,5> wh = {"LA","SA","LM","SM","BE"};
  which = wh[int(params.eigenvalues_select)];

  // Check ncv
  ncv = params.ncv;
  if(ncv == -1) ncv = std::min(2*int(params.n_eigenvalues)+2, N);
  else if(ncv <= params.n_eigenvalues || ncv > N)
   TRIQS_RUNTIME_ERROR << "arpack_worker: ncv must be within ]" << params.n_eigenvalues
                       << ";" << N << "]";
  v.resize(N,ncv);

  // Eigenvectors
  rvec = params.compute_eigenvectors;
  select.resize(ncv);

  // Tolerance
  tol = std::max(.0, params.tolerance);
  if(params.tolerance < 0) {
   std::cerr << "arpack_worker: negative tolerance " << params.tolerance
             << " is interpreted as machine epsilon." << std::endl;
  }

  // Initial residual vector
  if(params.random_residual_vector) {
   info = 0;
  } else {
   info = 1;
   if(params.init_residual_vector.size() != N)
    TRIQS_RUNTIME_ERROR << "arpack_worker: initial residual vector of a wrong size "
                        << params.init_residual_vector.size()
                        << " (must be " << N << ")";
    resid = params.init_residual_vector;
  }

  iparam[2] = int(params.max_iter); // Max number of iterations
  if(iparam[2] <= 0)
   TRIQS_RUNTIME_ERROR << "arpack_worker: maximum number of Arnoldi update iterations must be positive";
 }

 struct trivial_shifts_f { void operator()(vector_view<double> shifts){}};

 /**********************************
 * Solve the standard eigenproblem
 * A*x = \lambda*x
 **********************************/

 // a: callable taking 4 arguments
 // a(vector_view<double> from, int from_n, vector_view<double> to, int to_n)
 // 'a' is expected to act on 'from' and write the result to 'to': to = A*from
 // 'from' is also indirectly available as this->workspace_vector(from_n)
 // 'to' is also indirectly available as this->workspace_vector(to_n)
 //
 // shifts_f: callable taking one argument
 // shifts_f(vector_view<double> shifts)
 // 'shifts_f' is expected to place the shifts for implicit restart into 'shifts'
 template<typename A, typename ShiftsF = trivial_shifts_f>
 void operator()(A a, params_t const& params, ShiftsF shifts_f = {}) {

  prepare(params);

  iparam[0] = (std::is_same<ShiftsF,trivial_shifts_f>::value ? 1 : 0);
  iparam[6] = 1; // Mode 1, standard eigenproblem

  vector<double> workl(ncv*ncv + 8*ncv);

  rci_flag ido = Init;
  do {
   f77::aupd<true>((int&)ido, "I", N, which,
                   nev, tol, resid.data_start(), ncv,
                   v.data_start(), first_dim(v),
                   iparam, ipntr, workd.data_start(),
                   workl.data_start(), workl.size(),
                   info);
   switch(ido) {
    case ApplyOpInit:
    case ApplyOp: {
      int from_n = ipntr[0]-1;
      int to_n = ipntr[1]-1;
      a(workd(range(from_n,from_n+N)),from_n/N,workd(range(to_n,to_n+N)),to_n/N);
     }
     break;
    case Shifts:
     shifts_f(workl(range(ipntr[10]-1,ipntr[10]-1+iparam[7])));
     break;
    case Done: break;
    default: TRIQS_RUNTIME_ERROR << "arpack_worker: reverse communication interface error";
   }
  } while(ido != Done);

  switch(info) {
   case 0: break;
   case 1: throw(maxiter_reached(iparam[2]));
   case 3: throw(ncv_insufficient(ncv));
   default: TRIQS_RUNTIME_ERROR << "arpack_worker: dsaupd failed with error code " << info;
  }

  d.resize(nev);
  double sigma;

  f77::eupd(rvec, "A", select.data_start(),
            d.data_start(), v.data_start(), first_dim(v),
            sigma, "I", N, which,
            nev, tol, resid.data_start(), ncv,
            v.data_start(), first_dim(v),
            iparam, ipntr, workd.data_start(),
            workl.data_start(), workl.size(),
            info);

  if(info) TRIQS_RUNTIME_ERROR << "arpack_worker: dseupd failed with error code " << info;
 }

 /***********************************
 * Solve the generalized eigenproblem
 * A*x = \lambda*M*x
 ************************************/

 enum Mode : int {Invert = 2,           // OP = inv[M]*A and B = M
                  ShiftAndInvert = 3,   // OP = (inv[A - sigma*M])*M and B = M
                  Buckling = 4,         // OP = (inv[A - sigma*M])*A and B = A
                  Cayley = 5};          // OP = inv[A - sigma*M]*[A + sigma*M] and B = M

 // op: callable taking 5 arguments
 // op(vector_view<double> from, int from_n, vector_view<double> to, int to_n, bool Bx_available)
 // 'op' is expected to act on 'from' and write the result to 'to': to = OP*from
 // 'from' is also indirectly available as this->workspace_vector(from_n)
 // 'to' is also indirectly available as this->workspace_vector(to_n)
 // 'Bx_available' indicate whether B*x has already been computed and available as this->Bx_vector()
 //
 // b: callable taking 4 arguments
 // 'b' is expected to act on 'from' and write the result to 'to': to = B*from
 // 'from' is also indirectly available as this->workspace_vector(from_n)
 // 'to' is also indirectly available as this->workspace_vector(to_n)
 //
 // shifts_f: callable taking one argument
 // shifts_f(vector_view<double> shifts)
 // 'shifts_f' is expected to place the shifts for implicit restart into 'shifts'
 template<typename OP, typename B, typename ShiftsF = trivial_shifts_f>
 void operator()(OP op, B b, Mode mode, params_t const& params, ShiftsF shifts_f = {}) {

  prepare(params);

  iparam[0] = (std::is_same<ShiftsF,trivial_shifts_f>::value ? 1 : 0);
  iparam[6] = mode; // Modes 2-5, generalized eigenproblem

  vector<double> workl(ncv*ncv + 8*ncv);

  rci_flag ido = Init;
  do {
   f77::aupd<true>((int&)ido, "G", N, which,
                   nev, tol, resid.data_start(), ncv,
                   v.data_start(), first_dim(v),
                   iparam, ipntr, workd.data_start(),
                   workl.data_start(), workl.size(),
                   info);
   switch(ido) {
    case ApplyOpInit: {
      int from_n = ipntr[0]-1;
      int to_n = ipntr[1]-1;
      op(workd(range(from_n,from_n+N)),from_n/N,workd(range(to_n,to_n+N)),to_n/N,false);
     }
     break;
    case ApplyOp: {
      int from_n = ipntr[0]-1;
      int to_n = ipntr[1]-1;
      // B*x is available via Bx_vector()
      op(workd(range(from_n,from_n+N)),from_n/N,workd(range(to_n,to_n+N)),to_n/N,true);
     }
     break;
    case ApplyB: {
      int from_n = ipntr[0]-1;
      int to_n = ipntr[1]-1;
      b(workd(range(from_n,from_n+N)),from_n/N,workd(range(to_n,to_n+N)),to_n/N);
     }
     break;
    case Shifts:
     shifts_f(workl(range(ipntr[10]-1,ipntr[10]-1+iparam[7])));
     break;
    case Done: break;
    default: TRIQS_RUNTIME_ERROR << "arpack_worker: reverse communication interface error";
   }
  } while(ido != Done);

  switch(info) {
   case 0: break;
   case 1: throw(maxiter_reached(iparam[2]));
   case 3: throw(ncv_insufficient(ncv));
   default: TRIQS_RUNTIME_ERROR << "arpack_worker: dsaupd failed with error code " << info;
  }

  d.resize(nev);
  double sigma;
  if(mode != Invert) sigma = params.sigma;

  f77::eupd(rvec, "A", select.data_start(),
            d.data_start(), v.data_start(), first_dim(v),
            sigma, "G", N, which,
            nev, tol, resid.data_start(), ncv,
            v.data_start(), first_dim(v),
            iparam, ipntr, workd.data_start(),
            workl.data_start(), workl.size(),
            info);

  if(info) TRIQS_RUNTIME_ERROR << "arpack_worker: dseupd failed with error code " << info;
 }

 // Get view of a workspace vector
 vector_view<double> workspace_vector(int n) const {
  if(n < 0 || n > 2)
   TRIQS_RUNTIME_ERROR << "arpack_worker: valid indices of workspace vectors are 0, 1 and 2"
                       << " (got " << n << ")";
  return workd(range(n*N,(n+1)*N));
 }

 // Access eigenvalues
 vector_view<double> eigenvalues() const { return d; }

 // Access eigenvectors
 matrix_view<double> eigenvectors() const { return v(range(),range(nev)); }

 // Access residual vector
 vector_view<double> residual_vector() const { return resid; }

 // Previously computed vector B*x
 vector_view<double> Bx_vector() const {
  int n = ipntr[2]-1;
  return workd(range(n,n+N));
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

}}}
