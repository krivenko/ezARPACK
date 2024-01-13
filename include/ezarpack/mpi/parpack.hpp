/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2024 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/
/// @file ezarpack/mpi/parpack.hpp
/// @brief Declarations of parallel ARPACK-NG FORTRAN 77 subroutines used by
/// ezARPACK and their C++ wrappers.
#pragma once

#include <mpi.h>

#include "../common.hpp"

namespace ezarpack {
namespace mpi {
namespace f77 {

extern "C" {

/// @name External ARPACK-NG subroutines p*aupd()
///
/// These subroutines implement the Reverse Communication Interface for the
/// Implicitly Restarted Arnoldi Iteration in MPI parallelized mode.

/// @{

/// @brief External ARPACK-NG FORTRAN 77 subroutine [pdsaupd()][pdsaupd].
/// [pdsaupd]:
/// http://github.com/opencollab/arpack-ng/blob/master/PARPACK/SRC/MPI/pdsaupd.f
void pdsaupd_(const MPI_Fint&, // COMM
              int&,            // IDO
              const char*,     // BMAT
              const int&,      // N
              const char*,     // WHICH
              const int&,      // NEV
              const double&,   // TOL
              double[],        // RESID
              const int&,      // NCV
              double[],        // V
              const int&,      // LDV
              int[],           // IPARAM
              int[],           // IPNTR
              double[],        // WORKD
              double[],        // WORKL
              const int&,      // LWORKL
              int&);           // INFO

/// @brief External ARPACK-NG FORTRAN 77 subroutine [pdnaupd()][pdnaupd].
/// [pdnaupd]:
/// http://github.com/opencollab/arpack-ng/blob/master/PARPACK/SRC/MPI/pdnaupd.f
void pdnaupd_(const MPI_Fint&, // COMM
              int&,            // IDO
              const char*,     // BMAT
              const int&,      // N
              const char*,     // WHICH
              const int&,      // NEV
              const double&,   // TOL
              double[],        // RESID
              const int&,      // NCV
              double[],        // V
              const int&,      // LDV
              int[],           // IPARAM
              int[],           // IPNTR
              double[],        // WORKD
              double[],        // WORKL
              const int&,      // LWORKL
              int&);           // INFO

/// @brief External ARPACK-NG FORTRAN 77 subroutine [pznaupd()][pznaupd].
/// [pznaupd]:
/// http://github.com/opencollab/arpack-ng/blob/master/PARPACK/SRC/MPI/pznaupd.f
void pznaupd_(const MPI_Fint&, // COMM
              int&,            // IDO
              const char*,     // BMAT
              const int&,      // N
              const char*,     // WHICH
              const int&,      // NEV
              const double&,   // TOL
              dcomplex*,       // RESID
              const int&,      // NCV
              dcomplex*,       // V
              const int&,      // LDV
              int[],           // IPARAM
              int[],           // IPNTR
              dcomplex*,       // WORKD
              dcomplex*,       // WORKL
              const int&,      // LWORKL
              double[],        // RWORK
              int&);           // INFO
/// @}

} // extern "C"

/// @brief C++ wrapper for the Reverse Communication Interface ARPACK-NG
/// subroutines `pdsaupd_()` and `pdnaupd_()`.
/// @tparam Symmetric If `true`, calls `pdsaupd_()`, otherwise calls
/// `pdnaupd_()`.
/// @param comm MPI communicator.
/// @param ido Reverse communication flag.
/// @param bmat Specifies the type of the semi-inner product matrix `B`.
/// @param n Dimension of the MPI rank-local vector block.
/// @param which Specifies which of the Ritz values of `OP` to compute.
/// @param nev Number of eigenvalues of `OP` to be computed.
/// @param tol Stopping criterion: relative tolerance of the Ritz value.
/// @param resid MPI rank-local block of the initial residual vector.
/// @param ncv How many Lanczos/Arnoldi vectors to generate at each iteration.
/// @param v Contains the final set of Lanczos/Arnoldi basis vectors
/// (MPI rank-local portion).
/// @param ldv Leading dimension of `v` as declared in the calling program.
/// @param iparam Array of input/output parameter flags.
/// @param ipntr Pointer to mark the starting locations in the `workd` and
/// `workl` arrays for matrices/vectors used by the Lanczos/Arnoldi iteration.
/// @param workd Array to be used in the basic Lanczos/Arnoldi iteration for
/// reverse communication (MPI rank-local portion).
/// @param workl Work array.
/// @param lworkl Dimension of `workl`.
/// @param info Initial residual vector type (input) and error codes (output).
template<bool Symmetric = false>
inline void paupd(const MPI_Comm& comm,
                  rci_flag& ido,
                  const char* bmat,
                  int n,
                  const char* which,
                  int nev,
                  double tol,
                  double* resid,
                  int ncv,
                  double* v,
                  int ldv,
                  int* iparam,
                  int* ipntr,
                  double* workd,
                  double* workl,
                  int lworkl,
                  int& info) {
  if(Symmetric)
    pdsaupd_(MPI_Comm_c2f(comm), (int&)ido, bmat, n, which, nev, tol, resid,
             ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
  else
    pdnaupd_(MPI_Comm_c2f(comm), (int&)ido, bmat, n, which, nev, tol, resid,
             ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
}

/// @brief C++ wrapper for the Reverse Communication Interface ARPACK-NG
/// subroutine `pznaupd_()`.
/// @param comm MPI communicator.
/// @param ido Reverse communication flag.
/// @param bmat Specifies the type of the semi-inner product matrix `B`.
/// @param n Dimension of the MPI rank-local vector block.
/// @param which Specifies which of the Ritz values of `OP` to compute.
/// @param nev Number of eigenvalues of `OP` to be computed.
/// @param tol Stopping criterion: relative tolerance of the Ritz value.
/// @param resid MPI rank-local block of the initial residual vector.
/// @param ncv How many Arnoldi vectors to generate at each iteration.
/// @param v Contains the final set of Arnoldi basis vectors (MPI rank-local
/// portion).
/// @param ldv Leading dimension of `v` as declared in the calling program.
/// @param iparam Array of input/output parameter flags.
/// @param ipntr Pointer to mark the starting locations in the `workd` and
/// `workl` arrays for matrices/vectors used by the Arnoldi iteration.
/// @param workd Array to be used in the basic Arnoldi iteration for
/// reverse communication (MPI rank-local portion).
/// @param workl Work array.
/// @param lworkl Dimension of `workl`.
/// @param rwork Work array of dimension `ncv`.
/// @param info Initial residual vector type (input) and error codes (output).
inline void paupd(const MPI_Comm& comm,
                  rci_flag& ido,
                  const char* bmat,
                  int n,
                  const char* which,
                  int nev,
                  double tol,
                  dcomplex* resid,
                  int ncv,
                  dcomplex* v,
                  int ldv,
                  int* iparam,
                  int* ipntr,
                  dcomplex* workd,
                  dcomplex* workl,
                  int lworkl,
                  double* rwork,
                  int& info) {
  pznaupd_(MPI_Comm_c2f(comm), (int&)ido, bmat, n, which, nev, tol, resid, ncv,
           v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, info);
}

extern "C" {

/// @name External ARPACK-NG subroutines p*eupd()
///
/// These MPI parallelized subroutines return the converged approximations to
/// eigenvalues of @f$ \hat O \mathbf{z} = \lambda \hat B\mathbf{z} @f$
/// and (optionally):
///   -# the corresponding approximate eigenvectors;
///   -# an orthonormal basis for the associated approximate invariant subspace;
///   -# both.

/// @{

/// @brief External ARPACK-NG FORTRAN 77 subroutine [pdseupd()][pdseupd].
/// [pdseupd]:
/// http://github.com/opencollab/arpack-ng/blob/master/PARPACK/SRC/MPI/pdseupd.f
void pdseupd_(const MPI_Fint&, // COMM
              const int&,      // RVEC
              const char*,     // HOWMNY
              const int[],     // SELECT
              double[],        // D
              double[],        // Z
              const int&,      // LDZ
              const double&,   // SIGMA
              const char*,     // BMAT
              const int&,      // N
              const char*,     // WHICH
              const int&,      // NEV
              const double&,   // TOL
              double[],        // RESID
              const int&,      // NCV
              double[],        // V
              const int&,      // LDV
              int[],           // IPARAM
              int[],           // IPNTR
              double[],        // WORKD
              double[],        // WORKL
              const int&,      // LWORKL
              int&);           // INFO

/// @brief External ARPACK-NG FORTRAN 77 subroutine [pdneupd()][pdneupd].
/// [pdneupd]:
/// http://github.com/opencollab/arpack-ng/blob/master/PARPACK/SRC/MPI/pdneupd.f
void pdneupd_(const MPI_Fint&, // COMM
              const int&,      // RVEC
              const char*,     // HOWMNY
              const int[],     // SELECT
              double[],        // DR
              double[],        // DI
              double[],        // Z
              const int&,      // LDZ
              const double&,   // SIGMAR
              const double&,   // SIGMAI
              double[],        // WORKEV
              const char*,     // BMAT
              const int&,      // N
              const char*,     // WHICH
              const int&,      // NEV
              const double&,   // TOL
              double[],        // RESID
              const int&,      // NCV
              double[],        // V
              const int&,      // LDV
              int[],           // IPARAM
              int[],           // IPNTR
              double[],        // WORKD
              double[],        // WORKL
              const int&,      // LWORKL
              int&);           // INFO

/// @brief External ARPACK-NG FORTRAN 77 subroutine [pzneupd()][pzneupd].
/// [pzneupd]:
/// http://github.com/opencollab/arpack-ng/blob/master/PARPACK/SRC/MPI/pzneupd.f
void pzneupd_(const MPI_Fint&, // COMM
              const int&,      // RVEC
              const char*,     // HOWMNY
              const int[],     // SELECT
              dcomplex*,       // D
              dcomplex*,       // Z
              const int&,      // LDZ
              const dcomplex&, // SIGMA
              dcomplex*,       // WORKEV
              const char*,     // BMAT
              const int&,      // N
              const char*,     // WHICH
              const int&,      // NEV
              const double&,   // TOL
              dcomplex*,       // RESID
              const int&,      // NCV
              dcomplex*,       // V
              const int&,      // LDV
              int[],           // IPARAM
              int[],           // IPNTR
              dcomplex*,       // WORKD
              dcomplex*,       // WORKL
              const int&,      // LWORKL
              double[],        // RWORK
              int&);           // INFO

/// @}

} // extern "C"

/// @brief C++ wrapper for the post-processing ARPACK-NG subroutine
/// `pdseupd_()`.
/// @param comm MPI communicator.
/// @param rvec Specifies whether to compute Ritz vectors.
/// @param howmny Specifies how many Ritz vectors are wanted.
/// @param select Selection of Ritz vectors to be computed.
/// @param d Approximate Ritz values (output).
/// @param z Approximate `B`-orthonormal Ritz vectors (output, MPI rank-local
/// portion).
/// @param ldz The leading dimension of the array `Z`.
/// @param sigma Shift for spectral transformations.
/// @param bmat Specifies the type of the semi-inner product matrix `B`.
/// @param n Dimension of the MPI rank-local vector block.
/// @param which Specifies which of the Ritz values of `OP` to compute.
/// @param nev Number of eigenvalues of `OP` to be computed.
/// @param tol Stopping criterion: relative tolerance of the Ritz value.
/// @param resid MPI rank-local block of the initial residual vector.
/// @param ncv How many Lanczos vectors to generate at each iteration.
/// @param v Contains the final set of Lanczos basis vectors (MPI rank-local
/// portion).
/// @param ldv Leading dimension of `v` as declared in the calling program.
/// @param iparam Array of input/output parameter flags.
/// @param ipntr Pointer to mark the starting locations in the `workd` and
/// `workl` arrays for matrices/vectors used by the Lanczos iteration.
/// @param workd Array to be used in the basic Lanczos iteration for
/// reverse communication (MPI rank-local portion).
/// @param workl Work array.
/// @param lworkl Dimension of `workl`.
/// @param info Initial residual vector type (input) and error codes (output).
inline void peupd(const MPI_Comm& comm,
                  int rvec,
                  const char* howmny,
                  const int* select,
                  double* d,
                  double* z,
                  int ldz,
                  double sigma,
                  const char* bmat,
                  int n,
                  const char* which,
                  int nev,
                  double tol,
                  double* resid,
                  int ncv,
                  double* v,
                  int ldv,
                  int* iparam,
                  int* ipntr,
                  double* workd,
                  double* workl,
                  int lworkl,
                  int& info) {
  pdseupd_(MPI_Comm_c2f(comm), rvec, howmny, select, d, z, ldz, sigma, bmat, n,
           which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl,
           lworkl, info);
}

/// @brief C++ wrapper for the post-processing ARPACK-NG subroutine
/// `pdneupd_()`.
/// @param comm MPI communicator.
/// @param rvec Specifies whether to compute Ritz/Schur vectors.
/// @param howmny Specifies the form of the basis for the invariant subspace.
/// @param select Selection of Ritz vectors to be computed.
/// @param dr Real parts of approximate Ritz values (output).
/// @param di Imaginary parts of approximate Ritz values (output).
/// @param z Approximate `B`-orthonormal Ritz vectors (output, MPI rank-local
/// portion).
/// @param ldz The leading dimension of the array `Z`.
/// @param sigmar Real part of the shift for spectral transformations.
/// @param sigmai Imaginary part of the shift for spectral transformations.
/// @param workev Work array.
/// @param bmat Specifies the type of the semi-inner product matrix `B`.
/// @param n Dimension of the MPI rank-local vector block.
/// @param which Specifies which of the Ritz values of `OP` to compute.
/// @param nev Number of eigenvalues of `OP` to be computed.
/// @param tol Stopping criterion: relative tolerance of the Ritz value.
/// @param resid MPI rank-local block of the initial residual vector.
/// @param ncv How many Arnoldi vectors to generate at each iteration.
/// @param v Contains the final set of Arnoldi basis vectors (MPI rank-local
/// portion).
/// @param ldv Leading dimension of `v` as declared in the calling program.
/// @param iparam Array of input/output parameter flags.
/// @param ipntr Pointer to mark the starting locations in the `workd` and
/// `workl` arrays for matrices/vectors used by the Arnoldi iteration.
/// @param workd Array to be used in the basic Arnoldi iteration for
/// reverse communication (MPI rank-local portion).
/// @param workl Work array.
/// @param lworkl Dimension of `workl`.
/// @param info Initial residual vector type (input) and error codes (output).
inline void peupd(const MPI_Comm& comm,
                  int rvec,
                  const char* howmny,
                  const int* select,
                  double* dr,
                  double* di,
                  double* z,
                  int ldz,
                  double sigmar,
                  double sigmai,
                  double* workev,
                  const char* bmat,
                  int n,
                  const char* which,
                  int nev,
                  double tol,
                  double* resid,
                  int ncv,
                  double* v,
                  int ldv,
                  int* iparam,
                  int* ipntr,
                  double* workd,
                  double* workl,
                  int lworkl,
                  int& info) {
  pdneupd_(MPI_Comm_c2f(comm), rvec, howmny, select, dr, di, z, ldz, sigmar,
           sigmai, workev, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam,
           ipntr, workd, workl, lworkl, info);
}

/// @brief C++ wrapper for the post-processing ARPACK-NG subroutine
/// `pzneupd_()`.
/// @param comm MPI communicator.
/// @param rvec Specifies whether to compute Ritz/Schur vectors.
/// @param howmny Specifies the form of the basis for the invariant subspace.
/// @param select Selection of Ritz vectors to be computed.
/// @param d Approximate Ritz values (output).
/// @param z Approximate `B`-orthonormal Ritz vectors (output, MPI rank-local
/// portion).
/// @param ldz The leading dimension of the array `Z`.
/// @param sigma Shift for spectral transformations.
/// @param workev Work array.
/// @param bmat Specifies the type of the semi-inner product matrix `B`.
/// @param n Dimension of the MPI rank-local vector block.
/// @param which Specifies which of the Ritz values of `OP` to compute.
/// @param nev Number of eigenvalues of `OP` to be computed.
/// @param tol Stopping criterion: relative tolerance of the Ritz value.
/// @param resid MPI rank-local block of the initial residual vector.
/// @param ncv How many Arnoldi vectors to generate at each iteration.
/// @param v Contains the final set of Arnoldi basis vectors (MPI rank-local
/// portion).
/// @param ldv Leading dimension of `v` as declared in the calling program.
/// @param iparam Array of input/output parameter flags.
/// @param ipntr Pointer to mark the starting locations in the `workd` and
/// `workl` arrays for matrices/vectors used by the Arnoldi iteration.
/// @param workd Array to be used in the basic Arnoldi iteration for
/// reverse communication (MPI rank-local portion).
/// @param workl Work array.
/// @param lworkl Dimension of `workl`.
/// @param rwork Work array of dimension `ncv`.
/// @param info Initial residual vector type (input) and error codes (output).
inline void peupd(const MPI_Comm& comm,
                  int rvec,
                  const char* howmny,
                  const int* select,
                  dcomplex* d,
                  dcomplex* z,
                  int ldz,
                  dcomplex sigma,
                  dcomplex* workev,
                  const char* bmat,
                  int n,
                  const char* which,
                  int nev,
                  double tol,
                  dcomplex* resid,
                  int ncv,
                  dcomplex* v,
                  int ldv,
                  int* iparam,
                  int* ipntr,
                  dcomplex* workd,
                  dcomplex* workl,
                  int lworkl,
                  double* rwork,
                  int& info) {
  pzneupd_(MPI_Comm_c2f(comm), rvec, howmny, select, d, z, ldz, sigma, workev,
           bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd,
           workl, lworkl, rwork, info);
}

} // namespace f77
} // namespace mpi
} // namespace ezarpack
