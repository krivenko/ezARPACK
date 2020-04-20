/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2020 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/
/// @file arpack.hpp
/// @brief Declarations of ARPACK-NG FORTRAN 77 subroutines used by ezARPACK
/// and their C++ wrappers.
#pragma once

#include <complex>

namespace ezarpack {

/*! The double precision complex type used in ARPACK-NG calls. */
using dcomplex = std::complex<double>;

/// The ARPACK-NG `*aupd_()` procedures set the output argument `IDO` to one of
/// these values to signal the state of the Reverse Communication Interface
/// (RCI).
enum rci_flag : int {
  Init = 0,           /**< First call of the RCI.                             */
  ApplyOpInit = -1,   /**< Compute `Y = OP * X`
                          (force the starting vector into the range of `OP`). */
  ApplyOp = 1,        /**< Compute `Y = OP * X`.                              */
  ApplyB = 2,         /**< Compute `Y = B * X`.                               */
  Shifts = 3,         /**< Compute and return the shifts for the Implicitly
                           Restarted Arnoldi Method.                          */
  Done = 99           /**< Done with the iterations.                          */
};

namespace f77 {

extern "C" {

/// @name External ARPACK-NG subroutines *aupd()
///
/// These subroutines implement the Reverse Communication Interface for the
/// Implicitly Restarted Arnoldi Iteration.

/// @{

/// @brief External ARPACK-NG FORTRAN 77 subroutine [dsaupd()][dsaupd].
/// [dsaupd]: https://github.com/opencollab/arpack-ng/blob/master/SRC/dsaupd.f
void dsaupd_(int&,          // IDO
             const char*,   // BMAT
             const int&,    // N
             const char*,   // WHICH
             const int&,    // NEV
             const double&, // TOL
             double[],      // RESID
             const int&,    // NCV
             double[],      // V
             const int&,    // LDV
             int[],         // IPARAM
             int[],         // IPNTR
             double[],      // WORKD
             double[],      // WORKL
             const int&,    // LWORKL
             int&);         // INFO

/// @brief External ARPACK-NG FORTRAN 77 subroutine [dnaupd()][dnaupd].
/// [dnaupd]: https://github.com/opencollab/arpack-ng/blob/master/SRC/dnaupd.f
void dnaupd_(int&,          // IDO
             const char*,   // BMAT
             const int&,    // N
             const char*,   // WHICH
             const int&,    // NEV
             const double&, // TOL
             double[],      // RESID
             const int&,    // NCV
             double[],      // V
             const int&,    // LDV
             int[],         // IPARAM
             int[],         // IPNTR
             double[],      // WORKD
             double[],      // WORKL
             const int&,    // LWORKL
             int&);         // INFO

/// @brief External ARPACK-NG FORTRAN 77 subroutine [znaupd()][znaupd].
/// [znaupd]: https://github.com/opencollab/arpack-ng/blob/master/SRC/znaupd.f
void znaupd_(int&,          // IDO
             const char*,   // BMAT
             const int&,    // N
             const char*,   // WHICH
             const int&,    // NEV
             const double&, // TOL
             dcomplex*,     // RESID
             const int&,    // NCV
             dcomplex*,     // V
             const int&,    // LDV
             int[],         // IPARAM
             int[],         // IPNTR
             dcomplex*,     // WORKD
             dcomplex*,     // WORKL
             const int&,    // LWORKL
             double[],      // RWORK
             int&);         // INFO
/// @}

} // extern "C"

/// @brief C++ wrapper for the Reverse Communication Interface ARPACK-NG
/// subroutines `dsaupd_()` and `dnaupd_()`.
/// @tparam Symmetric If `true`, calls `dsaupd_()`, otherwise calls `dnaupd_()`.
/// @param ido Reverse communication flag.
/// @param bmat Specifies the type of the semi-inner product matrix `B`.
/// @param n Dimension of the eigenproblem.
/// @param which Specifies which of the Ritz values of `OP` to compute.
/// @param nev Number of eigenvalues of `OP` to be computed.
/// @param tol Stopping criterion: relative tolerance of the Ritz value.
/// @param resid Initial residual vector.
/// @param ncv How many Lanczos/Arnoldi vectors to generate at each iteration.
/// @param v Contains the final set of Lanczos/Arnoldi basis vectors.
/// @param ldv Leading dimension of `v` as declared in the calling program.
/// @param iparam Array of input/output parameter flags.
/// @param ipntr Pointer to mark the starting locations in the `workd` and
/// `workl` arrays for matrices/vectors used by the Lanczos/Arnoldi iteration.
/// @param workd Array to be used in the basic Lanczos/Arnoldi iteration for
/// reverse communication.
/// @param workl Work array.
/// @param lworkl Dimension of `workl`.
/// @param info Initial residual vector type (input) and error codes (output).
template<bool Symmetric = false>
inline void aupd(rci_flag& ido,
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
    dsaupd_((int&)ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam,
            ipntr, workd, workl, lworkl, info);
  else
    dnaupd_((int&)ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam,
            ipntr, workd, workl, lworkl, info);
}

/// @brief C++ wrapper for the Reverse Communication Interface ARPACK-NG
/// subroutine `znaupd_()`.
/// @param ido Reverse communication flag.
/// @param bmat Specifies the type of the semi-inner product matrix `B`.
/// @param n Dimension of the eigenproblem.
/// @param which Specifies which of the Ritz values of `OP` to compute.
/// @param nev Number of eigenvalues of `OP` to be computed.
/// @param tol Stopping criterion: relative tolerance of the Ritz value.
/// @param resid Initial residual vector.
/// @param ncv How many Arnoldi vectors to generate at each iteration.
/// @param v Contains the final set of Arnoldi basis vectors.
/// @param ldv Leading dimension of `v` as declared in the calling program.
/// @param iparam Array of input/output parameter flags.
/// @param ipntr Pointer to mark the starting locations in the `workd` and
/// `workl` arrays for matrices/vectors used by the Arnoldi iteration.
/// @param workd Array to be used in the basic Arnoldi iteration for
/// reverse communication.
/// @param workl Work array.
/// @param lworkl Dimension of `workl`.
/// @param rwork Work array of dimension `ncv`.
/// @param info Initial residual vector type (input) and error codes (output).
inline void aupd(rci_flag& ido,
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
  znaupd_((int&)ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam,
          ipntr, workd, workl, lworkl, rwork, info);
}

extern "C" {

/// @name External ARPACK-NG subroutines *eupd()
///
/// These subroutines return the converged approximations to eigenvalues
/// of `A*z = lambda*B*z` and (optionally):
///   -# the corresponding approximate eigenvectors;
///   -# an orthonormal basis for the associated approximate invariant subspace;
///   -# both.

/// @{

/// @brief External ARPACK-NG FORTRAN 77 subroutine [dseupd()][dseupd].
/// [dseupd]: https://github.com/opencollab/arpack-ng/blob/master/SRC/dseupd.f
void dseupd_(const int&,    // RVEC
             const char*,   // HOWMNY
             const int[],   // SELECT
             double[],      // D
             double[],      // Z
             const int&,    // LDZ
             const double&, // SIGMA
             const char*,   // BMAT
             const int&,    // N
             const char*,   // WHICH
             const int&,    // NEV
             const double&, // TOL
             double[],      // RESID
             const int&,    // NCV
             double[],      // V
             const int&,    // LDV
             int[],         // IPARAM
             int[],         // IPNTR
             double[],      // WORKD
             double[],      // WORKL
             const int&,    // LWORKL
             int&);         // INFO

/// @brief External ARPACK-NG FORTRAN 77 subroutine [dneupd()][dneupd].
/// [dneupd]: https://github.com/opencollab/arpack-ng/blob/master/SRC/dneupd.f
void dneupd_(const int&,    // RVEC
             const char*,   // HOWMNY
             const int[],   // SELECT
             double[],      // DR
             double[],      // DI
             double[],      // Z
             const int&,    // LDZ
             const double&, // SIGMAR
             const double&, // SIGMAI
             double[],      // WORKEV
             const char*,   // BMAT
             const int&,    // N
             const char*,   // WHICH
             const int&,    // NEV
             const double&, // TOL
             double[],      // RESID
             const int&,    // NCV
             double[],      // V
             const int&,    // LDV
             int[],         // IPARAM
             int[],         // IPNTR
             double[],      // WORKD
             double[],      // WORKL
             const int&,    // LWORKL
             int&);         // INFO

/// @brief External ARPACK-NG FORTRAN 77 subroutine [zneupd()][zneupd].
/// [zneupd]: https://github.com/opencollab/arpack-ng/blob/master/SRC/zneupd.f
void zneupd_(const int&,      // RVEC
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

/// @brief C++ wrapper for the post-processing ARPACK-NG subroutine `dseupd_()`.
/// @param rvec Specifies whether to compute Ritz vectors.
/// @param howmny Specifies how many Ritz vectors are wanted.
/// @param select Selection of Ritz vectors to be computed.
/// @param d Approximate Ritz values (output).
/// @param z Approximate `B`-orthonormal Ritz vectors (output).
/// @param ldz The leading dimension of the array `Z`.
/// @param sigma Shift for spectral transformations.
/// @param bmat Specifies the type of the semi-inner product matrix `B`.
/// @param n Dimension of the eigenproblem.
/// @param which Specifies which of the Ritz values of `OP` to compute.
/// @param nev Number of eigenvalues of `OP` to be computed.
/// @param tol Stopping criterion: relative tolerance of the Ritz value.
/// @param resid Initial residual vector.
/// @param ncv How many Lanczos vectors to generate at each iteration.
/// @param v Contains the final set of Lanczos basis vectors.
/// @param ldv Leading dimension of `v` as declared in the calling program.
/// @param iparam Array of input/output parameter flags.
/// @param ipntr Pointer to mark the starting locations in the `workd` and
/// `workl` arrays for matrices/vectors used by the Lanczos iteration.
/// @param workd Array to be used in the basic Lanczos iteration for
/// reverse communication.
/// @param workl Work array.
/// @param lworkl Dimension of `workl`.
/// @param info Initial residual vector type (input) and error codes (output).
inline void eupd(int rvec,
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
  dseupd_(rvec, howmny, select, d, z, ldz, sigma, bmat, n, which, nev, tol,
          resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
}

/// @brief C++ wrapper for the post-processing ARPACK-NG subroutine `dneupd_()`.
/// @param rvec Specifies whether to compute Ritz/Schur vectors.
/// @param howmny Specifies the form of the basis for the invariant subspace.
/// @param select Selection of Ritz vectors to be computed.
/// @param dr Real parts of approximate Ritz values (output).
/// @param di Imaginary parts of approximate Ritz values (output).
/// @param z Approximate `B`-orthonormal Ritz vectors (output).
/// @param ldz The leading dimension of the array `Z`.
/// @param sigmar Real part of the shift for spectral transformations.
/// @param sigmai Imaginary part of the shift for spectral transformations.
/// @param workev Work array.
/// @param bmat Specifies the type of the semi-inner product matrix `B`.
/// @param n Dimension of the eigenproblem.
/// @param which Specifies which of the Ritz values of `OP` to compute.
/// @param nev Number of eigenvalues of `OP` to be computed.
/// @param tol Stopping criterion: relative tolerance of the Ritz value.
/// @param resid Initial residual vector.
/// @param ncv How many Arnoldi vectors to generate at each iteration.
/// @param v Contains the final set of Arnoldi basis vectors.
/// @param ldv Leading dimension of `v` as declared in the calling program.
/// @param iparam Array of input/output parameter flags.
/// @param ipntr Pointer to mark the starting locations in the `workd` and
/// `workl` arrays for matrices/vectors used by the Arnoldi iteration.
/// @param workd Array to be used in the basic Arnoldi iteration for
/// reverse communication.
/// @param workl Work array.
/// @param lworkl Dimension of `workl`.
/// @param info Initial residual vector type (input) and error codes (output).
inline void eupd(int rvec,
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
  dneupd_(rvec, howmny, select, dr, di, z, ldz, sigmar, sigmai, workev, bmat, n,
          which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl,
          lworkl, info);
}

/// @brief C++ wrapper for the post-processing ARPACK-NG subroutine `zneupd_()`.
/// @param rvec Specifies whether to compute Ritz/Schur vectors.
/// @param howmny Specifies the form of the basis for the invariant subspace.
/// @param select Selection of Ritz vectors to be computed.
/// @param d Approximate Ritz values (output).
/// @param z Approximate `B`-orthonormal Ritz vectors (output).
/// @param ldz The leading dimension of the array `Z`.
/// @param sigma Shift for spectral transformations.
/// @param workev Work array.
/// @param bmat Specifies the type of the semi-inner product matrix `B`.
/// @param n Dimension of the eigenproblem.
/// @param which Specifies which of the Ritz values of `OP` to compute.
/// @param nev Number of eigenvalues of `OP` to be computed.
/// @param tol Stopping criterion: relative tolerance of the Ritz value.
/// @param resid Initial residual vector.
/// @param ncv How many Arnoldi vectors to generate at each iteration.
/// @param v Contains the final set of Arnoldi basis vectors.
/// @param ldv Leading dimension of `v` as declared in the calling program.
/// @param iparam Array of input/output parameter flags.
/// @param ipntr Pointer to mark the starting locations in the `workd` and
/// `workl` arrays for matrices/vectors used by the Arnoldi iteration.
/// @param workd Array to be used in the basic Arnoldi iteration for
/// reverse communication.
/// @param workl Work array.
/// @param lworkl Dimension of `workl`.
/// @param rwork Work array of dimension `ncv`.
/// @param info Initial residual vector type (input) and error codes (output).
inline void eupd(int rvec,
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
  zneupd_(rvec, howmny, select, d, z, ldz, sigma, workev, bmat, n, which, nev,
          tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork,
          info);
}

} // namespace f77

} // namespace ezarpack
