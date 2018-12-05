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

#include <complex>

namespace ezarpack {

 using dcomplex = std::complex<double>;

 namespace f77 {

  // Reverse communication interface for the Implicitly Restarted Arnoldi Iteration

  extern "C" {
   // https://github.com/opencollab/arpack-ng/blob/master/SRC/dsaupd.f
   void dsaupd_(int &,            // IDO
                const char *,     // BMAT
                const int &,      // N
                const char *,     // WHICH
                const int &,      // NEV
                const double &,   // TOL
                double [],        // RESID
                const int &,      // NCV
                double [],        // V
                const int &,      // LDV
                int [],           // IPARAM
                int [],           // IPNTR
                double [],        // WORKD
                double [],        // WORKL
                const int &,      // LWORKL
                int &);           // INFO

   // https://github.com/opencollab/arpack-ng/blob/master/SRC/dnaupd.f
   void dnaupd_(int &,            // IDO
                const char *,     // BMAT
                const int &,      // N
                const char *,     // WHICH
                const int &,      // NEV
                const double &,   // TOL
                double [],        // RESID
                const int &,      // NCV
                double [],        // V
                const int &,      // LDV
                int [],           // IPARAM
                int [],           // IPNTR
                double [],        // WORKD
                double [],        // WORKL
                const int &,      // LWORKL
                int &);           // INFO

   // https://github.com/opencollab/arpack-ng/blob/master/SRC/znaupd.f
   void znaupd_(int &,                  // IDO
                const char *,           // BMAT
                const int &,            // N
                const char *,           // WHICH
                const int &,            // NEV
                const double &,         // TOL
                dcomplex *, // RESID
                const int &,            // NCV
                dcomplex *, // V
                const int &,            // LDV
                int [],                 // IPARAM
                int [],                 // IPNTR
                dcomplex *, // WORKD
                dcomplex *, // WORKL
                const int &,            // LWORKL
                double [],              // RWORK
                int &);                 // INFO

  } // extern "C"

  template<bool Symmetric = false>
  inline void aupd(int & ido, const char * bmat, int n, const char * which, int nev, double tol,
                   double * resid, int ncv, double * v, int ldv, int * iparam, int * ipntr,
                   double * workd, double * workl, int lworkl, int & info) {
   if(Symmetric)
    dsaupd_(ido,bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info);
   else
    dnaupd_(ido,bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info);
  }
  inline void aupd(int & ido, const char * bmat, int n, const char * which, int nev, double tol,
                   dcomplex * resid, int ncv, dcomplex * v, int ldv, int * iparam, int * ipntr,
                   dcomplex * workd, dcomplex * workl, int lworkl, double * rwork, int & info) {
   znaupd_(ido,bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,rwork,info);
  }

  /* This subroutine returns the converged approximations to eigenvalues
   * of A*z = lambda*B*z and (optionally):
   *
   * (1) the corresponding approximate eigenvectors,
   *
   * (2) an orthonormal (Lanczos) basis for the associated approximate invariant subspace,
   *
   * (3) Both.
   */

  extern "C" {

  // https://github.com/opencollab/arpack-ng/blob/master/SRC/dseupd.f
  void dseupd_(const int &,      // RVEC
               const char *,     // HOWMNY
               const int [],     // SELECT
               double [],        // D
               double [],        // Z
               const int &,      // LDZ
               const double &,   // SIGMA
               const char *,     // BMAT
               const int &,      // N
               const char *,     // WHICH
               const int &,      // NEV
               const double &,   // TOL
               double [],        // RESID
               const int &,      // NCV
               double [],        // V
               const int &,      // LDV
               int [],           // IPARAM
               int [],           // IPNTR
               double [],        // WORKD
               double [],        // WORKL
               const int &,      // LWORKL
               int &);           // INFO

  // https://github.com/opencollab/arpack-ng/blob/master/SRC/dneupd.f
  void dneupd_(const int &,      // RVEC
               const char *,     // HOWMNY
               const int [],     // SELECT
               double [],        // DR
               double [],        // DI
               double [],        // Z
               const int &,      // LDZ
               const double &,   // SIGMAR
               const double &,   // SIGMAI
               double [],        // WORKEV
               const char *,     // BMAT
               const int &,      // N
               const char *,     // WHICH
               const int &,      // NEV
               const double &,   // TOL
               double [],        // RESID
               const int &,      // NCV
               double [],        // V
               const int &,      // LDV
               int [],           // IPARAM
               int [],           // IPNTR
               double [],        // WORKD
               double [],        // WORKL
               const int &,      // LWORKL
               int &);           // INFO

  // https://github.com/opencollab/arpack-ng/blob/master/SRC/zneupd.f
  void zneupd_(const int &,      // RVEC
               const char *,     // HOWMNY
               const int [],     // SELECT
               dcomplex *,       // D
               dcomplex *,       // Z
               const int &,      // LDZ
               const dcomplex &, // SIGMA
               dcomplex *,       // WORKEV
               const char *,     // BMAT
               const int &,      // N
               const char *,     // WHICH
               const int &,      // NEV
               const double &,   // TOL
               dcomplex *,       // RESID
               const int &,      // NCV
               dcomplex *,       // V
               const int &,      // LDV
               int [],           // IPARAM
               int [],           // IPNTR
               dcomplex *,       // WORKD
               dcomplex *,       // WORKL
               const int &,      // LWORKL
               double [],        // RWORK
               int &);           // INFO
  } // extern "C"

  inline void eupd(int rvec, const char * howmny, const int * select, double * d, double * z, int ldz, double sigma,
                   const char * bmat, int n, const char * which, int nev, double tol, double * resid, int ncv,
                   double * v, int ldv, int * iparam, int * ipntr, double * workd, double * workl, int lworkl,
                   int & info) {
    dseupd_(rvec,howmny,select,d,z,ldz,sigma,
            bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info);
  }
  inline void eupd(int rvec, const char * howmny, const int * select, double * dr, double * di, double * z, int ldz,
                   double sigmar, double sigmai, double * workev, const char * bmat, int n, const char * which, int nev,
                   double tol, double * resid, int ncv, double * v, int ldv, int * iparam, int * ipntr, double * workd,
                   double * workl, int lworkl, int & info) {
    dneupd_(rvec,howmny,select,dr,di,z,ldz,sigmar,sigmai,workev,
            bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info);
  }
  inline void eupd(int rvec, const char * howmny, const int * select, dcomplex * d, dcomplex * z, int ldz, dcomplex sigma,
                   dcomplex * workev, const char * bmat, int n, const char * which, int nev, double tol, dcomplex * resid,
                   int ncv, dcomplex * v, int ldv, int * iparam, int * ipntr, dcomplex * workd, dcomplex * workl,
                   int lworkl, double * rwork, int & info) {
    zneupd_(rvec,howmny,select,d,z,ldz,sigma,workev,
            bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,rwork,info);
  }

} // f77

} // ezarpack
