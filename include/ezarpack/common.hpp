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
/// @file ezarpack/common.hpp
/// @brief Declarations of common types and exceptions.
#pragma once

#include <complex>
#include <stdexcept>
#include <string>

namespace ezarpack {

/*! The double precision complex type used in ARPACK-NG calls. */
using dcomplex = std::complex<double>;

/// Kind of square matrix (linear operator) to solve an eigenproblem for.
enum operator_kind {
  Symmetric,  /**< Symmetric real matrix. */
  Asymmetric, /**< General real matrix. */
  Complex     /**< General complex matrix. */
};

/// The ARPACK-NG `*aupd_()` procedures set the output argument `IDO` to one of
/// these values to signal the state of the Reverse Communication Interface
/// (RCI).
enum rci_flag : int {
  Init = 0,         /**< First call of the RCI.                             */
  ApplyOpInit = -1, /**< Compute @f$ \mathbf{y} = \hat O \mathbf{x} @f$
                        (force the starting vector into the range of
                        @f$ \hat O @f$).                                    */
  ApplyOp = 1,      /**< Compute @f$ \mathbf{y} = \hat O \mathbf{x} @f$.    */
  ApplyB = 2,       /**< Compute @f$ \mathbf{y} = \hat B \mathbf{x} @f$.    */
  Shifts = 3,       /**< Compute and return the shifts for the Implicitly
                         Restarted Lanczos/Arnoldi Method.                  */
  Done = 99         /**< Done with the iterations.                          */
};

#ifndef DOXYGEN_IGNORE
#define ARPACK_SOLVER_ERROR(MSG) std::runtime_error("arpack_solver: " MSG)
#endif

/// @brief Exception: Maximum number of Implicitly restarted Arnoldi iterations
/// has been reached.
///
/// This exception can be thrown by `operator()` of `arpack_solver`
/// specializations.
struct maxiter_reached : public std::runtime_error {
  /// Maximum number of IRLM/IRAM iterations allowed.
  int maxiter;
  /// @param maxiter Maximum number of IRLM/IRAM iterations allowed.
  maxiter_reached(int maxiter)
      : ARPACK_SOLVER_ERROR("Maximum number of iterations (" +
                            std::to_string(maxiter) + ") reached"),
        maxiter(maxiter) {}
};

/// @brief Exception: No shifts could be applied during a cycle of
/// the Implicitly restarted Arnoldi iteration. Consider increasing the number
/// of Lanczos/Arnoldi vectors generated at each iteration.
///
/// This exception can be thrown by `operator()` of `arpack_solver`
/// specializations.
struct ncv_insufficient : public std::runtime_error {
  int ncv; /**< Number of Lanczos/Arnoldi vectors to be generated. */
  /// @param ncv Number of Lanczos/Arnoldi vectors to be generated.
  ncv_insufficient(int ncv)
      : ARPACK_SOLVER_ERROR("No shifts could be applied during a cycle "
                            "of the Implicitly restarted Arnoldi iteration. "
                            "Try increasing ncv (currently ncv = " +
                            std::to_string(ncv) + ")"),
        ncv(ncv) {}
};

} // namespace ezarpack
