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
/// @file ezarpack/arpack_worker.hpp
/// @brief Definition of ezARPACK's main class template `arpack_worker` and
/// related types.
///
/// This header file also includes all specializations of `arpack_worker`.
#pragma once

#include <array>
#include <limits.h>
#include <stdexcept>
#include <string>
#include <utility>

#include "arpack.hpp"

#include "storages/base.hpp"

namespace ezarpack {

/// Kind of square matrix (linear operator) to solve an eigenproblem for.
enum operator_kind {
  Symmetric,    /**< Symmetric real matrix. */
  Asymmetric,   /**< General real matrix. */
  Complex       /**< General complex matrix. */
};

/// @class arpack_worker
/// @brief Main class template of ezARPACK's API.
///
/// Instances of this class are used to solve all kinds of eigenproblems
/// supported by ezARPACK.
/// @tparam OpKind Kind of eigenproblem to be solved.
/// @tparam Backend Tag type specifying what *storage backend* (matrix/vector
/// algebra library) must be used by `arpack_worker`. The storage backend
/// determines types of internally stored data arrays and input/output view
/// objects exposed by methods of the class.
template<operator_kind OpKind, typename Backend> class arpack_worker;

#ifndef DOXYGEN_IGNORE
#define ARPACK_WORKER_ERROR(MSG) std::runtime_error("arpack_worker: " MSG)
#endif

/// @brief Exception: Maximum number of Implicitly restarted Arnoldi iterations
/// has been reached.
///
/// This exception can be thrown by `operator()` of `arpack_worker`
/// specializations.
struct maxiter_reached : public std::runtime_error {
  int maxiter; /**< Maximum number of iterations.*/
  maxiter_reached(int maxiter)
      : ARPACK_WORKER_ERROR("Maximum number of iterations (" +
                            std::to_string(maxiter) + ") reached"),
        maxiter(maxiter) {}
};

/// @brief Exception: No shifts could be applied during a cycle of
/// the Implicitly restarted Arnoldi iteration. Consider increasing the number
/// of Lanczos/Arnoldi vectors generated at each iteration.
///
/// This exception can be thrown by `operator()` of `arpack_worker`
/// specializations.
struct ncv_insufficient : public std::runtime_error {
  int ncv; /**< Number of Lanczos/Arnoldi vectors generated at each iteration.*/
  ncv_insufficient(int ncv)
      : ARPACK_WORKER_ERROR("No shifts could be applied during a cycle "
                            "of the Implicitly restarted Arnoldi iteration. "
                            "Try increasing ncv (currently ncv = " +
                            std::to_string(ncv) + ")"),
        ncv(ncv) {}
};

} // namespace ezarpack

#include "worker_asymmetric.hpp"
#include "worker_complex.hpp"
#include "worker_symmetric.hpp"

#undef ARPACK_WORKER_ERROR
