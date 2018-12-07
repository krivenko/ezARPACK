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

#include <array>
#include <iostream>
#include <limits.h>
#include <string>
#include <exception>

#include "arpack.hpp"

namespace ezarpack {

enum operator_kind {Symmetric, Asymmetric, Complex};

template<operator_kind OpKind, typename Backend> class arpack_worker;

// Reverse Communuication Interface Flag
enum rci_flag : int {Init = 0, ApplyOpInit = -1, ApplyOp = 1, ApplyB = 2, Shifts = 3, Done = 99};

// Exceptions
struct maxiter_reached : public std::runtime_error {
  int maxiter;
  maxiter_reached(int maxiter)
    : std::runtime_error("arpack_worker: maximum number of iterations ("  + std::to_string(maxiter) + ") reached"),
      maxiter(maxiter) {}
};

struct ncv_insufficient : public std::runtime_error {
  int ncv;
  ncv_insufficient(int ncv)
    : std::runtime_error("arpack_worker: No shifts could be applied during a cycle "
                         "of the Implicitly restarted Arnoldi iteration. "
                         "Try increasing ncv (currently ncv = " + std::to_string(ncv) + ")."),
      ncv(ncv) {}
};

}

#include "worker_symmetric.hpp"
#include "worker_asymmetric.hpp"
#include "worker_complex.hpp"
