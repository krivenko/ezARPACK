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

#include <array>
#include <iostream>
#include <limits.h>
#include <stdexcept>
#include <string>

#include "arpack.hpp"

#include "storages/base.hpp"

namespace ezarpack {

enum operator_kind { Symmetric, Asymmetric, Complex };

template<operator_kind OpKind, typename Backend> class arpack_worker;

// Exceptions

#define ARPACK_WORKER_ERROR(MSG) std::runtime_error("arpack_worker: " MSG)

struct maxiter_reached : public std::runtime_error {
  int maxiter;
  maxiter_reached(int maxiter)
      : ARPACK_WORKER_ERROR("Maximum number of iterations (" +
                            std::to_string(maxiter) + ") reached"),
        maxiter(maxiter) {}
};

struct ncv_insufficient : public std::runtime_error {
  int ncv;
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
