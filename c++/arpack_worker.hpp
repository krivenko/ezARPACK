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

#include <string>
#include <array>
#include <type_traits>
#include <triqs/arrays/vector.hpp>
#include <triqs/arrays/matrix.hpp>
#include <triqs/utility/c14.hpp>
#include <triqs/utility/exceptions.hpp>

#include "arpack.hpp"

namespace triqs { namespace arrays { namespace arpack {

enum operator_kind {Symmetric, Asymmetric, Complex};

template<operator_kind OpKind> class arpack_worker;

// Reverse Communuication Interface Flag
enum rci_flag : int {Init = 0, ApplyOpInit = -1, ApplyOp = 1, ApplyB = 2, Shifts = 3, Done = 99};

// Exceptions
struct maxiter_reached : public triqs::runtime_error {
 int maxiter;
 maxiter_reached(int maxiter) : maxiter(maxiter) {
  *this << "arpack_worker: maximum number of iterations ("  << maxiter << ") taken.";
 }
 maxiter_reached(maxiter_reached const& ) = default;
 maxiter_reached& operator=(maxiter_reached const&) = default;
};

struct ncv_insufficient : public triqs::runtime_error {
 int ncv;
 ncv_insufficient(int ncv) : ncv(ncv) {
  *this << "arpack_worker: No shifts could be applied during a cycle "
           "of the Implicitly restarted Arnoldi iteration. "
           "Try increasing ncv (currently ncv = " << ncv << ").";
 }
 ncv_insufficient(ncv_insufficient const& ) = default;
 ncv_insufficient& operator=(ncv_insufficient const&) = default;
};

}}}

#include "worker_symmetric.hpp"
#include "worker_asymmetric.hpp"
#include "worker_complex.hpp"
