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
struct maxiter_reached : public std::exception {
 int maxiter;
 std::string msg;
 maxiter_reached(int maxiter) : maxiter(maxiter) {
  msg = "arpack_worker: maximum number of iterations (" + std::to_string(maxiter) + ") taken.";
 }
 maxiter_reached(maxiter_reached const& ) = default;
 maxiter_reached& operator=(maxiter_reached const&) = default;
 const char* what() const noexcept { return msg.c_str(); }
};

struct ncv_insufficient : public std::exception {
 int ncv;
 std::string msg;
 ncv_insufficient(int ncv) : ncv(ncv) {
  msg = "arpack_worker: No shifts could be applied during a cycle "
        "of the Implicitly restarted Arnoldi iteration. "
        "Try increasing ncv (currently ncv = " + std::to_string(ncv) + ").";
 }
 ncv_insufficient(ncv_insufficient const& ) = default;
 ncv_insufficient& operator=(ncv_insufficient const&) = default;
 const char* what() const noexcept { return msg.c_str(); }
};

}}}

#include "worker_symmetric.hpp"
#include "worker_asymmetric.hpp"
#include "worker_complex.hpp"
