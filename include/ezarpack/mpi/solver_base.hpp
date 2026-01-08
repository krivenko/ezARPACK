/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2026 Igor Krivenko
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/
/// @file ezarpack/mpi/solver_base.hpp
/// @brief Definition of ezARPACK's primary class template `mpi::arpack_solver`.
#pragma once

#include <array>
#include <limits.h>
#include <utility>

#include "parpack.hpp"

#include "../storages/base.hpp"

namespace ezarpack {
namespace mpi {

/// @brief Main class template of ezARPACK's API (MPI version).
///
/// Instances of this class are used to solve all kinds of eigenproblems
/// supported by ezARPACK in MPI parallelized mode.
/// @tparam OpKind Kind of eigenproblem to be solved.
/// @tparam Backend Tag type specifying what *storage backend* (matrix/vector
/// algebra library) must be used by `arpack_solver`. The storage backend
/// determines types of internally stored data arrays and input/output view
/// objects exposed by methods of the class.
template<operator_kind OpKind, typename Backend> class arpack_solver {};

} // namespace mpi
} // namespace ezarpack
