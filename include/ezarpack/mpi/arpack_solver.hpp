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
/// @file ezarpack/mpi/arpack_solver.hpp
/// @brief This header file includes all specializations of
/// `mpi::arpack_solver`.
#pragma once

// clang-format off
#include "solver_base.hpp"
#include "solver_asymmetric.hpp"
#include "solver_complex.hpp"
#include "solver_symmetric.hpp"
// clang-format on
