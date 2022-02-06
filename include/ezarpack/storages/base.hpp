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
#pragma once

namespace ezarpack {

template<typename T> constexpr bool unsupportedStorageBackend() {
  return false;
}

/// Specializations of this traits structure describe the way
/// @ref ezarpack::arpack_solver class makes use of a particular vector/matrix
/// algebra library. The primary template cannot be instantiated, unless there
/// is a valid specialization for the requested tag type `Backend`.
/// @tparam Backend A tag type denoting a storage backend (vector/matrix algebra
/// library).
template<typename Backend> struct storage_traits {
  static_assert(unsupportedStorageBackend<Backend>(),
                "Storage backend is unsupported");
};

} // namespace ezarpack
