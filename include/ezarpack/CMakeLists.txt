#
# This file is part of ezARPACK, an easy-to-use C++ wrapper for
# the ARPACK-NG FORTRAN library.
#
# Copyright (C) 2016-2018 Igor Krivenko <igor.s.krivenko@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# Install C++ headers
install(FILES arpack.hpp
        arpack_worker.hpp
        worker_asymmetric.hpp
        worker_complex.hpp
        worker_symmetric.hpp
        DESTINATION include/ezarpack)

# TODO
# Generate and install version.hpp
#configure_file(version.hpp.in version.hpp)
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.hpp DESTINATION include/triqs/arrays/arpack)
