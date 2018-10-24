#
# This file is part of ezARPACK, an easy-to-use C++ wrapper for
# the ARPACK-NG FORTRAN library.
#
# Copyright (C) 2016-2018 Igor Krivenko <igor.s.krivenko@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

set(EXAMPLES example)

foreach(t ${EXAMPLES})
    set(s ${CMAKE_CURRENT_SOURCE_DIR}/${t}.cpp)
    add_executable(${t} ${s})
    target_link_libraries(${t} triqs ${ARPACK_LIB})
    triqs_set_rpath_for_target(${t})
endforeach()
