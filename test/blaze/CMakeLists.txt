#
# This file is part of ezARPACK, an easy-to-use C++ wrapper for
# the ARPACK-NG FORTRAN library.
#
# Copyright (C) 2016-2020 Igor Krivenko <igor.s.krivenko@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

find_package(LAPACK)
if(LAPACK_FOUND)
  foreach(t symmetric asymmetric complex)
    set(s ${CMAKE_CURRENT_SOURCE_DIR}/${t}.cpp)
    set(t blaze.worker_${t})
    add_executable(${t} ${s})
    target_include_directories(${t} PRIVATE ${CMAKE_SOURCE_DIR}/test)
    target_link_libraries(${t} PRIVATE ezarpack catch2
                                       blaze::blaze ${LAPACK_LIBRARIES})
    set_property(TARGET ${t} PROPERTY CXX_STANDARD 11)
    add_test(NAME ${t} COMMAND ${t})
  endforeach()
else(LAPACK_FOUND)
  message(STATUS "Cannot build Blaze storage backend tests without LAPACK")
endif(LAPACK_FOUND)
