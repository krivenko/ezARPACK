#
# This file is part of ezARPACK, an easy-to-use C++ wrapper for
# the ARPACK-NG FORTRAN library.
#
# Copyright (C) 2016-2020 Igor Krivenko <igor.s.krivenko@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

foreach(t symmetric asymmetric complex)
  set(s ${CMAKE_CURRENT_SOURCE_DIR}/${t}.cpp)
  set(t triqs.solver_${t})
  add_executable(${t} ${s})
  # FIXME: Setting CXX_STANDARD to 17 would require CMake >= 3.8,
  # which is too new to ask for.
  target_compile_options(${t} PRIVATE "-std=c++17")
  target_include_directories(${t} PRIVATE ${CMAKE_SOURCE_DIR}/test)
  target_link_libraries(${t} PRIVATE ezarpack catch2 triqs)
  triqs_set_rpath_for_target(${t})
  add_test(NAME ${t} COMMAND ${t})
endforeach()
