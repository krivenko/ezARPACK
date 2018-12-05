#
# This file is part of ezARPACK, an easy-to-use C++ wrapper for
# the ARPACK-NG FORTRAN library.
#
# Copyright (C) 2016-2018 Igor Krivenko <igor.s.krivenko@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

set(EXAMPLES)

# TRIQS example
if(TRIQS_FOUND)
  set(t "triqs.example")
  add_executable(${t} triqs.example.cpp)
  target_compile_options(${t} PRIVATE "-std=c++17")
  target_link_libraries(${t} triqs ${arpack_ng_LIBRARIES})
  list(APPEND EXAMPLES ${t})
endif(TRIQS_FOUND)

if(EXAMPLES)
  message(STATUS "Building examples:")
  foreach(example ${EXAMPLES})
    message(STATUS "  ${example}")
  endforeach(example ${EXAMPLES})
endif(EXAMPLES)
