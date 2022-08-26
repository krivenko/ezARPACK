#
# This file is part of ezARPACK, an easy-to-use C++ wrapper for
# the ARPACK-NG FORTRAN library.
#
# Copyright (C) 2016-2022 Igor Krivenko <igor.s.krivenko@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

foreach(t symmetric asymmetric complex)
  set(s ${t}.cpp)
  set(t nda.solver_${t})
  add_nda_executable(${t} ${s})
  target_link_libraries(${t} PRIVATE catch2 ${ARPACK_LIBRARIES})
  add_test(NAME ${t} COMMAND ${t})
endforeach()

# MPI tests
if(MPI_FOUND)
  foreach(t symmetric asymmetric complex)
    set(s mpi/${t}.cpp)
    set(t nda.solver_${t}.mpi)
    add_nda_executable(${t} ${s})
    target_link_libraries(${t} PRIVATE catch2_mpi ${PARPACK_LIBRARIES})
    add_mpi_test(${t} 1 2 3 4)
  endforeach()
endif()
