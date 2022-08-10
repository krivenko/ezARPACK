#
# This file is part of ezARPACK, an easy-to-use C++ wrapper for
# the ARPACK-NG FORTRAN library.
#
# Copyright (C) 2016-2022 Igor Krivenko <igor.s.krivenko@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

macro(add_xtensor_executable name source)
  add_executable(${name} ${source})
  target_include_directories(${name} PRIVATE ${CMAKE_SOURCE_DIR}/test
                                             ${xtensor_INCLUDE_DIRS}
                                             ${xtensor_blas_INCLUDE_DIR})
  target_link_libraries(${name} PRIVATE ezarpack ${arpack_ng_LIBRARIES})
  target_compile_definitions(${name} PRIVATE -DXTENSOR_USE_FLENS_BLAS)
  set_property(TARGET ${name} PROPERTY CXX_STANDARD 14)
endmacro()

foreach(t symmetric asymmetric complex)
  set(s ${CMAKE_CURRENT_SOURCE_DIR}/${t}.cpp)
  set(t xtensor.solver_${t})
  add_xtensor_executable(${t} ${s})
  target_link_libraries(${t} PRIVATE catch2)
  add_test(NAME ${t} COMMAND ${t})
endforeach()

# MPI tests
if(MPI_FOUND)
  foreach(t symmetric asymmetric complex)
    set(s ${CMAKE_CURRENT_SOURCE_DIR}/mpi/${t}.cpp)
    set(t xtensor.solver_${t}.mpi)
    add_xtensor_executable(${t} ${s})
    target_link_libraries(${t} PRIVATE catch2_mpi)
    add_mpi_test(${t} 1 2 3 4)
  endforeach()
endif()
