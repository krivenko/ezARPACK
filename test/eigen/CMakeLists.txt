#
# This file is part of ezARPACK, an easy-to-use C++ wrapper for
# the ARPACK-NG FORTRAN library.
#
# Copyright (C) 2016-2022 Igor Krivenko <igor.s.krivenko@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

macro(add_eigen_executable name source)
  add_executable(${name} ${source})
  target_include_directories(${name} PRIVATE ${CMAKE_SOURCE_DIR}/test)
  target_link_libraries(${name} PRIVATE ezarpack catch2 ${arpack_ng_LIBRARIES})
  if(TARGET Eigen3::Eigen)
    target_link_libraries(${name} PRIVATE Eigen3::Eigen)
  else(TARGET Eigen3::Eigen)
    target_compile_definitions(${name} PRIVATE ${EIGEN3_DEFINITIONS})
    target_include_directories(${name} PRIVATE ${EIGEN3_INCLUDE_DIRS})
  endif(TARGET Eigen3::Eigen)
  set_property(TARGET ${name} PROPERTY CXX_STANDARD 11)
endmacro()

foreach(t symmetric asymmetric complex)
  set(s ${CMAKE_CURRENT_SOURCE_DIR}/${t}.cpp)
  set(t eigen.solver_${t})
  add_eigen_executable(${t} ${s})
  add_test(NAME ${t} COMMAND ${t})
endforeach()

# MPI tests
if(MPI_FOUND)
  foreach(t symmetric asymmetric complex)
    set(s ${CMAKE_CURRENT_SOURCE_DIR}/${t}_mpi.cpp)
    set(t eigen.solver_${t}.mpi)
    add_eigen_executable(${t} ${s})
    add_mpi_test(${t} 1 2 4)
  endforeach()
endif()
