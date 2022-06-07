#
# This file is part of ezARPACK, an easy-to-use C++ wrapper for
# the ARPACK-NG FORTRAN library.
#
# Copyright (C) 2016-2022 Igor Krivenko <igor.s.krivenko@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# Build Catch2 main object file
add_library(catch2 STATIC catch2/catch2-main.cpp)
set_property(TARGET catch2 PROPERTY CXX_STANDARD 11)

if(MPI_FOUND)
  macro(add_mpi_test name)
    foreach(NP ${ARGN})
      add_test(NAME "${name}.np${NP}"
               COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NP}
                       ${MPIEXEC_PREFLAGS} ${name} ${MPIEXEC_POSTFLAGS})
    endforeach()
  endmacro()
endif()

# Tests of raw memory storage backend
add_subdirectory(raw)

# Tests of Eigen3 storage backend
if(Eigen3_FOUND)
  add_subdirectory(eigen)
endif(Eigen3_FOUND)

# Tests of Blaze storage backend
if(blaze_FOUND)
  add_subdirectory(blaze)
endif(blaze_FOUND)

# Tests of Armadillo storage backend
if(Armadillo_FOUND)
  add_subdirectory(armadillo)
endif(Armadillo_FOUND)

# Tests of uBLAS storage backend
if(Boost_FOUND)
  add_subdirectory(ublas)
endif(Boost_FOUND)

# Tests of TRIQS storage backend
if(TRIQS_FOUND)
  add_subdirectory(triqs)
endif(TRIQS_FOUND)

# Tests of nda storage backend
if(nda_FOUND)
  add_subdirectory(nda)
endif(nda_FOUND)

# Tests of xtensor storage backend
if(xtensor_FOUND AND xtensor-blas_FOUND)
  add_subdirectory(xtensor)
endif(xtensor_FOUND AND xtensor-blas_FOUND)
