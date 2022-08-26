#
# This file is part of ezARPACK, an easy-to-use C++ wrapper for
# the ARPACK-NG FORTRAN library.
#
# Copyright (C) 2016-2022 Igor Krivenko <igor.s.krivenko@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

set(EXAMPLES)

# Set the DT_RPATH attribute instead of DT_RUNPATH on executables
add_link_options(-Wl,--disable-new-dtags)

# Raw memory example
set(t "raw.example")
add_raw_executable(${t} ${t}.cpp)
target_link_libraries(${t} PRIVATE ${ARPACK_LIBRARIES})
list(APPEND EXAMPLES ${t})

# Eigen3 example
if(EIGEN3_FOUND)
  set(t "eigen.example")
  add_eigen_executable(${t} ${t}.cpp)
  target_link_libraries(${t} PRIVATE ${ARPACK_LIBRARIES})
  list(APPEND EXAMPLES ${t})
endif(EIGEN3_FOUND)

# Blaze example
if(blaze_FOUND)
  set(t "blaze.example")
  add_blaze_executable(${t} ${t}.cpp)
  target_link_libraries(${t} PRIVATE ${ARPACK_LIBRARIES})
  list(APPEND EXAMPLES ${t})
endif(blaze_FOUND)

# Armadillo example
if(Armadillo_FOUND)
  set(t "armadillo.example")
  add_armadillo_executable(${t} ${t}.cpp)
  target_link_libraries(${t} PRIVATE ${ARPACK_LIBRARIES})
  list(APPEND EXAMPLES ${t})
endif(Armadillo_FOUND)

# uBLAS example
if(Boost_FOUND)
  set(t "ublas.example")
  add_ublas_executable(${t} ${t}.cpp)
  target_link_libraries(${t} PRIVATE ${ARPACK_LIBRARIES})
  list(APPEND EXAMPLES ${t})
endif(Boost_FOUND)

# TRIQS example
if(TRIQS_FOUND)
  set(t "triqs.example")
  add_triqs_executable(${t} ${t}.cpp)
  target_link_libraries(${t} PRIVATE ${ARPACK_LIBRARIES})
  list(APPEND EXAMPLES ${t})
endif(TRIQS_FOUND)

# TRIQS/nda example
if(nda_FOUND)
  set(t "nda.example")
  add_nda_executable(${t} ${t}.cpp)
  target_link_libraries(${t} PRIVATE ${ARPACK_LIBRARIES})
  list(APPEND EXAMPLES ${t})
endif(nda_FOUND)

# xtensor example
if(xtensor_FOUND AND xtensor-blas_FOUND)
  set(t "xtensor.example")
  add_xtensor_executable(${t} ${t}.cpp)
  target_link_libraries(${t} PRIVATE ${ARPACK_LIBRARIES})
  list(APPEND EXAMPLES ${t})
endif(xtensor_FOUND AND xtensor-blas_FOUND)

# MPI examples
if(MPI_FOUND)
  add_subdirectory(mpi)
endif(MPI_FOUND)

if(EXAMPLES)
  message(STATUS "Building examples:")
  foreach(example ${EXAMPLES})
    message(STATUS "  ${example}")
  endforeach(example ${EXAMPLES})
endif(EXAMPLES)
