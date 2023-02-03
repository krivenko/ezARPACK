#
# This file is part of ezARPACK, an easy-to-use C++ wrapper for
# the ARPACK-NG FORTRAN library.
#
# Copyright (C) 2016-2023 Igor Krivenko <igor.s.krivenko@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

include_directories(${MPI_CXX_INCLUDE_PATH})
link_libraries(${MPI_CXX_LIBRARIES})

# Raw memory example
set(t "raw.mpi.example")
add_raw_executable(${t} ${t}.cpp)
target_link_libraries(${t} PRIVATE ${PARPACK_LIBRARIES})
list(APPEND EXAMPLES ${t})

# Eigen3 example
if(EIGEN3_FOUND)
  set(t "eigen.mpi.example")
  add_eigen_executable(${t} ${t}.cpp)
  target_link_libraries(${t} PRIVATE ${PARPACK_LIBRARIES})
  list(APPEND EXAMPLES ${t})
endif(EIGEN3_FOUND)

# Blaze example
if(blaze_FOUND)
  set(t "blaze.mpi.example")
  add_blaze_executable(${t} ${t}.cpp)
  target_link_libraries(${t} PRIVATE ${PARPACK_LIBRARIES})
  list(APPEND EXAMPLES ${t})
endif(blaze_FOUND)

# Armadillo example
if(Armadillo_FOUND)
  set(t "armadillo.mpi.example")
  add_armadillo_executable(${t} ${t}.cpp)
  target_link_libraries(${t} PRIVATE ${PARPACK_LIBRARIES})
  list(APPEND EXAMPLES ${t})
endif(Armadillo_FOUND)

# uBLAS example
if(Boost_FOUND)
  set(t "ublas.mpi.example")
  add_ublas_executable(${t} ${t}.cpp)
  target_link_libraries(${t} PRIVATE ${PARPACK_LIBRARIES})
  list(APPEND EXAMPLES ${t})
endif(Boost_FOUND)

# TRIQS example
if(TRIQS_FOUND)
  set(t "triqs.mpi.example")
  add_triqs_executable(${t} ${t}.cpp)
  target_link_libraries(${t} PRIVATE ${PARPACK_LIBRARIES})
  list(APPEND EXAMPLES ${t})
endif(TRIQS_FOUND)

# TRIQS/nda example
if(nda_FOUND)
  set(t "nda.mpi.example")
  add_nda_executable(${t} ${t}.cpp)
  target_link_libraries(${t} PRIVATE ${PARPACK_LIBRARIES})
  list(APPEND EXAMPLES ${t})
endif(nda_FOUND)

# xtensor example
if(xtensor_FOUND AND xtensor-blas_FOUND)
  set(t "xtensor.mpi.example")
  add_xtensor_executable(${t} ${t}.cpp)
  target_link_libraries(${t} PRIVATE ${PARPACK_LIBRARIES})
  list(APPEND EXAMPLES ${t})
endif(xtensor_FOUND AND xtensor-blas_FOUND)
