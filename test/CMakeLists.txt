#
# This file is part of ezARPACK, an easy-to-use C++ wrapper for
# the ARPACK-NG FORTRAN library.
#
# Copyright (C) 2016-2025 Igor Krivenko
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# Build Catch2 main object file
add_library(catch2 STATIC catch2/catch2-main.cpp)
set_property(TARGET catch2 PROPERTY CXX_STANDARD 11)

# Optionally use LLVM sanitizers
option(ASAN "Compile unit tests with LLVM Address Sanitizer" OFF)
option(UBSAN "Compile unit tests with LLVM Undefined Behavior Sanitizer" OFF)
mark_as_advanced(ASAN)
mark_as_advanced(UBSAN)

include_directories(${CMAKE_SOURCE_DIR}/test)

# Set the DT_RPATH attribute instead of DT_RUNPATH on executables
add_link_options(-Wl,--disable-new-dtags)

if(ASAN)
  if(NOT TARGET asan)
    find_package(sanitizer REQUIRED "asan")
  endif()
  link_libraries(asan)
endif()
if(UBSAN)
  if(NOT TARGET ubsan)
    find_package(sanitizer REQUIRED "ubsan")
  endif()
  link_libraries(ubsan)
endif()

# MPI unit tests
if(MPI_FOUND)
  # Build Catch2 object file with a custom main() that initializes and
  # finalizes MPI.
  add_library(catch2_mpi STATIC catch2/catch2-main-mpi.cpp)
  set_property(TARGET catch2_mpi PROPERTY CXX_STANDARD 11)
  target_include_directories(catch2_mpi PUBLIC ${MPI_CXX_INCLUDE_PATH})
  target_link_libraries(catch2_mpi PUBLIC ${MPI_CXX_LIBRARIES})

  # Add an MPI-enabled test
  macro(add_mpi_test name)
    foreach(NP ${ARGN})
      add_test(NAME "${name}.np${NP}"
               COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NP}
                       ${MPIEXEC_PREFLAGS} ${name} ${MPIEXEC_POSTFLAGS})
    endforeach(NP ${ARGN})
  endmacro(add_mpi_test name)
endif(MPI_FOUND)

# arpack_solver test variants corresponding to different OpKind
set(OPERATOR_KINDS symmetric asymmetric complex)
# mpi::arpack_solver test variants corresponding to different OpKind
set(OPERATOR_KINDS_MPI symmetric asymmetric)

# General complex matrix eigensolvers of PARPACK are known to be unstable in
# ARPACK-NG versions prior to 3.8.0. For this reason, we conditionally disable
# the MPI 'complex' unit tests.
#
# C.f. https://github.com/opencollab/arpack-ng/pull/245
if(${arpack-ng_VERSION} VERSION_GREATER_EQUAL 3.8.0)
  set(OPERATOR_KINDS_MPI ${OPERATOR_KINDS_MPI} complex)
else(${arpack-ng_VERSION} VERSION_GREATER_EQUAL 3.8.0)
  message(STATUS "ARPACK-NG is too old, disabling complex PARPACK tests")
endif(${arpack-ng_VERSION} VERSION_GREATER_EQUAL 3.8.0)

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
