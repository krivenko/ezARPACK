#
# This file is part of ezARPACK, an easy-to-use C++ wrapper for
# the ARPACK-NG FORTRAN library.
#
# Copyright (C) 2016-2020 Igor Krivenko <igor.s.krivenko@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

message(STATUS "Detecting matrix algebra backends")

# Find Eigen3
if(NOT POLICY CMP0074)
  set(Eigen3_DIR ${Eigen3_ROOT}/share/eigen3/cmake)
endif(NOT POLICY CMP0074)
find_package(Eigen3 CONFIG)
if(Eigen3_FOUND)
  if(Eigen3_VERSION)
    message(STATUS "Found Eigen3 version ${Eigen3_VERSION}")
    if(NOT ${Eigen3_VERSION} VERSION_LESS "3.3.0")
      add_compile_definitions(EIGEN_CAN_MIX_REAL_COMPLEX_EXPR)
    endif(NOT ${Eigen3_VERSION} VERSION_LESS "3.3.0")
  else(Eigen3_VERSION)
    message(STATUS "Found Eigen3 version ${EIGEN3_VERSION_STRING}")
    if(NOT ${EIGEN3_VERSION_STRING} VERSION_LESS "3.3.0")
      add_compile_definitions(EIGEN_CAN_MIX_REAL_COMPLEX_EXPR)
    endif(NOT ${EIGEN3_VERSION_STRING} VERSION_LESS "3.3.0")
  endif(Eigen3_VERSION)
endif(Eigen3_FOUND)

# Find Blaze
if(NOT POLICY CMP0074)
  set(blaze_DIR ${blaze_ROOT}/share/blaze/cmake)
endif(NOT POLICY CMP0074)
find_package(blaze 3.0 QUIET CONFIG)
if(blaze_FOUND)
  message(STATUS "Found Blaze version ${blaze_VERSION}")
endif(blaze_FOUND)

# Armadillo
if(NOT POLICY CMP0074)
  set(Armadillo_DIR ${Armadillo_ROOT}/share/Armadillo/CMake)
endif(NOT POLICY CMP0074)
find_package(Armadillo QUIET CONFIG)
if(NOT Armadillo_FOUND)
  if(Armadillo_ROOT)
    find_library(ARMADILLO_LIBRARY NAMES armadillo
                 PATHS ${Armadillo_ROOT}/lib
                 NO_DEFAULT_PATH)
  else(Armadillo_ROOT)
    find_library(ARMADILLO_LIBRARY NAMES armadillo)
  endif(Armadillo_ROOT)

  if(ARMADILLO_LIBRARY)
    set(Armadillo_FOUND TRUE)
    get_filename_component(ARMADILLO_PREFIX ${ARMADILLO_LIBRARY} DIRECTORY)
    add_library(armadillo INTERFACE)
    target_include_directories(armadillo INTERFACE ${ARMADILLO_PREFIX}/include)
    target_link_libraries(armadillo INTERFACE ${ARMADILLO_LIBRARY})
  endif(ARMADILLO_LIBRARY)
endif(NOT Armadillo_FOUND)
if(Armadillo_FOUND)
  message(STATUS "Found Armadillo")
endif(Armadillo_FOUND)

# Boost uBLAS
find_package(Boost 1.58)

# Find TRIQS
if(NOT POLICY CMP0074)
  set(Cpp2Py_DIR ${TRIQS_ROOT}/lib/cmake/cpp2py)
  set(TRIQS_DIR ${TRIQS_ROOT}/lib/cmake/triqs)
endif(NOT POLICY CMP0074)
find_package(Cpp2Py CONFIG)
find_package(TRIQS CONFIG)

# Find xtensor
if(NOT POLICY CMP0074)
  set(xtl_DIR ${xtensor_ROOT}/lib/cmake/xtl)
  set(xtensor_DIR ${xtensor_ROOT}/lib/cmake/xtensor)
  set(xtensor-blas_DIR ${xtensor-blas_ROOT}/lib/cmake/xtensor-blas)
endif(NOT POLICY CMP0074)
find_package(xtensor CONFIG 0.20)
find_package(xtensor-blas CONFIG 0.16)
if(xtensor_FOUND AND xtensor-blas_FOUND)
  message(STATUS "Found xtensor version ${xtensor_VERSION}")
  message(STATUS "Found xtensor-blas version ${xtensor-blas_VERSION}")
endif(xtensor_FOUND AND xtensor-blas_FOUND)
