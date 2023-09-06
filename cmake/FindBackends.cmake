#
# This file is part of ezARPACK, an easy-to-use C++ wrapper for
# the ARPACK-NG FORTRAN library.
#
# Copyright (C) 2016-2023 Igor Krivenko <igor.s.krivenko@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

message(STATUS "Detecting matrix algebra backends")

macro(add_raw_executable name source)
  add_executable(${name} ${source})
  target_link_libraries(${name} PRIVATE ezarpack)
  set_property(TARGET ${name} PROPERTY CXX_STANDARD 11)
endmacro()

# Find Eigen3
find_package(Eigen3 CONFIG)
if(Eigen3_FOUND)
  if(Eigen3_VERSION)
    message(STATUS "Found Eigen3 version ${Eigen3_VERSION}")
    if(NOT ${Eigen3_VERSION} VERSION_LESS "3.3.0")
      set(EIGEN_CAN_MIX_REAL_COMPLEX_EXPR TRUE)
    endif(NOT ${Eigen3_VERSION} VERSION_LESS "3.3.0")
  else(Eigen3_VERSION)
    message(STATUS "Found Eigen3 version ${EIGEN3_VERSION_STRING}")
    if(NOT ${EIGEN3_VERSION_STRING} VERSION_LESS "3.3.0")
      set(EIGEN_CAN_MIX_REAL_COMPLEX_EXPR TRUE)
    endif(NOT ${EIGEN3_VERSION_STRING} VERSION_LESS "3.3.0")
  endif(Eigen3_VERSION)

  macro(add_eigen_executable name source)
    add_executable(${name} ${source})
    target_link_libraries(${name} PRIVATE ezarpack)
    if(TARGET Eigen3::Eigen)
      target_link_libraries(${name} PRIVATE Eigen3::Eigen)
    else(TARGET Eigen3::Eigen)
      target_compile_definitions(${name} PRIVATE ${EIGEN3_DEFINITIONS})
      target_include_directories(${name} PRIVATE ${EIGEN3_INCLUDE_DIRS})
    endif(TARGET Eigen3::Eigen)
    set_property(TARGET ${name} PROPERTY CXX_STANDARD 11)
    if(EIGEN_CAN_MIX_REAL_COMPLEX_EXPR)
      target_compile_definitions(${name} PRIVATE
                                 EIGEN_CAN_MIX_REAL_COMPLEX_EXPR)
    endif(EIGEN_CAN_MIX_REAL_COMPLEX_EXPR)
  endmacro()
endif(Eigen3_FOUND)

# Find Blaze
find_package(blaze 3.0 QUIET CONFIG)
if(blaze_FOUND)
  message(STATUS "Found Blaze version ${blaze_VERSION}")

  macro(add_blaze_executable name source)
    add_executable(${name} ${source})
    target_link_libraries(${name} PRIVATE ezarpack blaze::blaze)
    set_property(TARGET ${name} PROPERTY CXX_STANDARD 14)
  endmacro()
endif(blaze_FOUND)

# Armadillo
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

  macro(add_armadillo_executable name source)
    add_executable(${name} ${source})
    target_link_libraries(${name} PRIVATE ezarpack armadillo)
    set_property(TARGET ${name} PROPERTY CXX_STANDARD 11)
  endmacro()
endif(Armadillo_FOUND)

# Boost uBLAS
find_package(Boost 1.58)
if(Boost_FOUND)
  macro(add_ublas_executable name source)
    add_executable(${name} ${source})
    target_include_directories(${name} PRIVATE ${Boost_INCLUDE_DIRS})
    target_link_libraries(${name} PRIVATE ezarpack)
    set_property(TARGET ${t} PROPERTY CXX_STANDARD 11)
  endmacro()
endif(Boost_FOUND)

# Find TRIQS
find_package(Cpp2Py CONFIG)
find_package(TRIQS CONFIG)
if(TRIQS_FOUND)
  if(TRIQS_VERSION VERSION_LESS "3.1")
    macro(add_triqs_executable name source)
      add_executable(${name} ${source})
      target_link_libraries(${name} PRIVATE ezarpack triqs)
      triqs_set_rpath_for_target(${name})
    endmacro()
  else()
    unset(TRIQS_FOUND)
    message(WARNING "TRIQS ${TRIQS_VERSION} (>3.0.x) now uses nda arrays. "
                    "Skipping building TRIQS unit tests and examples.")
  endif()
endif(TRIQS_FOUND)

# Find TRIQS nda
find_package(nda 1.1.0 CONFIG)
if(nda_FOUND)
  macro(add_nda_executable name source)
    add_executable(${name} ${source})
    target_link_libraries(${name} PRIVATE ezarpack nda::nda_c)
    if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      set_target_properties(${name} PROPERTIES
                            INSTALL_NAME_DIR "${nda_ROOT}/lib")
      set_target_properties(${name} PROPERTIES
                            INSTALL_RPATH    "${nda_ROOT}/lib")
      set_target_properties(${name} PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
    endif(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(${name} PROPERTIES BUILD_WITH_INSTALL_RPATH FALSE)
    set_target_properties(${name} PROPERTIES SKIP_BUILD_RPATH FALSE)
    set_target_properties(${name} PROPERTIES SKIP_INSTALL_RPATH FALSE)
  endmacro()
endif(nda_FOUND)

# Find xtensor
find_package(xtensor CONFIG 0.20)
find_package(xtensor-blas CONFIG 0.16)
if(xtensor_FOUND AND xtensor-blas_FOUND)
  message(STATUS "Found xtensor version ${xtensor_VERSION}")
  message(STATUS "Found xtensor-blas version ${xtensor-blas_VERSION}")

  macro(add_xtensor_executable name source)
    add_executable(${name} ${source})
    target_include_directories(${name} PRIVATE ${xtensor_INCLUDE_DIRS}
                                               ${xtensor_blas_INCLUDE_DIR})
    target_link_libraries(${name} PRIVATE ezarpack)
    target_compile_definitions(${name} PRIVATE -DXTENSOR_USE_FLENS_BLAS)
    set_property(TARGET ${name} PROPERTY CXX_STANDARD 14)
  endmacro()
endif(xtensor_FOUND AND xtensor-blas_FOUND)
