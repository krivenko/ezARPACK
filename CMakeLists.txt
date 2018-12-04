#
# This file is part of ezARPACK, an easy-to-use C++ wrapper for
# the ARPACK-NG FORTRAN library.
#
# Copyright (C) 2016-2018 Igor Krivenko <igor.s.krivenko@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# Start configuration
cmake_minimum_required(VERSION 3.0.2 FATAL_ERROR)

# Project name and version
set(EZARPACK_VERSION 0.1)
project(ezARPACK VERSION ${EZARPACK_VERSION} LANGUAGES CXX)
message(STATUS "Configuring ${PROJECT_NAME} version ${PROJECT_VERSION}")

# Prohibit in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds are not allowed. Please make a separate build directory.")
endif(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})

# find_package() should take into account <PackageName>_ROOT
# CMake variables and environment variables
cmake_policy(SET CMP0074 NEW)

# CMake options
option(Tests "Build unit tests" ON)
option(Examples "Build examples" ON)

# Are we building any components?
if(Tests OR Examples)

  # Set default build type to Release
  if(NOT DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE)
  endif(NOT DEFINED CMAKE_BUILD_TYPE)

  message(STATUS "Detecting ARPACK-NG headers and libraries")
  set(arpack-ng_DIR ${ARPACK_NG_ROOT}/lib/cmake)
  find_package(arpack-ng 3.6.0 REQUIRED)
  message(STATUS "  arpack_ng_INCLUDE_DIRS: ${arpack_ng_INCLUDE_DIRS}")
  message(STATUS "  arpack_ng_LIBRARIES: ${arpack_ng_LIBRARIES}")

  # Find arpack.hpp
  find_file(ARPACK_HPP "arpack.hpp" "${arpack_ng_INCLUDE_DIRS}")
  if(NOT ARPACK_HPP)
    message(FATAL_ERROR "Could not locate arpack.hpp in ${arpack_ng_INCLUDE_DIRS}.\n"
                        "Did you build ARPACK-NG with C++ bindings enabled?")
  endif(NOT ARPACK_HPP)

  # Add ARPACK-NG and ezARPACK include directories to all targets being built
  include_directories(${CMAKE_SOURCE_DIR}/include ${arpack_ng_INCLUDE_DIRS})

  # Find all storage backends
  include(cmake/find_backends.cmake)
endif()

include_directories(include)

# Install C++ headers
add_subdirectory(include/ezarpack)

# TODO: Install CMake configuration files for ezARPACK

# Build unit tests
if(Tests)
  # FIXME: local build of Google Test conflicts with TRIQS (issue #654)
  #include(cmake/add_gtest.cmake)
  enable_testing()
  add_subdirectory(test)
endif(Tests)

# Build examples
if(Examples)
  add_subdirectory(example)
endif(Examples)
