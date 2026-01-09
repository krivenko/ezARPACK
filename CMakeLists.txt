#
# This file is part of ezARPACK, an easy-to-use C++ wrapper for
# the ARPACK-NG FORTRAN library.
#
# Copyright (C) 2016-2026 Igor Krivenko
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# Start configuration
cmake_minimum_required(VERSION 3.18.0 FATAL_ERROR)

# Project name and version
set(EZARPACK_VERSION 1.1)
project(ezARPACK VERSION ${EZARPACK_VERSION} LANGUAGES CXX)
message(STATUS "Configuring ${PROJECT_NAME} version ${PROJECT_VERSION}")

# Prohibit in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds are not allowed."
                      "Please make a separate build directory.")
endif(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})

# CMake options
option(Documentation "Build documentation" OFF)
option(EnableMPI "Build with MPI support (PARPACK)" ON)
option(Examples "Build examples" ON)
option(Tests "Build unit tests" ON)

# Are we building any executables?
if(Tests OR Examples)

  # Detect if linker supports --disable-new-dtags
  include(CheckLinkerFlag)
  check_linker_flag(CXX "LINKER:--disable-new-dtags"
                    LINKER_HAS_DISABLE_NEW_TAGS)

  # Detect an MPI implementation for MPI-tests and/or examples
  if (EnableMPI)
    find_package(MPI 3.0)
    if(MPI_FOUND)
      message(STATUS "MPI_CXX_INCLUDE_PATH: ${MPI_CXX_INCLUDE_PATH}")
      message(STATUS "MPI_CXX_LIBRARIES: ${MPI_CXX_LIBRARIES}")
    endif(MPI_FOUND)
  endif(EnableMPI)

  # Detect ARPACK-NG
  include(cmake/FindARPACKNG.cmake)

  if(arpack-ng_FOUND)
    # Set default build type to Release
    if(NOT DEFINED CMAKE_BUILD_TYPE)
      set(CMAKE_BUILD_TYPE Release CACHE)
    endif(NOT DEFINED CMAKE_BUILD_TYPE)

    # Find all storage backends
    include(cmake/FindBackends.cmake)

    # Add ezARPACK include directory to all targets that require building
    include_directories(${CMAKE_SOURCE_DIR}/include)
  endif(arpack-ng_FOUND)

endif(Tests OR Examples)

# Install C++ headers
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include
        DESTINATION ${CMAKE_INSTALL_PREFIX})

# Define interface library target
add_library(ezarpack INTERFACE)
target_include_directories(ezarpack INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

# Write config version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/ezARPACKConfigVersion.cmake"
  VERSION ${EZARPACK_VERSION}
  COMPATIBILITY AnyNewerVersion
)

# Install CMake files
install(TARGETS ezarpack
  EXPORT ezARPACKTargets
  PUBLIC_HEADER DESTINATION include COMPONENT Development
)
configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/ezARPACKConfig.cmake.in"
    "${PROJECT_BINARY_DIR}/ezARPACKConfig.cmake"
    INSTALL_DESTINATION lib/cmake
)
install(EXPORT ezARPACKTargets DESTINATION lib/cmake)
install(FILES "${PROJECT_BINARY_DIR}/ezARPACKConfigVersion.cmake"
              "${PROJECT_BINARY_DIR}/ezARPACKConfig.cmake"
              "${PROJECT_SOURCE_DIR}/cmake/ezARPACKMacros.cmake"
        DESTINATION lib/cmake)

# Build documentation
if(Documentation)
  message(STATUS "Building documentation")

  # Detect Doxygen, Sphinx, Breathe and MathJax
  list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
  find_package(Doxygen REQUIRED)
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  find_package(Sphinx 2.0.0 REQUIRED)
  find_package(Breathe 4.13.0 REQUIRED)
  find_package(MySTParser REQUIRED)
  find_package(MathJax REQUIRED)

  add_subdirectory(doc)
endif(Documentation)

# Build unit tests
if(arpack-ng_FOUND AND Tests)
  enable_testing()
  add_subdirectory(test)
endif(arpack-ng_FOUND AND Tests)

# Build examples
if(arpack-ng_FOUND AND Examples)
  add_subdirectory(example)
endif(arpack-ng_FOUND AND Examples)
