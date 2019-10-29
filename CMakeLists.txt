#
# This file is part of ezARPACK, an easy-to-use C++ wrapper for
# the ARPACK-NG FORTRAN library.
#
# Copyright (C) 2016-2020 Igor Krivenko <igor.s.krivenko@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# Start configuration
cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

# Project name and version
set(EZARPACK_VERSION 0.8)
project(ezARPACK VERSION ${EZARPACK_VERSION} LANGUAGES CXX)
message(STATUS "Configuring ${PROJECT_NAME} version ${PROJECT_VERSION}")

# Prohibit in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds are not allowed."
                      "Please make a separate build directory.")
endif(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})

# find_package() should take into account <PackageName>_ROOT
# CMake variables and environment variables. Requires CMake 3.12 or newer.
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif(POLICY CMP0074)

# CMake options
option(Tests "Build unit tests" ON)
option(Examples "Build examples" ON)
option(Documentation "Build documentation" OFF)

# Detect ARPACK-NG as an optional dependency
message(STATUS "Detecting ARPACK-NG libraries")
set(arpack-ng_DIR ${ARPACK_NG_ROOT}/lib/cmake)
find_package(arpack-ng CONFIG 3.6.0)
if(arpack-ng_FOUND)
  message(STATUS "Found ARPACK-NG version ${arpack-ng_VERSION}")
  message(STATUS "  arpack_ng_LIBRARIES: ${arpack_ng_LIBRARIES}")
else(arpack-ng_FOUND)
  message(STATUS "ARPACK-NG not found."
                 "Disabling compilation of examples and unit tests.")
endif(arpack-ng_FOUND)

# Are we building any components?
if(arpack-ng_FOUND AND (Tests OR Examples))

  # Set default build type to Release
  if(NOT DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE)
  endif(NOT DEFINED CMAKE_BUILD_TYPE)

  # Add ezARPACK include directory to all targets that require building
  include_directories(${CMAKE_SOURCE_DIR}/include)

  # Find all storage backends
  include(cmake/FindBackends.cmake)
endif()

# Install C++ headers
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include
        DESTINATION ${CMAKE_INSTALL_PREFIX})

# Define interface library target
add_library(ezarpack INTERFACE)
target_include_directories(ezarpack INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

# Optionally add the ARPACK-NG libraries to our interface library
if(arpack-ng_FOUND)
  target_link_libraries(ezarpack INTERFACE ${arpack_ng_LIBRARIES})
endif(arpack-ng_FOUND)

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
        DESTINATION lib/cmake)

# Build documentation
if(Documentation)
  message(STATUS "Building documentation")

  # Detect Doxygen, Sphinx and Breathe
  list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
  find_package(Doxygen REQUIRED)
  find_package(PythonInterp 3 REQUIRED)
  find_package(Sphinx 2.0.0 REQUIRED)
  find_package(Breathe 4.13.0 REQUIRED)

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
