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

# Raw memory example
add_executable(raw.example raw.example.cpp)
target_link_libraries(raw.example ezarpack ${arpack_ng_LIBRARIES})
set_property(TARGET raw.example PROPERTY CXX_STANDARD 11)
list(APPEND EXAMPLES raw.example)

# Eigen3 example
if(EIGEN3_FOUND)
  set(t "eigen.example")
  add_executable(${t} eigen.example.cpp)
  target_link_libraries(${t} PRIVATE ezarpack ${arpack_ng_LIBRARIES})
  if(TARGET Eigen3::Eigen)
    target_link_libraries(${t} PRIVATE Eigen3::Eigen)
  else(TARGET Eigen3::Eigen)
    target_compile_definitions(${t} PRIVATE ${EIGEN3_DEFINITIONS})
    target_include_directories(${t} PRIVATE ${EIGEN3_INCLUDE_DIRS})
  endif(TARGET Eigen3::Eigen)
  set_property(TARGET ${t} PROPERTY CXX_STANDARD 11)
  list(APPEND EXAMPLES ${t})
endif(EIGEN3_FOUND)

# Blaze example
if(blaze_FOUND)
  set(t "blaze.example")
  add_executable(${t} blaze.example.cpp)
  target_link_libraries(${t} ezarpack blaze::blaze ${arpack_ng_LIBRARIES})
  set_property(TARGET ${t} PROPERTY CXX_STANDARD 11)
  list(APPEND EXAMPLES ${t})
endif(blaze_FOUND)

# Armadillo example
if(Armadillo_FOUND)
  set(t "armadillo.example")
  add_executable(${t} armadillo.example.cpp)
  target_link_libraries(${t} ezarpack armadillo ${arpack_ng_LIBRARIES})
  set_property(TARGET ${t} PROPERTY CXX_STANDARD 11)
  list(APPEND EXAMPLES ${t})
endif(Armadillo_FOUND)

# uBLAS example
if(Boost_FOUND)
  set(t "ublas.example")
  add_executable(${t} ublas.example.cpp)
  target_include_directories(${t} PRIVATE ${Boost_INCLUDE_DIRS})
  target_link_libraries(${t} ezarpack ${arpack_ng_LIBRARIES})
  set_property(TARGET ${t} PROPERTY CXX_STANDARD 11)
  list(APPEND EXAMPLES ${t})
endif(Boost_FOUND)

# TRIQS example
if(TRIQS_FOUND)
  set(t "triqs.example")
  add_executable(${t} triqs.example.cpp)
  # FIXME: Setting CXX_STANDARD to 17 would require CMake >= 3.8,
  # which is too new to ask for.
  target_compile_options(${t} PRIVATE "-std=c++17")
  target_link_libraries(${t} ezarpack triqs ${arpack_ng_LIBRARIES})
  list(APPEND EXAMPLES ${t})
endif(TRIQS_FOUND)

# TRIQS/nda example
if(nda_FOUND)
  set(t "nda.example")
  add_executable(${t} nda.example.cpp)
  # FIXME: Setting CXX_STANDARD to 20 would require CMake >= 3.12,
  # which is too new to ask for.
  target_compile_options(${t} PRIVATE "-std=c++20")
  target_link_libraries(${t} ezarpack nda::nda_c ${arpack_ng_LIBRARIES})
  list(APPEND EXAMPLES ${t})
endif(nda_FOUND)

# xtensor example
if(xtensor_FOUND AND xtensor-blas_FOUND)
  set(t "xtensor.example")
  add_executable(${t} xtensor.example.cpp)
  target_include_directories(${t} PRIVATE ${xtensor_INCLUDE_DIRS})
  target_link_libraries(${t} ezarpack xtensor ${arpack_ng_LIBRARIES})
  set_property(TARGET ${t} PROPERTY CXX_STANDARD 14)
  list(APPEND EXAMPLES ${t})
endif(xtensor_FOUND AND xtensor-blas_FOUND)

if(EXAMPLES)
  message(STATUS "Building examples:")
  foreach(example ${EXAMPLES})
    message(STATUS "  ${example}")
  endforeach(example ${EXAMPLES})
endif(EXAMPLES)
