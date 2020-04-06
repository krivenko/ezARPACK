#
# This file is part of ezARPACK, an easy-to-use C++ wrapper for
# the ARPACK-NG FORTRAN library.
#
# Copyright (C) 2016-2020 Igor Krivenko <igor.s.krivenko@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# Build Catch2 main object file
add_library(catch2 STATIC catch2/catch2-main.cpp)
set_property(TARGET catch2 PROPERTY CXX_STANDARD 11)

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

# Tests of xtensor storage backend
if(xtensor_FOUND AND xtensor-blas_FOUND)
  add_subdirectory(xtensor)
endif(xtensor_FOUND AND xtensor-blas_FOUND)
