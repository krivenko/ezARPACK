#
# This file is part of ezARPACK, an easy-to-use C++ wrapper for
# the ARPACK-NG FORTRAN library.
#
# Copyright (C) 2016-2023 Igor Krivenko <igor.s.krivenko@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

set(arpack-ng_REQUIRED_VERSION 3.6.0)

message(STATUS "Detecting ARPACK-NG libraries")

if((NOT arpack-ng_DIR) AND ARPACK_NG_ROOT)
  find_package(arpack-ng CONFIG ${arpack-ng_REQUIRED_VERSION}
                         NAMES arpackng arpack-ng
                         PATHS ${ARPACK_NG_ROOT}/lib/cmake
                               ${ARPACK_NG_ROOT}/lib/cmake/arpack-ng
                         NO_DEFAULT_PATH)
else((NOT arpack-ng_DIR) AND ARPACK_NG_ROOT)
  find_package(arpack-ng CONFIG ${arpack-ng_REQUIRED_VERSION}
                NAMES arpackng arpack-ng)
endif((NOT arpack-ng_DIR) AND ARPACK_NG_ROOT)

if(arpack-ng_FOUND)
  message(STATUS "Found ARPACK-NG version ${arpack-ng_VERSION}")

  # Find LAPACK
  find_package(LAPACK)
  if(LAPACK_FOUND)

    # ARPACK-NG prior to version 3.8.0 sets ${arpack_ng_LIBRARIES}
    # It contains full paths to both ARPACK and PARPACK libraries
    if(arpack_ng_LIBRARIES)

      set(ARPACK_LIBRARIES ${arpack_ng_LIBRARIES} ${LAPACK_LIBRARIES})
      if(MPI_FOUND)
        set(PARPACK_LIBRARIES ${arpack_ng_LIBRARIES} ${LAPACK_LIBRARIES})
      endif(MPI_FOUND)

    # ARPACK-NG 3.8.0 and newer export targets ARPACK::ARPACK and
    # PARPACK::PARPACK. Versions >= 3.9.0 export the latter target only if
    # PARPACK was actually built.
    else(arpack_ng_LIBRARIES)

      set(ARPACK_LIBRARIES ARPACK::ARPACK ${LAPACK_LIBRARIES})
      # ARPACK-NG 3.8.0 sets variable ${libdir}
      if(${arpack-ng_VERSION} VERSION_EQUAL 3.8.0)
        set_target_properties(ARPACK::ARPACK PROPERTIES
                              INTERFACE_LINK_DIRECTORIES ${libdir})
      endif(${arpack-ng_VERSION} VERSION_EQUAL 3.8.0)

      if(MPI_FOUND AND (TARGET PARPACK::PARPACK))
        set(PARPACK_LIBRARIES ARPACK::ARPACK PARPACK::PARPACK
            ${LAPACK_LIBRARIES})

        # ARPACK-NG 3.8.0 sets variable ${libdir}
        if(${arpack-ng_VERSION} VERSION_EQUAL 3.8.0)
          set_target_properties(PARPACK::PARPACK PROPERTIES
                                INTERFACE_LINK_DIRECTORIES ${libdir})
        endif(${arpack-ng_VERSION} VERSION_EQUAL 3.8.0)

      endif(MPI_FOUND AND (TARGET PARPACK::PARPACK))

    endif(arpack_ng_LIBRARIES)

  else(LAPACK_FOUND)
    message(STATUS "LAPACK not found. "
                   "Disabling compilation of examples and unit tests.")
    unset(arpack-ng_FOUND)
  endif(LAPACK_FOUND)

else(arpack-ng_FOUND)
  message(STATUS "ARPACK-NG not found. "
                 "Disabling compilation of examples and unit tests.")
endif(arpack-ng_FOUND)

if(arpack-ng_FOUND)
  message(STATUS "ARPACK_LIBRARIES: ${ARPACK_LIBRARIES}")
  message(STATUS "PARPACK_LIBRARIES: ${PARPACK_LIBRARIES}")
endif(arpack-ng_FOUND)
