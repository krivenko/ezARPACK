#
# This file is part of ezARPACK, an easy-to-use C++ wrapper for
# the ARPACK-NG FORTRAN library.
#
# Copyright (C) 2016-2024 Igor Krivenko <igor.s.krivenko@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# This macro finds a working installation of ARPACK-NG while dealing with
# version-to-version differences of ARPACK-NG's CMake interface.
#
# When called as find_arpackng([VERSION] [REQUIRED]), it forwards its optional
# arguments to the underlying `find_package()` calls.
#
# Setting CMake variable `ARPACK_NG_ROOT` will instruct `find_arpackng()` to
# look for ARPACK-NG at a specific installation prefix first.
#
# Upon successful detection of ARPACK-NG, `find_arpackng()` sets two variables
# that can be later passed to `target_link_libraries()`,
# - `ARPACK_LIBRARIES` - list of ARPACK libraries and linker flags
# - `PARPACK_LIBRARIES` - list of Parallel ARPACK libraries and linker flags
macro(find_arpackng)
  # Try ARPACK_NG_ROOT first ...
  set(ARGV_NO_REQUIRED ${ARGV})
  list(REMOVE_ITEM ARGV_NO_REQUIRED "REQUIRED")
  find_package(arpack-ng ${ARGV_NO_REQUIRED} QUIET CONFIG
               NAMES arpackng arpack-ng
               PATHS ${ARPACK_NG_ROOT}/lib/cmake
                     ${ARPACK_NG_ROOT}/lib/cmake/arpack-ng
               NO_DEFAULT_PATH)

  # If failed, try default CMake locations
  if(NOT arpack-ng_FOUND)
    find_package(arpack-ng ${ARGV} CONFIG NAMES arpackng arpack-ng)
  endif(NOT arpack-ng_FOUND)

  if(arpack-ng_FOUND)
    message(STATUS "Found ARPACK-NG version ${arpack-ng_VERSION}: "
                   "${arpack-ng_DIR}")

    # ARPACK-NG 3.8.0 and newer exports CMake targets
    if(TARGET ARPACK::ARPACK)
      set(ARPACK_LIBRARIES ARPACK::ARPACK)

      # Has target PARPACK::PARPACK been also exported by ARPACK-NG?
      if(TARGET PARPACK::PARPACK)
        set(PARPACK_LIBRARIES ARPACK::ARPACK PARPACK::PARPACK)
      endif(TARGET PARPACK::PARPACK)

      # ARPACK-NG 3.8.0 sets variable ${libdir}
      if(${arpack-ng_VERSION} VERSION_EQUAL 3.8.0)
        set_target_properties(ARPACK::ARPACK PROPERTIES
                              INTERFACE_LINK_DIRECTORIES ${libdir})
        if(TARGET PARPACK::PARPACK)
          set_target_properties(PARPACK::PARPACK PROPERTIES
                                INTERFACE_LINK_DIRECTORIES ${libdir})
        endif(TARGET PARPACK::PARPACK)
      endif(${arpack-ng_VERSION} VERSION_EQUAL 3.8.0)

    # ARPACK-NG prior to version 3.8.0 sets ${arpack_ng_LIBRARIES}
    # It contains full paths to both ARPACK and PARPACK libraries
    elseif(arpack_ng_LIBRARIES)
      set(ARPACK_LIBRARIES ${arpack_ng_LIBRARIES} ${LAPACK_LIBRARIES})
      set(PARPACK_LIBRARIES ${arpack_ng_LIBRARIES} ${LAPACK_LIBRARIES})
    else()
      message(WARNING "Could not deduce a list of ARPACK-NG libraries.")
    endif(TARGET ARPACK::ARPACK)

  endif(arpack-ng_FOUND)

endmacro(find_arpackng)
