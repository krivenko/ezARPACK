#
# This file is part of ezARPACK, an easy-to-use C++ wrapper for
# the ARPACK-NG FORTRAN library.
#
# Copyright (C) 2016-2020 Igor Krivenko <igor.s.krivenko@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

if(NOT MATHJAX_DIR)
  find_path(MATHJAX_DIR
    NAMES MathJax.js
    PATHS "/usr/share/javascript/mathjax"
    DOC "Path to MathJax.js script"
  )

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(MathJax
    FOUND_VAR MathJax_FOUND
    REQUIRED_VARS MATHJAX_DIR
    FAIL_MESSAGE "Failed to find MathJax"
  )
endif(NOT MATHJAX_DIR)
