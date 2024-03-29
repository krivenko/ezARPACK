#
# This file is part of ezARPACK, an easy-to-use C++ wrapper for
# the ARPACK-NG FORTRAN library.
#
# Copyright (C) 2016-2024 Igor Krivenko <igor.s.krivenko@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

file(GLOB_RECURSE ALL_HEADERS ${CMAKE_SOURCE_DIR}/include/*.hpp)

set(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR}/include)
set(DOXYGEN_INDEX_FILE ${CMAKE_CURRENT_BINARY_DIR}/doxygen/xml/index.xml)
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/doxygen)
set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

# Generate Doxygen config
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

# Create Doxygen output directory
file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR}/xml)

# Run Doxygen
add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
                   DEPENDS ${ALL_HEADERS}
                   COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                   MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
                   COMMENT "Parsing sources with Doxygen")

add_custom_target(Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})

set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx)
set(SPHINX_CONFIG ${CMAKE_CURRENT_BINARY_DIR}/conf.py)
set(SPHINX_INDEX_FILE ${SPHINX_BUILD}/index.html)
set(BREATHE_PROJECT ${DOXYGEN_OUTPUT_DIR}/xml)

file(GLOB_RECURSE ALL_RST ${CMAKE_CURRENT_SOURCE_DIR}/*.rst)

if(EXISTS "${MATHJAX_DIR}/MathJax.js")        # MathJax 2
  set(MATHJAX_FILE_PATH "${MATHJAX_DIR}/MathJax.js?config=default")
elseif(EXISTS "${MATHJAX_DIR}/tex-chtml.js")  # MathJax 3
  set(MATHJAX_FILE_PATH "${MATHJAX_DIR}/tex-chtml.js")
else()
  message(FATAL_ERROR "Could not find a MathJax component file")
endif()

# Generate Sphinx config
configure_file(conf.py.in conf.py @ONLY)

# Run Sphinx
add_custom_command(OUTPUT ${SPHINX_INDEX_FILE}
                   COMMAND
                   ${SPHINX_COMMAND} -c ${CMAKE_CURRENT_BINARY_DIR} -b html
                   ${SPHINX_SOURCE} ${SPHINX_BUILD}
                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                   DEPENDS ${ALL_RST} ${DOXYGEN_INDEX_FILE}
                   MAIN_DEPENDENCY ${CMAKE_CURRENT_BINARY_DIR}/conf.py
                   COMMENT "Generating documentation with Sphinx")

add_custom_target(Sphinx ALL DEPENDS ${SPHINX_INDEX_FILE})

# Install docs
include(GNUInstallDirs)
install(DIRECTORY ${SPHINX_BUILD} DESTINATION ${CMAKE_INSTALL_DOCDIR})
