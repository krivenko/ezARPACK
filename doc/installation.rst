.. _installation:

Installation instructions
=========================

ezARPACK is a header-only C++ library, which means it is usable without
installation. You can obtain the latest version of the library by making a local
clone of its `GitHub repository <https://github.com/krivenko/ezARPACK>`_.
Assuming that `git <https://git-scm.com/>`_ package is installed and visible via
the `PATH` environment variable, the following shell command will create a
directory `ezARPACK.src` containing the latest ezARPACK source files.

.. code-block:: shell

    $ git clone https://github.com/krivenko/ezARPACK.git ezARPACK.src

You can then use ezARPACK in your own code by passing
'`-I/path/to/ezARPACK.src/include`' to the compiler command line and adding
'`-L/path/to/ARPACK-NG/lib -larpack`' to the linker command line.

If your project uses `CMake <https://cmake.org/download/>`_ as its build system,
you might want to install ezARPACK header files along with CMake configuration
scripts. This will allow you to import ezARPACK targets, as well as to run unit
tests and build example programs. The minimum required version of CMake is
**3.1.0**.

The following sequence of shell commands will build unit tests and examples.

.. code-block:: shell

    $ mkdir ezARPACK.build && cd ezARPACK.build
    $ cmake ../ezARPACK.src                                   \
    $ -DCMAKE_INSTALL_PREFIX=<ezARPACK_installation_prefix>   \
      -DARPACK_NG_ROOT=<ARPACK-NG_installation_prefix>        \
      -DEigen3_ROOT=<Eigen3_installation_prefix>              \
      -Dblaze_ROOT=<Blaze_installation_prefix>                \
      -DArmadillo_ROOT=<Armadillo_installation_prefix>        \
      -DBOOST_ROOT=<Boost_installation_prefix>                \
      -DTRIQS_ROOT=<TRIQS_installation_prefix>                \
      -Dxtensor_ROOT=<xtensor_installation_prefix>            \
      -Dxtensor-blas_ROOT=<xtensor-blas_installation_prefix>  \
      -DExamples=ON                                           \
      -DTests=ON
    $ make

Compilation of the tests can be disabled with CMake flag `-DTests=OFF`
*(not recommended)*. Examples are compiled by default and can be disabled
with `-DExamples=OFF`.

It is generally recommended to run the unit tests before installing the library
to detect potential problems early on.

.. code-block:: shell

    $ make test

Failed tests (having anything different from `Passed` in the test result column)
signal host-specific compilation/linking problems or bugs in the library itself.
The following command completes installation of library's files into
`<ezARPACK_installation_prefix>`.

.. code-block:: shell

    $ make install

All unit tests and examples require
`ARPACK-NG <https://github.com/opencollab/arpack-ng>`_ to be built. Most of them
also depend on external matrix manipulation libraries and will be automatically
disabled if a usable installation of a respective library cannot be found.
The `*_ROOT` CMake options allow to manually specify installation locations of
said libraries.

Documentation of ezARPACK can optionally be built and installed using the
`Decumentation` CMake flag (requires `Doxygen <https://www.doxygen.nl/>`_,
`Sphinx <https://www.sphinx-doc.org>`_,
`Breathe <https://breathe.readthedocs.io>`_ and
`MathJax <https://www.mathjax.org/>`_).

The table below gives a complete lists of supported CMake options with their
meaning.

+-------------------------+----------------------------------------------------+
| Option name             | Description                                        |
+=========================+====================================================+
| CMAKE_INSTALL_PREFIX    | Path to the directory ezARPACK will be installed   |
|                         | into.                                              |
+-------------------------+----------------------------------------------------+
| CMAKE_BUILD_TYPE        | CMake build type (`Release`, `Debug` or            |
|                         | `RelWithDebInfo`) used to compile unit tests and   |
|                         | examples.                                          |
+-------------------------+----------------------------------------------------+
| Tests=[ON|OFF]          | Enable/disable compilation of unit tests.          |
+-------------------------+----------------------------------------------------+
| Examples=[ON|OFF]       | Enable/disable compilation of example programs.    |
+-------------------------+----------------------------------------------------+
| ARPACK_NG_ROOT          | Path to ARPACK-NG installation.                    |
+-------------------------+----------------------------------------------------+
| Eigen3_ROOT             | Path to Eigen 3 installation.                      |
+-------------------------+----------------------------------------------------+
| blaze_ROOT              | Path to Blaze installation.                        |
+-------------------------+----------------------------------------------------+
| Armadillo_ROOT          | Path to Armadillo installation.                    |
+-------------------------+----------------------------------------------------+
| BOOST_ROOT              | Path to Boost (uBLAS) installation.                |
+-------------------------+----------------------------------------------------+
| TRIQS_ROOT              | Path to TRIQS installation.                        |
+-------------------------+----------------------------------------------------+
| xtensor_ROOT            | Path to xtensor installation.                      |
+-------------------------+----------------------------------------------------+
| xtensor-blas_ROOT       | Path to xtensor-blas installation.                 |
+-------------------------+----------------------------------------------------+
| Documentation=[ON|OFF]  | Enable/disable generation of ezARPACK's            |
|                         | Doxygen and Sphinx documentation.                  |
+-------------------------+----------------------------------------------------+
| Doxygen_ROOT            | Path to Doxygen installation.                      |
+-------------------------+----------------------------------------------------+
| Sphinx_ROOT             | Path to Sphinx installation.                       |
+-------------------------+----------------------------------------------------+
| MathJax_ROOT            | Path to MathJax installation (directory containing |
|                         | `MathJax.js`).                                     |
+-------------------------+----------------------------------------------------+
