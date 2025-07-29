# Changelog

All notable changes to ezARPACK will be documented in this file.

## [1.1] - Unreleased

* Compilation of MPI tests and examples is made opt-out via a CMake option
  `EnableMPI` (enabled by default). Thanks to @vincentmr for the contribution.
* Fixed a build failure caused by
  [opencollab/arpack-ng#377](https://github.com/opencollab/arpack-ng/pull/377).
* Applied compatibility fixes necessary for xtensor>=0.26.
* Check that Armadillo has been built with LAPACK support, skip the
  corresponding unit tests and examples otherwise.
* Do not build TRIQS unit tests and examples when a newer TRIQS version (>=3.1)
  is found (contribution from @thenoursehorse,
  [#7](https://github.com/krivenko/ezARPACK/pull/7)).
* Removed unused CMake option `Skip_arpack-ng`.
* Fixed flaky asymmetric MPI tests
  (issue [#8](https://github.com/krivenko/ezARPACK/issues/8)).
* Fixed deprecation warnings about default constructor of `nda::range()`.
* Set CMake policy `CMP0167` to `OLD` to silence warnings about
  `FindBoost.cmake`.
* Switched from the deprecated `FindPythonInterp.cmake` to `FindPython3.cmake`.
* Minor documentation fixes and improvements. In particular, compatibility with
  Sphinx>=7.2 has been improved.
* Both MathJax 2 and MathJax 3 are now supported for documentation compilation.
* Added a directory with developer tools (`tools`).
  
## [1.0] - 2022-09-04

* Wrappers for Parallel ARPACK with MPI message passing layer have been added.
  New wrapper classes `arpack_worker` are defined in a nested namespace
  `ezarpack::mpi` and have an API very similar to that of their serial
  counterparts.
* New accessor `arpack_solver::dim()` that returns dimension of the eigenproblem.
* Fixed a serious bug in the Blaze storage backend. By default, Blaze adds
  padding elements to data arrays when storing matrices. This fact was
  overlooked, which resulted in ARPACK procedures being called with wrong `LDV`
  and `LDZ` arguments. Resolving the issue required adding a new function,
  `storage_traits<Backend>::get_col_spacing()`.
* Drop support for the CMake option `EZARPACK_LINK_TO_ARPACK_NG`:
  `find_package(ezARPACK)` will never try to detect ARPACK-NG.
* Export a new CMake macro `find_arpackng()`. It finds a working installation of
  ARPACK-NG while dealing with version-to-version differences of ARPACK-NG's
  CMake interface.

## [0.10] - 2022-02-08

A new storage backend for the [TRIQS/nda](https://github.com/TRIQS/nda) library
has been added. This is the first release that is considered stable.

## [0.9] - 2020-07-04

Version 0.9 comes with a complete set of documentation and a few breaking
changes.

* A more common term 'solver' is now used instead of 'worker' everywhere in the
  code. In particular, `arpack_worker` has been renamed to `arpack_solver`.
* Method `arpack_solver::from_vector_n()` has been renamed to
  `arpack_solver::in_vector_n()`.
* Method `arpack_solver::to_vector_n()` has been renamed to
  `arpack_solver::out_vector_n()`.
* Computational mode `Invert` has been renamed to `Inverse`.
* Field `n_converged` has been removed from the `stats_t` structures.
  The number of converged Ritz values is now returned by method
  `arpack_solver::nconv()`.
