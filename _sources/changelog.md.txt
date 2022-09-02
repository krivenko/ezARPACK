# Changelog

All notable changes to ezARPACK will be documented in this file.

## [1.0] - Unreleased

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
