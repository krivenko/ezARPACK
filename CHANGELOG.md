# Changelog

All notable changes to ezARPACK will be documented in this file.

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
