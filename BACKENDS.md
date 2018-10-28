Backend-dependent parts of `worker_*`
=====================================

Data members
------------

  * `vector<double>`
  * `vector<dcomplex>`
  * `matrix<double>`
  * `matrix<dcomplex>`
  * `vector<int>`

Members of `params_t`
---------------------

  * `vector<double>`
  * `vector<dcomplex>`

Arguments of `A::operator()`
----------------------------

  * `vector_view<double>`
  * `vector_const_view<double>`
  * `vector_view<dcomplex>`
  * `vector_const_view<dcomplex>`

Arguments of `ShiftsF::operator()`
----------------------------------

  * `vector_view<double>`
  * `vector_view<dcomplex>`

Return types of accessors
-------------------------

  * `vector_view<double>`
  * `vector_view<dcomplex>`
  * `matrix_view<double>`
  * `matrix_view<dcomplex>`

**Matrices have to be stored in the column-major (FORTRAN) order.**

Backends
========

Raw arrays
----------

  * `vector<T>` -> `T*`
  * `matrix<T>` -> `T*`
  * `vector_view<T>` -> `T*`
  * `vector_const_view<T>` -> `T const *`

TRIQS arrays
------------

  * `vector<T>` -> `triqs::arrays::vector<T>`
  * `matrix<T>` -> `triqs::arrays::matrix<T>`
  * `vector_view<T>` -> `triqs::arrays::vector_view<T>`
  * `vector_const_view<T>` -> `triqs::arrays::vector_const_view<T>`

Boost.MultiArray
----------------

  * `vector<T>` -> `boost::multi_array<T, 1>`
  * `matrix<T>` -> `boost::multi_array<T, 2>`
  * `vector_view<T>` -> `boost::multi_array_ref<T, 1>`
  * `vector_const_view<T>` -> `boost::const_multi_array_ref<T, 1>`

All containers are to be constructed with `store = fortran_storage_order()`.

Eigen 3
-------

  * `vector<double>` -> `Eigen::VectorXd`
  * `vector<dcomplex>` -> `Eigen::VectorXcd`
  * `vector<int>` -> `Eigen::VectorXi`
  * `matrix<double>` -> `Eigen::MatrixXd`
  * `matrix<dcomplex>` -> `Eigen::MatrixXcd`
  * `vector_view<double>` -> `Eigen::Block<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic, true>`
  * `vector_view<dcomplex>` -> `Eigen::Block<Eigen::VectorXcd, Eigen::Dynamic, Eigen::Dynamic, true>`
  * `vector_const_view<double>` -> `Eigen::Block<const Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic, true>`
  * `vector_const_view<dcomplex>` -> `Eigen::Block<const Eigen::VectorXcd, Eigen::Dynamic, Eigen::Dynamic, true>`

Blaze
-----

  * `vector<T>` -> `blaze::DynamicVector<T>`
  * `matrix<T>` -> `blaze::DynamicMatrix<T, blaze::columnMajor>`
  * `vector_view<T>` -> `decltype(blaze::submatrix(...))`
  * `vector_const_view<T>` -> `const decltype(blaze::submatrix(...))`

Armadillo
---------

  * `vector<double>` -> `arma::vec`
  * `vector<dcomplex>` -> `arma::cx_vec`
  * `vector<int>` -> `arma::Col<int>`
  * `matrix<double>` -> `arma::mat`
  * `matrix<dcomplex>` -> `arma::cx_mat`
  * `vector_view<double>` -> `arma::vec` constructed as `arma::vec(ptr_aux_mem, number_of_elements, false, true)`
  * `vector_view<dcomplex>` -> `arma::cx_vec` constructed as `arma::cx_vec(ptr_aux_mem, number_of_elements, false, true)`
  * `vector_const_view<double>` -> `const arma::vec` constructed as `arma::vec(ptr_aux_mem, number_of_elements, false, true)`
  * `vector_const_view<dcomplex>` -> `const arma::cx_vec` constructed as `arma::cx_vec(ptr_aux_mem, number_of_elements, false, true)`
