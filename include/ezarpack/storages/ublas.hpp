/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2019 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/
#pragma once

#include <cmath>
#include <complex>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include "base.hpp"

namespace ezarpack {

/// uBLAS storage backend tag
struct ublas_storage {};

/// Traits of the uBLAS storage backend
template<> struct storage_traits<ublas_storage> {

  template<typename T> using vector = boost::numeric::ublas::vector<T>;
  template<typename T>
  using matrix =
      boost::numeric::ublas::matrix<T, boost::numeric::ublas::column_major>;

  template<typename T>
  using vector_view = boost::numeric::ublas::vector_range<vector<T>>;
  template<typename T>
  using vector_const_view =
      const boost::numeric::ublas::vector_range<const vector<T>>;
  template<typename T>
  using matrix_view = boost::numeric::ublas::matrix_range<matrix<T>>;
  template<typename T>
  using matrix_const_view =
      const boost::numeric::ublas::matrix_range<const matrix<T>>;

  using dcomplex = std::complex<double>;

  // Storage types
  using real_vector_type = vector<double>;
  using complex_vector_type = vector<dcomplex>;
  using int_vector_type = vector<int>;

  using real_matrix_type = matrix<double>;
  using complex_matrix_type = matrix<dcomplex>;

  // View types
  using real_vector_view_type = vector_view<double>;
  using real_vector_const_view_type = vector_const_view<double>;
  using complex_vector_view_type = vector_view<dcomplex>;
  using complex_vector_const_view_type = vector_const_view<dcomplex>;

  using real_matrix_view_type = matrix_view<double>;
  using real_matrix_const_view_type = matrix_const_view<double>;
  using complex_matrix_view_type = matrix_view<dcomplex>;
  using complex_matrix_const_view_type = matrix_const_view<dcomplex>;

  // Factories
  inline static real_vector_type make_real_vector(int size) {
    return real_vector_type(size);
  }
  inline static complex_vector_type make_complex_vector(int size) {
    return complex_vector_type(size);
  }
  inline static int_vector_type make_int_vector(int size) {
    return int_vector_type(size);
  }
  inline static real_matrix_type make_real_matrix(int rows, int cols) {
    return real_matrix_type(rows, cols);
  }
  inline static complex_matrix_type make_complex_matrix(int rows, int cols) {
    return complex_matrix_type(rows, cols);
  }

  // Destructors (No-op)
  template<typename T> inline static void destroy(vector<T>&) {}
  template<typename T> inline static void destroy(matrix<T>&) {}

  // Resize
  template<typename T> inline static void resize(vector<T>& v, int size) {
    v.resize(size, false);
  }
  template<typename T>
  inline static void resize(matrix<T>& m, int rows, int cols) {
    m.resize(rows, cols, false);
  }

  // Get pointer to data array
  template<typename T> inline static T* get_data_ptr(vector<T>& v) {
    return &v(0);
  }
  template<typename T> inline static T* get_data_ptr(matrix<T>& m) {
    return &m(0, 0);
  }

  // Make vector view
  template<typename T>
  inline static vector_view<T> make_vector_view(vector<T>& v) {
    return boost::numeric::ublas::subrange(v, 0, v.size());
  }
  template<typename T>
  inline static vector_const_view<T>
  make_vector_const_view(vector<T> const& v) {
    return boost::numeric::ublas::subrange(v, 0, v.size());
  }

  // Make subvector view
  template<typename T>
  inline static vector_view<T>
  make_vector_view(vector<T>& v, int start, int size) {
    return boost::numeric::ublas::subrange(v, start, start + size);
  }
  template<typename T>
  inline static vector_const_view<T>
  make_vector_const_view(vector<T> const& v, int start, int size) {
    return boost::numeric::ublas::subrange(v, start, start + size);
  }

  // Make matrix view
  template<typename T>
  inline static matrix_view<T> make_matrix_view(matrix<T>& m) {
    using namespace boost::numeric::ublas;
    return project(m, range(0, m.size1()), range(0, m.size2()));
  }
  template<typename T>
  inline static matrix_const_view<T>
  make_matrix_const_view(matrix<T> const& m) {
    using namespace boost::numeric::ublas;
    return project(m, range(0, m.size1()), range(0, m.size2()));
  }

  // Make submatrix view including 'cols' leftmost columns
  template<typename T>
  inline static matrix_view<T>
  make_matrix_view(matrix<T>& m, int /* rows */, int cols) {
    using namespace boost::numeric::ublas;
    return project(m, range(0, m.size1()), range(0, cols));
  }
  template<typename T>
  inline static matrix_const_view<T>
  make_matrix_const_view(matrix<T> const& m, int /* rows */, int cols) {
    using namespace boost::numeric::ublas;
    return project(m, range(0, m.size1()), range(0, cols));
  }

  // worker_asymmetric: Extract Ritz values from 'dr' and 'di' vectors
  inline static complex_vector_type
  make_asymm_eigenvalues(real_vector_type const& dr,
                         real_vector_type const& di,
                         int nev) {
    complex_vector_type res(subrange(di, 0, nev));
    res *= dcomplex(0, 1);
    res += subrange(dr, 0, nev);
    return res;
  }

  // worker_asymmetric: Extract Ritz/Schur vectors from 'z' matrix
  inline static complex_matrix_type
  make_asymm_eigenvectors(real_matrix_type const& z,
                          real_vector_type const& di,
                          int N,
                          int nev) {
    complex_matrix_type res(N, nev);
    dcomplex I(0, 1);
    for(int i = 0; i < nev; ++i) {
      if(di(i) == 0) {
        column(res, i) = column(z, i);
      } else {
        column(res, i) = std::copysign(1.0, di(i)) * column(z, i + 1);
        column(res, i) *= I;
        column(res, i) += column(z, i);
        if(i < nev - 1) {
          column(res, i + 1) = conj(column(res, i));
          ++i;
        }
      }
    }
    return res;
  }
};

} // namespace ezarpack
