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

#include <complex>
#include <blaze/Math.h>

#include "base.hpp"

namespace ezarpack {

/// Blaze storage backend tag
struct blaze_storage {};

/// Traits of the Blaze storage backend
template<> struct storage_traits<blaze_storage> {

  template<typename T> using vector = blaze::DynamicVector<T>;
  template<typename T> using matrix =
    blaze::DynamicMatrix<T, blaze::columnMajor>;

  template<typename T> using vector_view = blaze::Subvector<vector<T>>;
  template<typename T> using vector_const_view =
    blaze::Subvector<const vector<T>>;
  template<typename T> using matrix_view = blaze::Submatrix<matrix<T>>;
  template<typename T> using matrix_const_view =
    blaze::Submatrix<const matrix<T>>;

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
  template<typename T> inline static void destroy(vector<T> &) {}
  template<typename T> inline static void destroy(matrix<T> &) {}

  // Resize
  template<typename T>
  inline static void resize(vector<T> &v, int size) { v.resize(size); }
  template<typename T>
  inline static void resize(matrix<T> &m, int rows, int cols) {
    m.resize(rows, cols);
  }

  // Get pointer to data array
  template<typename T>
  inline static T* get_data_ptr(vector<T> &v) { return v.data(); }
  template<typename T>
  inline static T* get_data_ptr(matrix<T> &m) { return m.data(); }

  // Make vector view
  template<typename T>
  inline static vector_view<T> make_vector_view(vector<T> & v) {
    return blaze::subvector(v, 0, v.size(), blaze::unchecked);
  }
  template<typename T>
  inline static
  vector_const_view<T> make_vector_const_view(vector<T> const& v) {
    return blaze::subvector(v, 0, v.size(), blaze::unchecked);
  }

  // Make subvector view
  template<typename T>
  inline static vector_view<T> make_vector_view(vector<T> & v,
                                                int start,
                                                int size) {
    return blaze::subvector(v, start, size, blaze::unchecked);
  }
  template<typename T>
  inline static vector_const_view<T> make_vector_const_view(vector<T> const& v,
                                                            int start,
                                                            int size) {
    return blaze::subvector(v, start, size, blaze::unchecked);
  }

  // Make matrix view
  template<typename T>
  inline static matrix_view<T> make_matrix_view(matrix<T> & m) {
    return blaze::submatrix(m, 0, 0, m.rows(), m.columns(), blaze::unchecked);
  }
  template<typename T>
  inline static
  matrix_const_view<T> make_matrix_const_view(matrix<T> const& m) {
    return blaze::submatrix(m, 0, 0, m.rows(), m.columns(), blaze::unchecked);
  }

  // Make submatrix view including 'cols' leftmost columns
  template<typename T>
  inline static
  matrix_view<T> make_matrix_view(matrix<T> & m, int rows, int cols) {
    return blaze::submatrix(m, 0, 0, rows, cols, blaze::unchecked);
  }
  template<typename T>
  inline static matrix_const_view<T> make_matrix_const_view(matrix<T> const& m,
                                                            int rows,
                                                            int cols) {
    return blaze::submatrix(m, 0, 0, rows, cols, blaze::unchecked);
  }

  // worker_asymmetric: Extract Ritz values from 'dr' and 'di' vectors
  inline static
  complex_vector_type make_asymm_eigenvalues(real_vector_type const& dr,
                                             real_vector_type const& di,
                                             int nev) {
    return blaze::subvector(dr, 0, nev) +
           dcomplex(0, 1) * blaze::subvector(di, 0, nev);
  }

  // worker_asymmetric: Extract Ritz/Schur vectors from 'z' matrix
  inline static
  complex_matrix_type make_asymm_eigenvectors(real_matrix_type const& z,
                                              real_vector_type const& di,
                                              int N,
                                              int nev) {
    complex_matrix_type res(N, nev);
    dcomplex I(0, 1);
    using namespace blaze;
    for(int i = 0; i < nev; ++i) {
      if(di[i] == 0) {
        column(res, i) = column(z, i);
      } else {
        column(res, i) = column(z, i) +
                         I*std::copysign(1.0, di[i])*column(z, i+1);
        if(i < nev-1) {
          column(res, i+1) = conj(column(res, i));
          ++i;
        }
      }
    }
    return res;
  }
};

} // namespace ezarpack
