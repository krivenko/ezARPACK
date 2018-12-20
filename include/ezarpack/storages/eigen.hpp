/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2018 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/
#pragma once

#include <complex>
#include <Eigen/Core>

#include "base.hpp"

namespace ezarpack {

/// Eigen3 storage backend tag
struct eigen_storage {};

/// Traits of the Eigen3 storage backend
template<> struct storage_traits<eigen_storage> {

  template<typename T> using vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;
  template<typename T> using matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

  template<typename T> using vector_view = Eigen::VectorBlock<vector<T>, Eigen::Dynamic>;
  template<typename T> using vector_const_view = Eigen::VectorBlock<const vector<T>, Eigen::Dynamic>;
  template<typename T> using matrix_view = Eigen::Block<matrix<T>, Eigen::Dynamic, Eigen::Dynamic, true>;
  template<typename T> using matrix_const_view = Eigen::Block<const matrix<T>, Eigen::Dynamic, Eigen::Dynamic, true>;

  // Storage types
  using real_vector_type = vector<double>;
  using complex_vector_type = vector<std::complex<double>>;
  using int_vector_type = vector<int>;

  using real_matrix_type = matrix<double>;
  using complex_matrix_type = matrix<std::complex<double>>;

  // View types
  using real_vector_view_type = vector_view<double>;
  using real_vector_const_view_type = vector_const_view<double>;
  using complex_vector_view_type = vector_view<std::complex<double>>;
  using complex_vector_const_view_type = vector_const_view<std::complex<double>>;

  using real_matrix_view_type = matrix_view<double>;
  using real_matrix_const_view_type = matrix_const_view<double>;
  using complex_matrix_view_type = matrix_view<std::complex<double>>;
  using complex_matrix_const_view_type = matrix_const_view<std::complex<double>>;

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
  inline static void resize(matrix<T> &m, int rows, int cols) { m.resize(rows, cols); }

  // Get pointer to data array
  template<typename T>
  inline static T* get_data_ptr(vector<T> &v) { return v.data(); }
  template<typename T>
  inline static T* get_data_ptr(matrix<T> &m) { return m.data(); }

  // Make vector view
  template<typename T>
  inline static vector_view<T> make_vector_view(vector<T> & v) {
    return v.head(v.size());
  }
  template<typename T>
  inline static vector_const_view<T> make_vector_const_view(vector<T> const& v) {
    return v.head(v.size());
  }

  // Make subvector view
  template<typename T>
  inline static vector_view<T> make_vector_view(vector<T> & v, int start, int size) {
    return v.segment(start, size);
  }
  template<typename T>
  inline static vector_const_view<T> make_vector_const_view(vector<T> const& v, int start, int size) {
    return v.segment(start, size);
  }

  // Make matrix view
  template<typename T>
  inline static matrix_view<T> make_matrix_view(matrix<T> & m) {
    return m.leftCols(m.cols());
  }
  template<typename T>
  inline static matrix_const_view<T> make_matrix_const_view(matrix<T> const& m) {
    return m.leftCols(m.cols());
  }

  // Make submatrix view including 'cols' leftmost columns
  template<typename T>
  inline static matrix_view<T> make_matrix_view(matrix<T> & m, int /* rows */, int cols) {
    return m.leftCols(cols);
  }
  template<typename T>
  inline static matrix_const_view<T> make_matrix_const_view(matrix<T> const& m, int /* rows */, int cols) {
    return m.leftCols(cols);
  }

  // Combine two real vectors to form one complex
  inline static complex_vector_type combine_re_im(real_vector_type const& re,
                                                  real_vector_type const& im,
                                                  int size) {
    return re.head(size) + std::complex<double>(0, 1) * im.head(size);
  }

  // Extract Ritz/Schur vectors from the 'z' matrix
  inline static complex_matrix_type extract_eigenvectors(real_matrix_type const& z,
                                                         real_vector_type const& di,
                                                         int N,
                                                         int nev) {
    complex_matrix_type res(N, nev);
    std::complex<double> I(0, 1);
    for(int i = 0; i < nev; ++i) {
      if(di(i) == 0) {
        res.col(i) = z.col(i);
      } else {
        res.col(i) = z.col(i) + I*std::copysign(1.0, di(i))*z.col(i+1);
        if(i < nev-1) {
          res.col(i+1) = res.col(i).conjugate();
          ++i;
        }
      }
    }
    return res;
  }
};

} // namespace ezarpack
