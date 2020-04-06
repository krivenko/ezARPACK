/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2020 Igor Krivenko <igor.s.krivenko@gmail.com>
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
  template<typename T>
  using matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

  template<typename T>
  using vector_view = Eigen::VectorBlock<vector<T>, Eigen::Dynamic>;
  template<typename T>
  using vector_const_view = Eigen::VectorBlock<const vector<T>, Eigen::Dynamic>;
  template<typename T>
  using matrix_view =
      Eigen::Block<matrix<T>, Eigen::Dynamic, Eigen::Dynamic, true>;
  template<typename T>
  using matrix_const_view =
      Eigen::Block<const matrix<T>, Eigen::Dynamic, Eigen::Dynamic, true>;

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
    v.resize(size);
  }
  template<typename T>
  inline static void resize(matrix<T>& m, int rows, int cols) {
    m.resize(rows, cols);
  }

  // Get pointer to data array
  template<typename T> inline static T* get_data_ptr(vector<T>& v) {
    return v.data();
  }
  template<typename T> inline static T* get_data_ptr(matrix<T>& m) {
    return m.data();
  }

  // Make vector view
  template<typename T>
  inline static vector_view<T> make_vector_view(vector<T>& v) {
    return v.head(v.size());
  }
  template<typename T>
  inline static vector_const_view<T>
  make_vector_const_view(vector<T> const& v) {
    return v.head(v.size());
  }

  // Make subvector view
  template<typename T>
  inline static vector_view<T>
  make_vector_view(vector<T>& v, int start, int size) {
    return v.segment(start, size);
  }
  template<typename T>
  inline static vector_const_view<T>
  make_vector_const_view(vector<T> const& v, int start, int size) {
    return v.segment(start, size);
  }

  // Make matrix view
  template<typename T>
  inline static matrix_view<T> make_matrix_view(matrix<T>& m) {
    return m.leftCols(m.cols());
  }
  template<typename T>
  inline static matrix_const_view<T>
  make_matrix_const_view(matrix<T> const& m) {
    return m.leftCols(m.cols());
  }

  // Make submatrix view including 'cols' leftmost columns
  template<typename T>
  inline static matrix_view<T>
  make_matrix_view(matrix<T>& m, int /* rows */, int cols) {
    return m.leftCols(cols);
  }
  template<typename T>
  inline static matrix_const_view<T>
  make_matrix_const_view(matrix<T> const& m, int /* rows */, int cols) {
    return m.leftCols(cols);
  }

  // worker_asymmetric: Extract Ritz values from 'dr' and 'di' vectors
  inline static complex_vector_type
  make_asymm_eigenvalues(real_vector_type const& dr,
                         real_vector_type const& di,
                         int nconv) {
#ifdef EIGEN_CAN_MIX_REAL_COMPLEX_EXPR
    return dr.head(nconv) + dcomplex(0, 1) * di.head(nconv);
#else
    return dr.head(nconv).cast<dcomplex>() + dcomplex(0, 1) * di.head(nconv);
#endif
  }

  // worker_asymmetric: Compute Ritz values from Ritz vectors
  template<typename A>
  inline static complex_vector_type
  make_asymm_eigenvalues(real_vector_type const& z,
                         real_vector_type const& di,
                         A&& a,
                         int N,
                         int nconv) {
    complex_vector_type lambda(nconv);
    real_vector_type Ax1(N), Ax2(N);
    auto to1 = make_vector_view(Ax1);
    auto to2 = make_vector_view(Ax2);
    for(int i = 0; i < nconv; ++i) {
      if(di(i) == 0) {
        auto from1 = make_vector_const_view(z, i * N, N);
        a(from1, to1);
        lambda(i) = from1.dot(to1);
      } else {
        auto from1 = make_vector_const_view(z, i * N, N);
        auto from2 = make_vector_const_view(z, (i + 1) * N, N);
        a(from1, to1);
        a(from2, to2);
        lambda(i) = dcomplex(from1.dot(to1) + from2.dot(to2),
                             from1.dot(to2) - from2.dot(to1));
        ++i;
        lambda(i) = std::conj(lambda(i - 1));
      }
    }
    return lambda;
  }

  // worker_asymmetric: Extract Ritz/Schur vectors from 'z' matrix
  inline static complex_matrix_type
  make_asymm_eigenvectors(real_vector_type const& z,
                          real_vector_type const& di,
                          int N,
                          int nconv) {
    complex_matrix_type res(N, nconv);
    dcomplex I(0, 1);
    for(int i = 0; i < nconv; ++i) {
      if(di(i) == 0) {
#ifdef EIGEN_CAN_MIX_REAL_COMPLEX_EXPR
        res.col(i) = z.segment(i * N, N);
#else
        res.col(i) = z.segment(i * N, N).cast<dcomplex>();
#endif
      } else {
#ifdef EIGEN_CAN_MIX_REAL_COMPLEX_EXPR
        res.col(i) = z.segment(i * N, N) +
                     I * std::copysign(1.0, di(i)) * z.segment((i + 1) * N, N);
#else
        res.col(i) = z.segment(i * N, N).cast<dcomplex>() +
                     I * std::copysign(1.0, di(i)) * z.segment((i + 1) * N, N);
#endif
        if(i < nconv - 1) {
          res.col(i + 1) = res.col(i).conjugate();
          ++i;
        }
      }
    }
    return res;
  }
};

} // namespace ezarpack
