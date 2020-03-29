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

#include <armadillo>

#include "base.hpp"

namespace ezarpack {

/// Armadillo storage backend tag
struct armadillo_storage {};

/// Traits of the Armadillo storage backend
template<> struct storage_traits<armadillo_storage> {

  template<typename T> using vector = arma::Col<T>;
  template<typename T> using matrix = arma::Mat<T>;

  template<typename T> using vector_view = arma::subview_col<T>;
  template<typename T> using vector_const_view = const arma::subview_col<T>;
  template<typename T> using matrix_view = arma::subview<T>;
  template<typename T> using matrix_const_view = const arma::subview<T>;

  using dcomplex = std::complex<double>;
  using span = arma::span;

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
    return v.memptr();
  }
  template<typename T> inline static T* get_data_ptr(matrix<T>& m) {
    return m.memptr();
  }

  // Make vector view
  template<typename T>
  inline static vector_view<T> make_vector_view(vector<T>& v) {
    return v(span::all);
  }
  template<typename T>
  inline static vector_const_view<T>
  make_vector_const_view(vector<T> const& v) {
    return v(span::all);
  }

  // Make subvector view
  template<typename T>
  inline static vector_view<T>
  make_vector_view(vector<T>& v, int start, int size) {
    return v.subvec(start, start + size - 1);
  }
  template<typename T>
  inline static vector_const_view<T>
  make_vector_const_view(vector<T> const& v, int start, int size) {
    return v.subvec(start, start + size - 1);
  }

  // Make matrix view
  template<typename T>
  inline static matrix_view<T> make_matrix_view(matrix<T>& m) {
    return m(span::all, span::all);
  }
  template<typename T>
  inline static matrix_const_view<T>
  make_matrix_const_view(matrix<T> const& m) {
    return m(span::all, span::all);
  }

  // Make submatrix view including 'cols' leftmost columns
  template<typename T>
  inline static matrix_view<T>
  make_matrix_view(matrix<T>& m, int /* rows */, int cols) {
    return m.head_cols(cols);
  }
  template<typename T>
  inline static matrix_const_view<T>
  make_matrix_const_view(matrix<T> const& m, int /* rows */, int cols) {
    return m.head_cols(cols);
  }

  // worker_asymmetric: Extract Ritz values from 'dr' and 'di' vectors
  inline static complex_vector_type
  make_asymm_eigenvalues(real_vector_type const& dr,
                         real_vector_type const& di,
                         int nconv) {
    return dr.head(nconv) + dcomplex(0, 1) * di.head(nconv);
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
      if(di[i] == 0) {
        auto from1 = make_vector_const_view(z, i * N, N);
        a(from1, to1);
        lambda(i) = dot(from1, to1);
      } else {
        auto from1 = make_vector_const_view(z, i * N, N);
        auto from2 = make_vector_const_view(z, (i + 1) * N, N);
        a(from1, to1);
        a(from2, to2);
        lambda[i] = dcomplex(dot(from1, to1) + dot(from2, to2),
                             dot(from1, to2) - dot(from2, to1));
        ++i;
        lambda[i] = std::conj(lambda[i - 1]);
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
    dcomplex one(1);
    dcomplex I(0, 1);
    for(int i = 0; i < nconv; ++i) {
      if(di[i] == 0) {
        res.col(i) = one * z.subvec(i * N, (i + 1) * N - 1);
      } else {
        res.col(i) = z.subvec(i * N, (i + 1) * N - 1) +
                     I * std::copysign(1.0, di[i]) *
                         z.subvec((i + 1) * N, (i + 2) * N - 1);
        if(i < nconv - 1) {
          res.col(i + 1) = arma::conj(res.col(i));
          ++i;
        }
      }
    }
    return res;
  }
};

} // namespace ezarpack
