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
#include <utility>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xcomplex.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xview.hpp>

#include "base.hpp"

namespace ezarpack {

/// xtensor storage backend tag
struct xtensor_storage {};

/// Traits of the xtensor storage backend
template<> struct storage_traits<xtensor_storage> {
private:
  // Implementation details

  template<typename T> using vector = xt::xtensor<T, 1>;
  template<typename T>
  using matrix = xt::xtensor<T, 2, xt::layout_type::column_major>;

  // xtensor documentation does not recommend to use xt::xview types explicitly.
  // Therefore, we have to manually detect the return type of xt::view() here.
  template<typename T>
  using vector_view =
      decltype(xt::view(std::declval<vector<T>&>(), xt::range(0, 0)));
  template<typename T>
  using vector_const_view =
      decltype(xt::view(std::declval<vector<T> const&>(), xt::range(0, 0)));
  template<typename T>
  using matrix_view = decltype(
      xt::view(std::declval<matrix<T>&>(), xt::range(0, 0), xt::range(0, 0)));
  template<typename T>
  using matrix_const_view = decltype(xt::view(std::declval<matrix<T> const&>(),
                                              xt::range(0, 0),
                                              xt::range(0, 0)));

  using dcomplex = std::complex<double>;

public:
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

  using real_matrix_const_view_type = matrix_const_view<double>;
  using complex_matrix_const_view_type = matrix_const_view<dcomplex>;

  // Factories
  inline static real_vector_type make_real_vector(int size) {
    return real_vector_type(real_vector_type::shape_type{size_t(size)});
  }
  inline static complex_vector_type make_complex_vector(int size) {
    return complex_vector_type(complex_vector_type::shape_type{size_t(size)});
  }
  inline static int_vector_type make_int_vector(int size) {
    return int_vector_type(int_vector_type::shape_type{size_t(size)});
  }
  inline static real_matrix_type make_real_matrix(int rows, int cols) {
    return real_matrix_type(
        real_matrix_type::shape_type{size_t(rows), size_t(cols)});
  }
  inline static complex_matrix_type make_complex_matrix(int rows, int cols) {
    return complex_matrix_type(
        complex_matrix_type::shape_type{size_t(rows), size_t(cols)});
  }

  // Destructors (No-op)
  template<typename T> inline static void destroy(vector<T>&) {}
  template<typename T> inline static void destroy(matrix<T>&) {}

  // Resize
  template<typename T> inline static void resize(vector<T>& v, int size) {
    v.resize(typename vector<T>::shape_type{size_t(size)});
  }
  template<typename T>
  inline static void resize(matrix<T>& m, int rows, int cols) {
    m.resize(typename matrix<T>::shape_type{size_t(rows), size_t(cols)});
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
    return xt::view(v, xt::range(0, v.size()));
  }

  // Make subvector view
  template<typename T>
  inline static vector_view<T>
  make_vector_view(vector<T>& v, int start, int size) {
    return xt::view(v, xt::range(start, start + size));
  }
  template<typename T>
  inline static vector_const_view<T>
  make_vector_const_view(vector<T> const& v, int start, int size) {
    return xt::view(v, xt::range(start, start + size));
  }

  // Make matrix view
  template<typename T>
  inline static matrix_const_view<T>
  make_matrix_const_view(matrix<T> const& m) {
    return xt::view(m, xt::range(0, m.shape(0)), xt::range(0, m.shape(1)));
  }

  // Make submatrix view including 'cols' leftmost columns
  template<typename T>
  inline static matrix_const_view<T>
  make_matrix_const_view(matrix<T> const& m, int /* rows */, int cols) {
    return xt::view(m, xt::range(0, m.shape(0)), xt::range(0, cols));
  }

  // worker_asymmetric: Extract Ritz values from 'dr' and 'di' vectors
  inline static complex_vector_type
  make_asymm_eigenvalues(real_vector_type const& dr,
                         real_vector_type const& di,
                         int nconv) {
    auto r = xt::range(0, nconv);
    return xt::view(dr, r) + dcomplex(0, 1) * xt::view(di, r);
  }

  // worker_asymmetric: Compute Ritz values from Ritz vectors
  template<typename A>
  inline static complex_vector_type
  make_asymm_eigenvalues(real_vector_type const& z,
                         real_vector_type const& di,
                         A&& a,
                         int N,
                         int nconv) {
    auto lambda = complex_vector_type::from_shape({size_t(nconv)});
    auto Ax1 = real_vector_type::from_shape({size_t(N)});
    auto Ax2 = real_vector_type::from_shape({size_t(N)});
    auto to1 = make_vector_view(Ax1);
    auto to2 = make_vector_view(Ax2);
    using xt::linalg::dot;
    for(int i = 0; i < nconv; ++i) {
      if(di(i) == 0) {
        auto from1 = make_vector_const_view(z, i * N, N);
        a(from1, to1);
        lambda(i) = dot(from1, to1)[0];
      } else {
        auto from1 = make_vector_const_view(z, i * N, N);
        auto from2 = make_vector_const_view(z, (i + 1) * N, N);
        a(from1, to1);
        a(from2, to2);
        lambda(i) = dcomplex(dot(from1, to1)[0] + dot(from2, to2)[0],
                             dot(from1, to2)[0] - dot(from2, to1)[0]);
        ++i;
        lambda(i) = std::conj(lambda(i - 1));
      }
    }
    return lambda;
  }

  // worker_asymmetric: Extract Ritz vectors from 'z' matrix
  inline static complex_matrix_type
  make_asymm_eigenvectors(real_vector_type const& z,
                          real_vector_type const& di,
                          int N,
                          int nconv) {
    auto res = complex_matrix_type::from_shape({size_t(N), size_t(nconv)});
    dcomplex I(0, 1);
    for(int i = 0; i < nconv; ++i) {
      if(di(i) == 0) {
        xt::view(res, xt::all(), i) =
            xt::view(z, xt::range(i * N, (i + 1) * N));
      } else {
        xt::view(res, xt::all(), i) =
            xt::view(z, xt::range(i * N, (i + 1) * N)) +
            I * std::copysign(1.0, di(i)) *
                xt::view(z, xt::range((i + 1) * N, (i + 2) * N));
        if(i < nconv - 1) {
          xt::view(res, xt::all(), i + 1) =
              xt::conj(xt::view(res, xt::all(), i));
          ++i;
        }
      }
    }
    return res;
  }
};

} // namespace ezarpack
