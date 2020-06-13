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

#include <cmath>
#include <complex>
#include <memory>

#include "base.hpp"

namespace ezarpack {

/// Raw memory storage backend tag
struct raw_storage {};

/// Traits of the raw memory storage backend
template<> struct storage_traits<raw_storage> {
private:
  // Implementation details

  using dcomplex = std::complex<double>;

public:
  // Storage types
  using real_vector_type = std::unique_ptr<double[]>;
  using complex_vector_type = std::unique_ptr<dcomplex[]>;
  using int_vector_type = std::unique_ptr<int[]>;

  using real_matrix_type = std::unique_ptr<double[]>;
  using complex_matrix_type = std::unique_ptr<dcomplex[]>;

  // View types
  using real_vector_view_type = double*;
  using real_vector_const_view_type = double const*;
  using complex_vector_view_type = dcomplex*;
  using complex_vector_const_view_type = dcomplex const*;

  using real_matrix_const_view_type = double const*;
  using complex_matrix_const_view_type = dcomplex const*;

  // Factories
  inline static real_vector_type make_real_vector(int size) {
    return real_vector_type(new double[size]);
  }
  inline static complex_vector_type make_complex_vector(int size) {
    return complex_vector_type(new dcomplex[size]);
  }
  inline static int_vector_type make_int_vector(int size) {
    return int_vector_type(new int[size]);
  }
  inline static real_matrix_type make_real_matrix(int rows, int cols) {
    return real_matrix_type(new double[rows * cols]);
  }
  inline static complex_matrix_type make_complex_matrix(int rows, int cols) {
    return complex_matrix_type(new dcomplex[rows * cols]);
  }

  // Destructors (No-op)
  template<typename T> inline static void destroy(std::unique_ptr<T[]>& p) {}

  // Resize
  template<typename T>
  inline static void resize(std::unique_ptr<T[]>& v, int size) {
    v.reset(new T[size]);
  }
  template<typename T>
  inline static void resize(std::unique_ptr<T[]>& m, int rows, int cols) {
    m.reset(new T[rows * cols]);
  }

  // Get pointer to data array
  template<typename T> inline static T* get_data_ptr(std::unique_ptr<T[]>& p) {
    return p.get();
  }

  // Make vector view
  template<typename T>
  inline static T* make_vector_view(std::unique_ptr<T[]>& v) {
    return v.get();
  }

  // Make subvector view
  template<typename T>
  inline static T*
  make_vector_view(std::unique_ptr<T[]>& v, int start, int /* size */) {
    return v.get() + start;
  }
  template<typename T>
  inline static T const* make_vector_const_view(std::unique_ptr<T[]> const& v,
                                                int start,
                                                int /* size */) {
    return v.get() + start;
  }

  // Make matrix view
  template<typename T>
  inline static T const* make_matrix_const_view(std::unique_ptr<T[]> const& m) {
    return m.get();
  }

  // Make submatrix view including 'cols' leftmost columns
  template<typename T>
  inline static T const* make_matrix_const_view(std::unique_ptr<T[]> const& m,
                                                int /* rows */,
                                                int /* cols */) {
    return m.get();
  }

  // worker_asymmetric: Extract Ritz values from 'dr' and 'di' vectors
  inline static complex_vector_type
  make_asymm_eigenvalues(real_vector_type const& dr,
                         real_vector_type const& di,
                         int nconv) {
    complex_vector_type res(new dcomplex[nconv]);
    for(int n = 0; n < nconv; ++n) {
      res[n] = dcomplex(dr[n], di[n]);
    }
    return res;
  }

  // worker_asymmetric: Compute Ritz values from Ritz vectors
  template<typename A>
  inline static complex_vector_type
  make_asymm_eigenvalues(real_vector_type const& z,
                         real_vector_type const& di,
                         A&& a,
                         int N,
                         int nconv) {
    complex_vector_type lambda(new dcomplex[nconv]);
    real_vector_type Ax1(new double[N]), Ax2(new double[N]);
    auto to1 = make_vector_view(Ax1);
    auto to2 = make_vector_view(Ax2);
    auto dot = [&](real_vector_const_view_type v1,
                   real_vector_const_view_type v2) {
      double s = 0;
      for(int i = 0; i < N; ++i)
        s += v1[i] * v2[i];
      return s;
    };
    for(int i = 0; i < nconv; ++i) {
      if(di[i] == 0) {
        auto from1 = make_vector_const_view(z, i * N, N);
        a(from1, to1);
        lambda[i] = dot(from1, to1);
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

  // worker_asymmetric: Extract Ritz vectors from 'z' matrix
  inline static complex_matrix_type
  make_asymm_eigenvectors(real_vector_type const& z,
                          real_vector_type const& di,
                          int N,
                          int nconv) {
    complex_matrix_type res(new dcomplex[N * nconv]);
    dcomplex I(0, 1);
    for(int i = 0; i < nconv; ++i) {
      if(di[i] == 0) {
        for(int n = 0; n < N; ++n)
          res[n + N * i] = z[n + N * i];
      } else {
        for(int n = 0; n < N; ++n) {
          res[n + N * i] =
              z[n + N * i] + I * std::copysign(1.0, di[i]) * z[n + N * (i + 1)];
        }
        if(i < nconv - 1) {
          for(int n = 0; n < N; ++n) {
            res[n + N * (i + 1)] = std::conj(res[n + N * i]);
          }
          ++i;
        }
      }
    }
    return res;
  }
};

} // namespace ezarpack
