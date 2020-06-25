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

/// Raw memory storage backend tag.
///
/// Passing this tag as the second template parameter of arpack_worker
/// instructs it to use `std::unique_ptr` and raw pointers as vector/matrix and
/// view types respectively.
struct raw_storage {};

/// Traits of the raw memory storage backend.
///
/// Member typedefs of this structure describe types of objects that will be
/// used by arpack_worker to store numerical arrays and expose their partial
/// views. Structure's static member functions are called to control arrays'
/// lifetime, create the partial views, and perform some data post-processing
/// operations.
///
/// @warning This backend is not meant for general use. Its sole purpose is to
/// make unit testing possible even when no external vector/matrix algebra
/// library could be found.
template<> struct storage_traits<raw_storage> {
private:

  // Implementation details

  using dcomplex = std::complex<double>;

public:

  /// @name Vector and matrix storage types
  /// @{

  /// @brief One-dimensional container owning a contiguous array of `double`.
  using real_vector_type = std::unique_ptr<double[]>;
  /// @brief One-dimensional container owning a contiguous array of
  /// `std::complex<double>`.
  using complex_vector_type = std::unique_ptr<dcomplex[]>;
  /// @brief One-dimensional container owning a contiguous array of `int`.
  using int_vector_type = std::unique_ptr<int[]>;

  /// @brief Two-dimensional container owning a contiguous array of `double`.
  /// The storage order is column-major.
  using real_matrix_type = std::unique_ptr<double[]>;
  /// @brief Two-dimensional container owning a contiguous array of
  /// `std::complex<double>`. The storage order is column-major.
  using complex_matrix_type = std::unique_ptr<dcomplex[]>;

  /// @}

  /// @name View types
  /// @{

  /// Contiguous partial view of a real vector (subvector).
  using real_vector_view_type = double*;
  /// Contiguous partial constant view of a real vector (subvector).
  using real_vector_const_view_type = double const*;
  /// Contiguous partial view of a complex vector (subvector).
  using complex_vector_view_type = dcomplex*;
  /// Contiguous partial constant view of a complex vector (subvector).
  using complex_vector_const_view_type = dcomplex const*;

  /// Contiguous partial constant view of a real matrix (matrix block) that
  /// includes a number of the leftmost columns.
  using real_matrix_const_view_type = double const*;
  /// Contiguous partial constant view of a complex matrix (matrix block) that
  /// includes a number of the leftmost columns.
  using complex_matrix_const_view_type = dcomplex const*;

  /// @}

  /// @name Functions to create/destroy/resize data containers
  /// @{

  /// Constructs a real vector container.
  /// @param size Size of the vector.
  /// @return Constructed vector.
  inline static real_vector_type make_real_vector(int size) {
    return real_vector_type(new double[size]);
  }
  /// Constructs a complex vector container.
  /// @param size Size of the vector.
  /// @return Constructed vector.
  inline static complex_vector_type make_complex_vector(int size) {
    return complex_vector_type(new dcomplex[size]);
  }
  /// Constructs an integer vector container.
  /// @param size Size of the vector.
  /// @return Constructed vector.
  inline static int_vector_type make_int_vector(int size) {
    return int_vector_type(new int[size]);
  }
  /// Constructs a real matrix container.
  /// @param rows Number of matrix rows.
  /// @param cols Number of matrix columns.
  /// @return Constructed matrix.
  inline static real_matrix_type make_real_matrix(int rows, int cols) {
    return real_matrix_type(new double[rows * cols]);
  }
  /// Constructs a complex matrix container.
  /// @param rows Number of matrix rows.
  /// @param cols Number of matrix columns.
  /// @return Constructed matrix.
  inline static complex_matrix_type make_complex_matrix(int rows, int cols) {
    return complex_matrix_type(new dcomplex[rows * cols]);
  }

  /// Destroys a vector/matrix container.
  /// @tparam T Vector/matrix element type.
  template<typename T> inline static void destroy(std::unique_ptr<T[]>& p) {}

  /// Resizes a vector container.
  /// @tparam T Vector element type.
  /// @param v Vector container to resize.
  /// @param size New vector size.
  template<typename T>
  inline static void resize(std::unique_ptr<T[]>& v, int size) {
    v.reset(new T[size]);
  }
  /// Resizes a matrix container.
  /// @tparam T Matrix element type.
  /// @param m Matrix container to resize.
  /// @param rows New number of matrix rows.
  /// @param cols New number of matrix columns.
  template<typename T>
  inline static void resize(std::unique_ptr<T[]>& m, int rows, int cols) {
    m.reset(new T[rows * cols]);
  }

  /// Returns a pointer to the underlying data array owned by a vector/matrix.
  /// @tparam T Vector/matrix element type.
  /// @param vm Vector/matrix to retrieve the data pointer from.
  /// @return Pointer to the data array.
  template<typename T> inline static T* get_data_ptr(std::unique_ptr<T[]>& vm) {
    return vm.get();
  }

  /// @}

  /// @name Functions to create vector/matrix views
  /// @{

  /// Makes a complete view of a vector.
  /// @tparam T Vector element type.
  /// @param v Vector container to make a view of.
  /// @return View of the full vector.
  template<typename T>
  inline static T* make_vector_view(std::unique_ptr<T[]>& v) {
    return v.get();
  }

  /// Makes a partial view of a vector.
  /// @tparam T Vector element type.
  /// @param v Vector container to make a view of.
  /// @param start Position of the starting element in the partial view.
  /// @param size **[ignored]** Number of elements in the partial view.
  /// @return Subvector view.
  template<typename T>
  inline static T*
  make_vector_view(std::unique_ptr<T[]>& v, int start, int size) {
    return v.get() + start;
  }
  /// Makes a constant partial view of a vector.
  /// @tparam T Vector element type.
  /// @param v Vector container to make a view of.
  /// @param start Position of the starting element in the partial view.
  /// @param size **[ignored]** Number of elements in the partial view.
  /// @return Subvector view.
  template<typename T>
  inline static T const* make_vector_const_view(std::unique_ptr<T[]> const& v,
                                                int start,
                                                int size) {
    return v.get() + start;
  }

  /// Makes a complete view of a matrix.
  /// @tparam T Matrix element type.
  /// @param m Matrix container to make a view of.
  /// @return View of the full matrix.
  template<typename T>
  inline static T const* make_matrix_const_view(std::unique_ptr<T[]> const& m) {
    return m.get();
  }

  /// @brief Makes a constant partial view of a matrix including a number of
  /// the leftmost columns.
  /// @tparam T Matrix element type.
  /// @param m Matrix container to make a view of.
  /// @param rows **[ignored]** Number of matrix rows.
  /// @param cols **[ignored]** Number of the leftmost columns in the resulting
  /// view.
  /// @return Submatrix view.
  template<typename T>
  inline static T const* make_matrix_const_view(std::unique_ptr<T[]> const& m,
                                                int rows,
                                                int cols) {
    return m.get();
  }

  /// @}

  /// @name Post-processing required to compute eigenvalues/eigenvectors
  /// @{

  /// @brief Combines real and imaginary parts of eigenvalues computed by
  /// ezarpack::arpack_worker<Asymmetric, Backend>.
  ///
  /// @param dr Real parts of the computed eigenvalues.
  /// @param di Imaginary parts of the computed eigenvalues.
  /// @param nconv Number of the converged eigenvalues.
  /// @return Complex vector of eigenvalues.
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

  /// @brief Given a linear operator @f$ \hat A @f$, computes `nconv`
  /// eigenvalues from Ritz vectors @f$ \mathbf{x} @f$ as Rayleigh quotients
  /// @f$ \frac{\mathbf{x}^\dagger \hat A \mathbf{x}}
  ///          {\mathbf{x}^\dagger \hat M \mathbf{x}} @f$.
  /// This function is called by
  /// ezarpack::arpack_worker<Asymmetric, Backend>::eigenvalues(A &&) const.
  ///
  /// @tparam A Type of the callable object representing the linear operator
  /// @f$ \hat A @f$.
  /// @param z Holds components of the Ritz vectors @f$ \mathbf{x} @f$ as
  /// a sequence of `nconv` length-`N` chunks. Meaning of each chunk depends on
  /// the corresponding component of `di`, see below.
  /// @param di If `di[i]` is zero, then the `i`-th chunk of `z` contains a
  /// real Ritz vector. Otherwise, `di[i] = -di[i+1] != 0`, in which case
  /// the `i`-th and `(i+1)`-th chunks of `z` are
  /// real and imaginary parts of a complex Ritz vector respectively. Every
  /// real Ritz vector corresponds to one real eigenvalue, while every
  /// complex vector gives rise to a complex conjugate pair of eigenvalues.
  /// @param a Callable object representing @f$ \hat A @f$.
  /// @param N Dimension of the eigenproblem.
  /// @param nconv Number of the converged Ritz vectors.
  /// @return Complex vector of eigenvalues.
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

  /// @brief Extracts `nconv` complex Ritz vectors from ARPACK-NG's internal
  /// representation. This function is called by
  /// ezarpack::arpack_worker<Asymmetric, Backend>::eigenvectors() const.
  ///
  /// @param z Holds components of the Ritz vectors @f$ \mathbf{x} @f$ as
  /// a sequence of `nconv` length-`N` chunks. Meaning of each chunk depends on
  /// the corresponding component of `di`, see below.
  /// @param di If `di[i]` is zero, then the `i`-th chunk of `z` contains a
  /// real Ritz vector. Otherwise, `di[i] = -di[i+1] != 0`, in which case
  /// the `i`-th and `(i+1)`-th chunks of `z` are
  /// real and imaginary parts of a complex Ritz vector respectively. Every such
  /// pair corresponds to a complex conjugate pair of Ritz vectors, so that the
  /// total amount of vectors stored in `z` is exactly `nconv`.
  /// @param N Dimension of the eigenproblem.
  /// @param nconv Number of the converged Ritz vectors.
  /// @return Complex matrix, whose columns are Ritz vectors (eigenvectors).
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

  /// @}
};

} // namespace ezarpack
