/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2024 Igor Krivenko
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

/// Blaze storage backend tag.
///
/// Passing this tag as the second template parameter of
/// @ref ezarpack::arpack_solver instructs it to use Blaze's vector/matrix/view
/// types.
struct blaze_storage {};

/// Traits of the Blaze storage backend.
///
/// Member typedefs of this structure describe types of objects that will be
/// used by arpack_solver to store numerical arrays and expose their partial
/// views. Structure's static member functions are called to control arrays'
/// lifetime, create the partial views, and perform some data post-processing
/// operations.
template<> struct storage_traits<blaze_storage> {
private:
  /// @name Private implementation details
  /// @{

  /// Blaze dynamic column vector.
  /// @tparam T Vector element type.
  template<typename T> using vector = blaze::DynamicVector<T>;
  /// Blaze dynamic column-major matrix.
  /// @tparam T Matrix element type.
  template<typename T>
  using matrix = blaze::DynamicMatrix<T, blaze::columnMajor>;

  /// Vector view type (subvector).
  /// @tparam T Vector element type.
  template<typename T> using vector_view = blaze::Subvector<vector<T>>;
  /// Constant vector view type (subvector of a constant vector).
  /// @tparam T Vector element type.
  template<typename T>
  using vector_const_view = blaze::Subvector<const vector<T>>;
  /// Matrix view type (submatrix).
  /// @tparam T Matrix element type.
  template<typename T> using matrix_view = blaze::Submatrix<matrix<T>>;
  /// Constant matrix view type (submatrix of a constant matrix).
  /// @tparam T Matrix element type.
  template<typename T>
  using matrix_const_view = blaze::Submatrix<const matrix<T>>;

  /// @}

  using dcomplex = std::complex<double>;

public:
  /// @name Vector and matrix storage types
  /// @{

  /// @brief One-dimensional container owning a contiguous array of `double`.
  using real_vector_type = vector<double>;
  /// @brief One-dimensional container owning a contiguous array of
  /// `std::complex<double>`.
  using complex_vector_type = vector<dcomplex>;
  /// @brief One-dimensional container owning a contiguous array of `int`.
  using int_vector_type = vector<int>;

  /// @brief Two-dimensional container owning a contiguous array of `double`.
  /// The storage order is column-major.
  using real_matrix_type = matrix<double>;
  /// @brief Two-dimensional container owning a contiguous array of
  /// `std::complex<double>`. The storage order is column-major.
  using complex_matrix_type = matrix<dcomplex>;

  /// @}

  /// @name View types
  /// @{

  /// Contiguous partial view of a real vector (subvector).
  using real_vector_view_type = vector_view<double>;
  /// Contiguous partial constant view of a real vector (subvector).
  using real_vector_const_view_type = vector_const_view<double>;
  /// Contiguous partial view of a complex vector (subvector).
  using complex_vector_view_type = vector_view<dcomplex>;
  /// Contiguous partial constant view of a complex vector (subvector).
  using complex_vector_const_view_type = vector_const_view<dcomplex>;

  /// Contiguous partial constant view of a real matrix (matrix block) that
  /// includes a number of the leftmost columns.
  using real_matrix_const_view_type = matrix_const_view<double>;
  /// Contiguous partial constant view of a complex matrix (matrix block) that
  /// includes a number of the leftmost columns.
  using complex_matrix_const_view_type = matrix_const_view<dcomplex>;

  /// @}

  /// @name Functions to create/destroy/resize data containers
  /// @{

  /// Constructs a real vector container.
  /// @param size Size of the vector.
  /// @return Constructed vector.
  inline static real_vector_type make_real_vector(int size) {
    return real_vector_type(size);
  }
  /// Constructs a complex vector container.
  /// @param size Size of the vector.
  /// @return Constructed vector.
  inline static complex_vector_type make_complex_vector(int size) {
    return complex_vector_type(size);
  }
  /// Constructs an integer vector container.
  /// @param size Size of the vector.
  /// @return Constructed vector.
  inline static int_vector_type make_int_vector(int size) {
    return int_vector_type(size);
  }
  /// Constructs a real matrix container.
  /// @param rows Number of matrix rows.
  /// @param cols Number of matrix columns.
  /// @return Constructed matrix.
  inline static real_matrix_type make_real_matrix(int rows, int cols) {
    return real_matrix_type(rows, cols);
  }
  /// Constructs a complex matrix container.
  /// @param rows Number of matrix rows.
  /// @param cols Number of matrix columns.
  /// @return Constructed matrix.
  inline static complex_matrix_type make_complex_matrix(int rows, int cols) {
    return complex_matrix_type(rows, cols);
  }

  /// Destroys a vector container.
  /// This function is, in fact, no-op, because memory owned by a Blaze vector
  /// will automatically be freed by its destructor.
  /// @tparam T Vector element type.
  template<typename T> inline static void destroy(vector<T>&) {}
  /// Destroys a matrix container.
  /// This function is, in fact, no-op, because memory owned by a Blaze matrix
  /// will automatically be freed by its destructor.
  /// @tparam T Matrix element type.
  template<typename T> inline static void destroy(matrix<T>&) {}

  /// Resizes a vector container.
  /// @tparam T Vector element type.
  /// @param v Vector container to resize.
  /// @param size New vector size.
  template<typename T> inline static void resize(vector<T>& v, int size) {
    v.resize(size);
  }
  /// Resizes a matrix container.
  /// @tparam T Matrix element type.
  /// @param m Matrix container to resize.
  /// @param rows New number of matrix rows.
  /// @param cols New number of matrix columns.
  template<typename T>
  inline static void resize(matrix<T>& m, int rows, int cols) {
    m.resize(rows, cols);
  }

  /// @}

  /// @name Access to underlying memory buffers
  /// @{

  /// Returns a pointer to the underlying data array owned by a vector.
  /// @tparam T Vector element type.
  /// @param v Vector to retrieve the data pointer from.
  /// @return Pointer to the data array.
  template<typename T> inline static T* get_data_ptr(vector<T>& v) {
    return v.data();
  }
  /// Returns a pointer to the underlying data array owned by a matrix.
  /// @tparam T Matrix element type.
  /// @param m Matrix to retrieve the data pointer from.
  /// @return Pointer to the data array.
  template<typename T> inline static T* get_data_ptr(matrix<T>& m) {
    return m.data();
  }

  /// Returns the spacing between the beginning of two columns of a matrix.
  ///
  /// The spacing can be different from the number of matrix rows if padding
  /// elements are added to the stored data array. Returning a negative value
  /// signals that the spacing equals the number of rows.
  /// @tparam T Matrix element type.
  /// @param m Matrix to retrieve the spacing from.
  /// @return Column spacing.
  template<typename T> inline static int get_col_spacing(matrix<T> const& m) {
    return m.spacing();
  }

  /// @}

  /// @name Functions to create vector/matrix views
  /// @{

  /// Makes a complete view of a vector.
  /// @tparam T Vector element type.
  /// @param v Vector container to make a view of.
  /// @return View of the full vector.
  template<typename T>
  inline static vector_view<T> make_vector_view(vector<T>& v) {
    return blaze::subvector(v, 0, v.size(), blaze::unchecked);
  }

  /// Makes a partial view of a vector.
  /// @tparam T Vector element type.
  /// @param v Vector container to make a view of.
  /// @param start Position of the starting element in the partial view.
  /// @param size Number of elements in the partial view.
  /// @return Subvector view.
  template<typename T>
  inline static vector_view<T>
  make_vector_view(vector<T>& v, int start, int size) {
    return blaze::subvector(v, start, size, blaze::unchecked);
  }
  /// Makes a constant partial view of a vector.
  /// @tparam T Vector element type.
  /// @param v Vector container to make a view of.
  /// @param start Position of the starting element in the partial view.
  /// @param size Number of elements in the partial view.
  /// @return Subvector view.
  template<typename T>
  inline static vector_const_view<T>
  make_vector_const_view(vector<T> const& v, int start, int size) {
    return blaze::subvector(v, start, size, blaze::unchecked);
  }

  /// Makes a complete view of a matrix.
  /// @tparam T Matrix element type.
  /// @param m Matrix container to make a view of.
  /// @return View of the full matrix.
  template<typename T>
  inline static matrix_const_view<T>
  make_matrix_const_view(matrix<T> const& m) {
    return blaze::submatrix(m, 0, 0, m.rows(), m.columns(), blaze::unchecked);
  }

  /// @brief Makes a constant partial view of a matrix including a number of
  /// the leftmost columns.
  /// @tparam T Matrix element type.
  /// @param m Matrix container to make a view of.
  /// @param rows Number of matrix rows.
  /// @param cols Number of the leftmost columns in the resulting view.
  /// @return Submatrix view.
  template<typename T>
  inline static matrix_const_view<T>
  make_matrix_const_view(matrix<T> const& m, int rows, int cols) {
    return blaze::submatrix(m, 0, 0, rows, cols, blaze::unchecked);
  }

  /// @}

  /// @name Post-processing required to compute eigenvalues/eigenvectors
  /// @{

  /// @brief Combines real and imaginary parts of eigenvalues computed by
  /// ezarpack::arpack_solver<Asymmetric, Backend>.
  ///
  /// @param dr Real parts of the computed eigenvalues.
  /// @param di Imaginary parts of the computed eigenvalues.
  /// @param nconv Number of the converged eigenvalues.
  /// @return Complex vector of eigenvalues.
  inline static complex_vector_type
  make_asymm_eigenvalues(real_vector_type const& dr,
                         real_vector_type const& di,
                         int nconv) {
    return blaze::subvector(dr, 0, nconv) +
           dcomplex(0, 1) * blaze::subvector(di, 0, nconv);
  }

  /// @brief Given a linear operator @f$ \hat A @f$, computes `nconv`
  /// eigenvalues from Ritz vectors @f$ \mathbf{x} @f$ as Rayleigh quotients
  /// @f$ \frac{\mathbf{x}^\dagger \hat A \mathbf{x}}
  ///          {\mathbf{x}^\dagger \hat M \mathbf{x}} @f$.
  /// This function is called by
  /// ezarpack::arpack_solver<Asymmetric, Backend>::eigenvalues(A &&) const.
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
    complex_vector_type lambda(nconv);
    real_vector_type Ax1(N), Ax2(N);
    auto out1 = make_vector_view(Ax1);
    auto out2 = make_vector_view(Ax2);
    for(int i = 0; i < nconv; ++i) {
      if(di[i] == 0) {
        auto in1 = make_vector_const_view(z, i * N, N);
        a(in1, out1);
        lambda[i] = dot(in1, out1);
      } else {
        auto in1 = make_vector_const_view(z, i * N, N);
        auto in2 = make_vector_const_view(z, (i + 1) * N, N);
        a(in1, out1);
        a(in2, out2);
        lambda[i] = dcomplex(dot(in1, out1) + dot(in2, out2),
                             dot(in1, out2) - dot(in2, out1));
        ++i;
        lambda[i] = std::conj(lambda[i - 1]);
      }
    }
    return lambda;
  }

  /// @brief Extracts `nconv` complex Ritz vectors from ARPACK-NG's internal
  /// representation. This function is called by
  /// ezarpack::arpack_solver<Asymmetric, Backend>::eigenvectors() const.
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
    complex_matrix_type res(N, nconv);
    dcomplex I(0, 1);
    using namespace blaze;
    for(int i = 0; i < nconv; ++i) {
      if(di[i] == 0) {
        column(res, i) = subvector(z, i * N, N);
      } else {
        column(res, i) =
            subvector(z, i * N, N) +
            I * std::copysign(1.0, di[i]) * subvector(z, (i + 1) * N, N);
        if(i < nconv - 1) {
          column(res, i + 1) = conj(column(res, i));
          ++i;
        }
      }
    }
    return res;
  }

  /// @}
};

} // namespace ezarpack
