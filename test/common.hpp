/*******************************************************************************
 *
 * This file is part of ezARPACK, an easy-to-use C++ wrapper for
 * the ARPACK-NG FORTRAN library.
 *
 * Copyright (C) 2016-2022 Igor Krivenko <igor.s.krivenko@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 ******************************************************************************/
#pragma once

#include <functional>
#include <numeric>
#include <type_traits>
#include <utility>
#include <vector>

#include <catch2/catch.hpp>

#include "ezarpack/solver_base.hpp"

using namespace ezarpack;

// Some commonly used functions
template<operator_kind MKind>
using scalar_t =
    typename std::conditional<MKind == Complex, dcomplex, double>::type;

template<operator_kind MKind> scalar_t<MKind> reflect_coeff(scalar_t<MKind> x);
template<> double reflect_coeff<Symmetric>(double x) { return x; }
template<> double reflect_coeff<Asymmetric>(double x) { return -x; }
template<> dcomplex reflect_coeff<Complex>(dcomplex x) { return -x; }

// Set initial residual vector
template<operator_kind OpKind, typename Backend>
void set_init_residual_vector(arpack_solver<OpKind, Backend>& ar) {
  int const N = ar.dim();
  for(int i = 0; i < N; ++i)
    ar.residual_vector()[i] = double(i) / N;
}

//
// Exact Shift Strategy functors
//

template<operator_kind OpKind, typename VectorConstView, typename VectorView>
class exact_shift_strategy;

// Symmetric version
template<typename VectorConstView, typename VectorView>
class exact_shift_strategy<Symmetric, VectorConstView, VectorView> {

  std::vector<int> p;
  std::function<int(VectorView)> size_f;

public:
  template<typename SizeF>
  exact_shift_strategy(SizeF size_f) : size_f(size_f) {}

  // Symmetric version
  void operator()(VectorConstView ritz_values,
                  VectorConstView ritz_bounds,
                  VectorView shifts) {
    int np = size_f(shifts);
    if(np == 0) return;

    p.resize(np);
    std::iota(p.begin(), p.end(), 0);
    // p will contain the permutation that puts ritz_bounds in
    // the descending order of magnitude
    std::sort(p.begin(), p.end(),
              [&](int i, int j) { return ritz_bounds[i] > ritz_bounds[j]; });
    // Apply permutation p to ritz_values and use the result to fill shifts
    for(int i = 0; i < np; ++i)
      shifts[i] = ritz_values[p[i]];
  }
};

// Asymmetric version
template<typename VectorConstView, typename VectorView>
class exact_shift_strategy<Asymmetric, VectorConstView, VectorView> {

  std::vector<int> p;
  std::function<int(VectorView)> size_f;

public:
  template<typename SizeF>
  exact_shift_strategy(SizeF size_f) : size_f(std::move(size_f)) {}

  // Asymmetric version
  void operator()(VectorConstView ritz_values_re,
                  VectorConstView ritz_values_im,
                  VectorConstView ritz_bounds,
                  VectorView shifts_re,
                  VectorView shifts_im) {
    int np = size_f(shifts_re);
    if(np == 0) return;

    p.resize(np);
    std::iota(p.begin(), p.end(), 0);
    // p will contain the permutation that puts ritz_bounds in
    // the descending order of magnitude
    std::sort(p.begin(), p.end(),
              [&](int i, int j) { return ritz_bounds[i] > ritz_bounds[j]; });
    // Apply permutation p to ritz_values_* and use the result to fill shifts
    for(int i = 0; i < np; ++i) {
      shifts_re[i] = ritz_values_re[p[i]];
      shifts_im[i] = ritz_values_im[p[i]];
    }
  }
};

// Complex version
template<typename VectorConstView, typename VectorView>
class exact_shift_strategy<Complex, VectorConstView, VectorView> {

  std::vector<int> p;
  std::function<int(VectorView)> size_f;

public:
  template<typename SizeF>
  exact_shift_strategy(SizeF size_f) : size_f(std::move(size_f)) {}

  // Complex version
  void operator()(VectorConstView ritz_values,
                  VectorConstView ritz_bounds,
                  VectorView shifts) {
    int np = size_f(shifts);
    if(np == 0) return;

    p.resize(np);
    std::iota(p.begin(), p.end(), 0);
    // p will contain the permutation that puts ritz_bounds in
    // the descending order of magnitude
    std::sort(p.begin(), p.end(), [&](int i, int j) {
      return std::abs(ritz_bounds[i]) > std::abs(ritz_bounds[j]);
    });
    // Apply permutation p to ritz_values and use the result to fill shifts
    for(int i = 0; i < np; ++i)
      shifts[i] = ritz_values[p[i]];
  }
};
