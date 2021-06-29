// Created by Petr Karnakov on 30.07.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <memory>

#include "geom/mesh.h"
#include "solver/cond.h"

template <class M_>
class UNormal {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Vect2 = generic::Vect<Scal, 2>;
  using Vect3 = generic::Vect<Scal, 3>;
  using Vect4 = generic::Vect<Scal, 4>;

  // Computes normal using combined Youngs scheme and height-functions
  // m: mesh
  // fcu: volume fraction
  // fci: interface mask (1: contains interface)
  // edim: effective dimension
  // Output: set to NaN if fci=0
  // fcn: normal with norm1()=1, antigradient of fcu [s]
  static void CalcNormal(
      M& m, const FieldCell<Scal>& fcu, const FieldCell<bool>& fci, size_t edim,
      FieldCell<Vect>& fcn);

  // Computes normal by Youngs scheme
  // m: mesh
  // fcu: volume fraction
  // fci: interface mask (1: contains interface)
  // Output: set to NaN if fci=0
  // fcn: normal with norm1()=1, antigradient of fcu [s]
  static void CalcNormalYoungs(
      M& m, const FieldCell<Scal>& fcu, const FieldCell<bool>& fci,
      FieldCell<Vect>& fcn);

  // u: volume fraction, array of size 3x3x3
  static Vect2 GetNormalYoungs(const std::array<Scal, 9>& u);
  static Vect3 GetNormalYoungs(const std::array<Scal, 27>& u);
  static Vect4 GetNormalYoungs(const std::array<Scal, 81>& u);
  // u: volume fraction, array of size 3x3x3
  // n: guess for normal, updated if heights give steeper estimate
  static void GetNormalHeight(const std::array<Scal, 9>& u, Vect2& n);
  static void GetNormalHeight(const std::array<Scal, 27>& u, Vect3& n);
  static void GetNormalHeight(const std::array<Scal, 81>& u, Vect4& n);

 public:
  struct Imp; // implementation
};
