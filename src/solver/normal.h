#pragma once

#include <memory>

#include "geom/mesh.h"

namespace solver {

template <class M_>
class UNormal {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  // Computes normal by combined Youngs scheme and height-functions
  // m: mesh
  // fcu: volume fraction
  // fci: interface mask (1: contains interface)
  // edim: effective dimension
  // Output: set to NaN if fci=0
  // fcn: normal with norm1()=1, antigradient of fcu [s]
  static void CalcNormal(
      M& m, const FieldCell<Scal>& fcu, const FieldCell<bool>& fci,
      size_t edim, FieldCell<Vect>& fcn);

  // Computes normal by combined Youngs scheme and height-functions
  // m: mesh
  // fcu: volume fraction [a]
  // fcud2: volume fraction difference double (xpp-xmm, ...) [a]
  // fcud4: volume fraction difference quad (xp4-xm4, ...) [a]
  // edim: effective dimension
  // Output: set to NaN if fci=0
  // fch: curvature [i]
  static void CalcHeight(
      M& m, const FieldCell<Scal>& fcu,
      const FieldCell<Vect>& fcud2,
      const FieldCell<Vect>& fcud4, const FieldCell<Vect>& fcud6,
      size_t edim, FieldCell<Vect>& fch);

  // Computes normal by combined Youngs scheme and height-functions
  // m: mesh
  // fcu: volume fraction [a]
  // fch: height function [a]
  // fcn: normal to defining the direction of heights as n.abs().argmax()
  // edim: effective dimension
  // Output: set to NaN if fci=0
  // fck: curvature [i]
  static void CalcCurvHeight(
      M& m, const FieldCell<Scal>& fcu, const FieldCell<Vect>& fch,
      const FieldCell<Vect>& fcn, size_t edim, FieldCell<Scal>& fck);

  // Computes normal by Youngs scheme
  // m: mesh
  // fcu: volume fraction
  // fci: interface mask (1: contains interface)
  // Output: set to NaN if fci=0
  // fcn: normal with norm1()=1, antigradient of fcu [s]
  static void CalcNormalYoung(
      M& m, const FieldCell<Scal>& fcu, const FieldCell<bool>& fci,
      FieldCell<Vect>& fcn);

 public:
  struct Imp; // implementation
};

} // namespace solver
