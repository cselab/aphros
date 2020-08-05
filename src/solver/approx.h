// Created by Petr Karnakov on 15.05.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <array>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "cond.h"
#include "embed.h"
#include "geom/mesh.h"
#include "solver.h"

template <class Idx, class Expr>
class Approx {
 public:
  virtual ~Approx() {}
  virtual Expr GetExpr(Idx) const = 0;
};

// sc: scheme:
//   - fou: first order upwind
//   - cd: central differences (mean value)
//   - sou: second order upwind
//   - quick: QUICK
template <class Scal>
std::array<Scal, 3> GetCoeff(ConvSc);

// Explicit interpolation to faces near inner cells.
// fc: field cell [s]
// fc: gradient [s]
// ffw: flow direction [i]
// sc: scheme:
// th: threshold for flow direction, ffw > th or ffw < -th
// Output:
// ff: face cell [i], resize if needed
template <class T, class M>
void InterpolateI(
    const FieldCell<T>& fc, const FieldCell<typename M::Vect>& fcgp,
    const FieldFace<T>& ffw, const M& m, ConvSc sc, typename M::Scal th,
    FieldFace<T>& ff);

// Implicit interpolation to inner faces with deferred correction.
// fc: field cell [s]
// fc: gradient [s]
// ffw: flow direction [i]
// sc: scheme:
//   - fou: first order upwind
//   - cd: central differences (mean value)
//   - sou: second order upwind
//   - quick: QUICK
// df: deferred correction factor (1. fully deferred)
// th: threshold for flow direction, ffw > th or ffw < -th
// Output:
// ff: face cell [i], resize if needed
template <class T, class M, class Expr>
void InterpolateI(
    const FieldCell<T>& fc, const FieldCell<typename M::Vect>& fcgp,
    const FieldFace<T>& ffw, FieldFace<Expr>& ff, const M& m, ConvSc sc,
    typename M::Scal df, typename M::Scal th);

// Explicit gradient on faces near inner cells.
// fc: field [s]
// Output:
// ff: normal gradient [i]
template <class M, class T>
void GradientI(const FieldCell<T>& fc, const M& m, FieldFace<T>& ff);

// Implicit gradient in inner faces.
// Output:
// ff: normal gradient [i]
template <class M, class Expr>
void GradientI(FieldFace<Expr>& ff, const M& m);

// Interpolates from nodes to faces
template <class T, class M>
FieldFace<T> Interpolate(const FieldNode<T>& fn, const M& m);

// Interpolation to inner faces.
// fc: field cell [s]
// Output:
// ff: face cell [i]
template <class T, class M>
void InterpolateI(const FieldCell<T>& fc, FieldFace<T>& ff, const M& m);

// Interpolation to support faces.
// fc: field cell [a]
// Output:
// ff: face cell [s]
template <class T, class M>
void InterpolateS(const FieldCell<T>& fc, FieldFace<T>& ff, const M& m);

// Linear extrapolation.
// xt: target
// x0,x1: points
// v0,v1: values
template <class T, class Scal>
T UExtrap(Scal xt, Scal x0, const T& v0, Scal x1, const T& v1) {
  return v0 + (v1 - v0) * ((xt - x0) / (x1 - x0));
}

// Returns average of fieldface.
// ff: fieldface [a]
// Output:
// fieldcell [a]
template <class T, class M>
FieldCell<T> Average(const FieldFace<T>& ff, const M& m);

// Smoothens fieldcell with node-based averaging.
// fc: fieldcell [s]
// rep: number of iterations
// Output:
// fc: smooth field [s]
template <class T, class M>
void SmoothenNode(FieldCell<T>& fc, M& m, size_t rep);

// Smoothens node field.
// fc: fieldcell [s]
// iters: number of iterations
// Output:
// fc: smooth field [s]
template <class T, class M>
void SmoothenNode(FieldNode<T>& fc, M& m, size_t iters);

// Returns gradient.
// ff: scalar fieldface [s]
// Output:
// gradient vector [s]
template <class M>
FieldCell<typename M::Vect> Gradient(
    const FieldFace<typename M::Scal>& ff, const M& m);

// Convention: Use Get/Set for fast procedures and Calc for those requiring
// computation: second(field, idx) vs GetNorm(field)

template <class Field, class M, class Scal = typename M::Scal>
Scal CalcDiff(const Field& fa, const Field& fb, const M& m);

template <class Idx, class M, class Scal = typename M::Scal>
Scal CalcDiff(
    const GField<typename M::Vect, Idx>& fa,
    const GField<typename M::Vect, Idx>& fb, const M& m);

// Coefficients for approximation of gradient with polynomial.
// x: target point
// z: stencil points
// Output:
// k: such that grad(x) = sum_i (ki * f(zi))
template <class Scal>
std::vector<Scal> GetGradCoeffs(Scal x, const std::vector<Scal>& z);

// Returns GetGradCoeffs(x,z[b:]) preceeded by b zeros.
template <class Scal>
std::vector<Scal> GetGradCoeffs(Scal x, const std::vector<Scal>& z, size_t b);

// Apply boundary conditions to halo cells
template <class T, class M>
void BcApply(FieldCell<T>& uc, const MapEmbed<BCond<T>>& me, const M& m);

// Apply reflection on all boundaries
// fill: value for other types that CondFaceReflect
template <class T, class M>
void BcReflectAll(FieldCell<T>& uc, const MapEmbed<BCond<T>>& me, const M& m);

template <class M, class Expr>
class FaceGrad : public Approx<IdxFace, Expr> {
 public:
  using Scal = typename M::Scal;
  explicit FaceGrad(const M& m) : m(m) {}
  Expr GetExpr(IdxFace f) const override {
    Expr e;
    IdxCell cm = m.GetCell(f, 0);
    IdxCell cp = m.GetCell(f, 1);
    // XXX: adhoc uniform
    Scal a = m.GetArea(f) / m.GetVolume(cp);
    e.InsertTerm(-a, cm);
    e.InsertTerm(a, cp);
    return e;
  }

 private:
  const M& m;
};

template <class Scal>
Scal Superbee(Scal p, Scal q) {
  if (p > 0. && q > 0.) {
    return std::max(std::min(2 * p, q), std::min(p, 2 * q));
  } else if (p < 0. && q < 0.) {
    return -std::max(std::min(-2 * p, -q), std::min(-p, -2 * q));
  }
  return 0.;
}

