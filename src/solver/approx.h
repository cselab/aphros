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

// Interpolates from cells to inner faces.
// T: value type (Scal or Vect)
// fc: field cell [s]
// fcgp: gradient [s]
// ffw: flow direction [i]
// sc: scheme:
// th: threshold for flow direction, ffw > th or ffw < -th
// Output:
// ff: face cell [i], resize if needed
template <class T, class M>
void Interpolate(
    const FieldCell<T>& fc, const FieldCell<typename M::Vect>& fcgp,
    const MapCondFace& mfc, const FieldFace<T>& ffw, const M& m, ConvSc sc,
    typename M::Scal th, FieldFace<T>& ff);

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

// Explicit gradient on boundary faces.
// fc: field [s]
// mfc: face conditions
// Output:
// ff: normal gradient [i]
template <class M, class T>
void GradientB(
    const FieldCell<T>& fc, const MapCondFace& mfc, const M& m,
    FieldFace<T>& ff);

// Explicit gradient on inner faces.
// fc: field [s]
// Output:
// ff: normal gradient [i]
template <class M, class T>
void Gradient(
    const FieldCell<T>& fc, const MapCondFace& mfc, const M& m,
    FieldFace<T>& ff);

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

// Interpolation to faces with defined conditions.
// fc: field cell [i]
// mfc: face cond
// Output:
// ff: values updated on faces defined in mfc
template <class T, class M>
void InterpolateB(
    const FieldCell<T>& fc, const MapCondFace& mfc, FieldFace<T>& ff,
    const M& m);

// Interpolates from cells to support faces.
// T: value type (Scal or Vect)
// fc: field cell [a]
// mfc: face cond
// Output:
// field face [s]
template <class T, class M>
FieldFace<T> Interpolate(
    const FieldCell<T>& fc, const MapCondFace& mfc, const M& m);

// Second order upwind interpolation with TVD Superbee limiter
// fc: fieldcell [a]
// fcg: gradient of field [a]
// mfc: face cond
// ffw: flow direction [s]
// Output:
// fieldface [s]
template <class M>
FieldFace<typename M::Scal> InterpolateSuperbee(
    const FieldCell<typename M::Scal>& fc,
    const FieldCell<typename M::Vect>& fcg, const MapCondFace& mfc,
    const FieldFace<typename M::Scal>& ffw, const M& m,
    typename M::Scal th = 1e-8);

// Returns average of fieldface.
// ff: fieldface [a]
// Output:
// fieldcell [a]
template <class T, class M>
FieldCell<T> Average(const FieldFace<T>& ff, const M& m);

// Smoothens fieldcell.
// fc: fieldcell [s]
// mfc: condface
// rep: number of iterations
// Output:
// fc: smooth field [s]
template <class T, class M>
void Smoothen(FieldCell<T>& fc, const MapCondFace& mfc, M& m, size_t rep);

// Smoothens fieldcell with node-based averaging.
// fc: fieldcell [s]
// rep: number of iterations
// Output:
// fc: smooth field [s]
template <class T, class M>
void SmoothenNode(FieldCell<T>& fc, M& m, size_t rep);

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
void BcApply(FieldCell<T>& uc, const MapCondFace& mfc, const M& m);
template <class T, class M>
void BcApply(FieldCell<T>& uc, const MapEmbed<BCond<T>>& me, const M& m);

// Apply reflection on all boundaries
// fill: value for other types that CondFaceReflect
template <class T, class M>
void BcReflectAll(FieldCell<T>& uc, const MapCondFace& mfc, const M& m);
template <class T, class M>
void BcReflectAll(FieldCell<T>& uc, const MapEmbed<BCond<T>>& me, const M& m);

template <class M, class Expr>
class FaceValB : public Approx<IdxFace, Expr> {
 public:
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  FaceValB(const M& m, const MapCondFace& mfc) : m(m), mfc_(mfc) {}
  Expr GetExpr(IdxFace f) const override {
    Expr e;
    if (auto cb = mfc_.find(f)) {
      IdxCell cm = m.GetCell(f, 0);
      IdxCell cp = m.GetCell(f, 1);
      e.InsertTerm(0, cm);
      e.InsertTerm(0, cp);
      if (auto cd = cb->template Get<CondFaceVal<Scal>>()) {
        e.SetConstant(cd->second());
      } else if (auto cd = cb->template Get<CondFaceGrad<Scal>>()) {
        size_t id = cd->GetNci();
        IdxCell c = m.GetCell(f, id);
        Scal g = (id == 0 ? 1. : -1.);
        Scal a = m.GetVolume(c) / m.GetArea(f) * 0.5 * g;
        e.SetConstant(a * cd->GetGrad());
        e.InsertTerm(1., c);
      } else {
        throw std::runtime_error("FaceValB: unknown cond");
      }
    } else {
      throw std::runtime_error("FaceValB: unset cond");
    }
    return e;
  }

 private:
  const M& m;
  const MapCondFace& mfc_;
};

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

template <class M, class Expr>
class FaceGradB : public Approx<IdxFace, Expr> {
 public:
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  explicit FaceGradB(const M& m, const MapCondFace& mfc) : m(m), mfc_(mfc) {}
  Expr GetExpr(IdxFace f) const override {
    Expr e;
    if (auto cb = mfc_.find(f)) {
      IdxCell cm = m.GetCell(f, 0);
      IdxCell cp = m.GetCell(f, 1);
      e.InsertTerm(0, cm);
      e.InsertTerm(0, cp);
      if (auto cd = cb->template Get<CondFaceGrad<Scal>>()) {
        e.SetConstant(cd->GetGrad());
      } else if (auto cd = cb->template Get<CondFaceVal<Scal>>()) {
        size_t id = cd->GetNci();
        IdxCell c = m.GetCell(f, id);
        Scal g = (id == 0 ? 1. : -1.);
        Scal hr = m.GetArea(f) / m.GetVolume(c);
        Scal a = hr * 2 * g;
        e.SetConstant(a * cd->second());
        e.InsertTerm(-a, c);
      } else {
        throw std::runtime_error("FaceGradB: unknown cond");
      }
    } else {
      throw std::runtime_error("FaceGradB: unset cond");
    }
    return e;
  }

 private:
  const M& m;
  const MapCondFace& mfc_;
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

