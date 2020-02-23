// Created by Petr Karnakov on 11.02.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <array>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "cond.h"
#include "embed.h"
#include "geom/mesh.h"
#include "parse/vars.h"
#include "solver.h"

template <class Scal_>
struct ULinear {
  static constexpr size_t dim = 3;
  using Scal = Scal_;
  using Vect = generic::Vect<Scal, dim>;

  template <size_t N>
  static std::array<Scal, N> Mul(
      std::array<Scal, N * N> a, std::array<Scal, N> x) {
    using Int = size_t;
    std::array<Scal, N> r;
    for (Int i = 0; i < N; ++i) {
      r[i] = 0;
      for (Int j = 0; j < N; ++j) {
        r[i] += a[i * N + j] * x[j];
      }
    }
    return r;
  }

  // Solves linear system a*x=b.
  template <size_t N>
  static std::array<Scal, N> SolveLinear(
      std::array<Scal, N * N> a, std::array<Scal, N> b) {
    using Int = size_t;
    auto aa = [&a](Int i, Int j) -> Scal& { return a[i * N + j]; };
    auto swaprows = [&aa, &b](Int i, Int ip) {
      if (i == ip) {
        return;
      }
      for (Int j = 0; j < N; ++j) {
        std::swap(aa(i, j), aa(ip, j));
      }
      std::swap(b[i], b[ip]);
    };
    auto ipivot = [&aa](const Int j) {
      Int imax = j;
      for (Int i = j + 1; i < N; ++i) {
        if (std::abs(aa(i, j)) > std::abs(aa(imax, j))) {
          imax = i;
        }
      }
      return imax;
    };
    auto addrow = [&aa, &b](Int i, Int ip, Scal ap) {
      for (Int j = 0; j < N; ++j) {
        aa(i, j) += aa(ip, j) * ap;
      }
      b[i] += b[ip] * ap;
    };
    for (Int j = 0; j < N; ++j) {
      const Int ip = ipivot(j);
      swaprows(ip, j);
      for (Int i = j + 1; i < N; ++i) {
        addrow(i, j, -aa(i, j) / aa(j, j));
      }
    }
    std::array<Scal, N> x;
    for (Int i = N; i > 0;) {
      --i;
      Scal t = b[i];
      for (Int j = i + 1; j < N; ++j) {
        t -= aa(i, j) * x[j];
      }
      x[i] = t / aa(i, i);
    }
    return x;
  }

  // Fits linear function to set of points and values
  //   u = g.dot(x) + u0
  // Returns {g, u0}.
  static std::pair<Vect, Scal> FitLinear(
      const std::vector<Vect>& xx, const std::vector<Scal>& uu);

  // Fits linear function to set of points and values
  //   u = g.dot(x) + u0
  // Returns {g, u0}.
  static std::pair<generic::Vect<Vect, dim>, Vect> FitLinear(
      const std::vector<Vect>& xx, const std::vector<Vect>& uu);

  template <class T>
  static T EvalLinear(
      const std::pair<generic::Vect<T, dim>, T>& p, const Vect& x) {
    auto& g = p.first;
    auto& u0 = p.second;
    return g[0] * x[0] + g[1] * x[1] + g[2] * x[2] + u0;
  }
};

template <class M_>
struct UEmbed {
  using M = M_;
  using EB = Embed<M>;
  static constexpr size_t dim = M::dim;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Type = typename EB::Type;

  struct CondEmbed {
    enum class Type { value, gradient };
    Scal u; // value or normal gradient
  };

  static FieldNode<Scal> InitEmbed(const M& m, const Vars& var, bool verb);

  // feu: field on embedded boundaries [a]
  // Returns:
  // field on cells [a]
  template <class T>
  static FieldCell<T> Interpolate(const FieldEmbed<T>& feu, const EB& eb);

  // Mid-point interpolation.
  // fcu: field [a]
  // bc: boundary conditions type, 0: value, 1: grad
  // bcv: value or normal gradient (grad dot GetNormal)
  // Returns:
  // field on embedded boundaries [s]
  template <class T>
  static FieldEmbed<T> Interpolate(
      const FieldCell<T>& fcu, const MapCondFace& mfc, size_t bc, T bcv,
      const EB& eb);

  // Upwind interpolation.
  // fcu: field [a]
  // ffv: volume flux
  // mfc: conditions on faces
  // bc: boundary conditions type, 0: value, 1: grad
  // bcv: value or normal gradient (grad dot GetNormal)
  // Returns:
  // field on embedded boundaries [s]
  template <class T>
  static FieldEmbed<T> InterpolateUpwind(
      const FieldCell<T>& fcu, const FieldEmbed<Scal>& fev,
      const MapCondFace& mfc, size_t bc, T bcv, const EB& eb);

  // Interpolates field from cells to faces with an upwind scheme,
  // bilinear in cut faces, and first order upwind in embed faces.
  // fcu: field [s]
  // fcg: gradient [s]
  // ffv: volume flux [i]
  // sc: interpolation scheme
  // Returns:
  // field on faces boundaries [s]
  static FieldFace<Scal> InterpolateUpwindBilinear(
      const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcg,
      const FieldFace<Scal>& ffv, ConvSc sc, const EB& eb);

  static FieldEmbed<Scal> InterpolateUpwindBilinear(
      const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcg,
      const MapCondFace& mfc, size_t bc, const MapCell<Scal>& mcu,
      const FieldEmbed<Scal>& fev, ConvSc sc, const EB& eb);

  static FieldEmbed<Scal> InterpolateUpwindBilinear(
      const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcg,
      const MapCondFace& mfc, size_t bc, Scal bcv, const FieldEmbed<Scal>& fev,
      ConvSc sc, const EB& eb);

  // Upwind interpolation.
  // fcu: field [s]
  // fcg: gradient [s]
  // mfc: conditions on faces
  // bc: boundary conditions type, 0: value, 1: grad
  // bcv: value or normal gradient (grad dot GetNormal)
  // ffv: volume flux [i]
  // sc: interpolation scheme
  // Returns:
  // field on embedded boundaries [s]
  // TODO: interpolation from upwind cell on embedded boundaries
  //       (see InterpolateUpwind() above)
  static FieldEmbed<Scal> InterpolateUpwind(
      const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcg,
      const MapCondFace& mfc, size_t bc, const MapCell<Scal>& mcu,
      const FieldEmbed<Scal>& fev, ConvSc sc, const EB& eb);

  static FieldEmbed<Scal> InterpolateUpwind(
      const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcg,
      const MapCondFace& mfc, size_t bc, Scal bcv, const FieldEmbed<Scal>& fev,
      ConvSc sc, const EB& eb);

  // feu: field on embedded boundaries [a]
  // Returns:
  // gradient on cells [a]
  static FieldCell<Vect> GradientGauss(
      const FieldEmbed<Scal>& feu, const EB& eb);

  // Gradient from linear fit to face centers.
  // feu: field on embedded boundaries [a]
  // Returns:
  // gradient on cells [a]
  static FieldCell<Vect> GradientLinearFit(
      const FieldEmbed<Scal>& feu, const EB& eb);

  template <class T>
  static FieldCell<T> AverageCutCells(const FieldCell<T>& fcu, const EB& eb);

  template <class T>
  static FieldCell<T> RedistributeCutCells(
      const FieldCell<T>& fcu, const EB& eb);

  // Gradient on faces with limited denominator.
  // fcu: field [a]
  // bc: boundary conditions type, 0: value, 1: grad
  // bcv: value or normal gradient (grad dot GetNormal)
  // Returns:
  // grad dot GetNormal on embedded boundaries [s]
  template <class T>
  static FieldEmbed<T> GradientLimited(
      const FieldCell<T>& fcu, const MapCondFace& mfc, size_t bc, T bcv,
      const EB& eb);

  // Updates flux in cut faces using bilinear interpolation from regular faces
  // (Schwartz,2006).
  // feu: field [a]
  // bc: boundary conditions type, 0: value, 1: grad
  // bcv: value or normal gradient (grad dot GetNormal)
  // Returns:
  // grad dot GetNormal on embedded boundaries [s]
  template <class T>
  static FieldFace<T> InterpolateBilinearFaces(
      const FieldFace<T>& ffu, const EB& eb);

  // Interpolates field from cells to faces, bilinear in cut faces.
  // fcu: field on cells [s]
  // Returns:
  // ffu: field on faces [i]
  template <class T>
  static FieldFace<T> InterpolateBilinear(
      const FieldCell<T>& fcu, const EB& eb);

  // Interpolates field from cells to embed faces.
  // fcu: field on cells [s]
  // mcu: value of on embed faces
  // Returns:
  // feu: field on embed faces [i]
  template <class T>
  static void InterpolateEmbedFaces(
      const FieldCell<T>& fcu, size_t bc, const MapCell<T>& mcu,
      FieldEmbed<T>& feu, const EB& eb);

  // Interpolates field from cells to embed faces.
  // fcu: field on cells [s]
  // mcu: value of on embed faces
  // Returns:
  // feu: field on embed faces [i]
  template <class T>
  static void InterpolateUpwindEmbedFaces(
      const FieldCell<T>& fcu, size_t bc, const MapCell<Scal>& mcu,
      const FieldEmbed<Scal>& fev, FieldEmbed<T>& feu, const EB& eb);

  // Interpolates field from cells to combined field, bilinear in cut faces.
  // fcu: field on cells [s]
  // mcu: value of on embed faces
  // Returns:
  // feu: field on embed faces [i]
  template <class T>
  static FieldEmbed<T> InterpolateBilinear(
      const FieldCell<T>& fcu, size_t bc, const MapCell<T>& mcu, const EB& eb);

  template <class T>
  static FieldEmbed<T> InterpolateBilinear(
      const FieldCell<T>& fcu, size_t bc, T bcv, const EB& eb);

  template <class T>
  static FieldEmbed<T> InterpolateBilinear(
      const FieldCell<T>& fcu, const MapCondFace& mfc, size_t bc, T bcv,
      const EB& eb);

  template <class T>
  static FieldEmbed<T> Interpolate(
      const FieldCell<T>& fcu, const MapEmbed<BCond<T>>& mebc, const EB& eb);
  template <class T>
  static FieldEmbed<T> Gradient(
      const FieldCell<T>& fcu, const MapEmbed<BCond<T>>& mebc, const EB& eb);
  static FieldEmbed<Scal> InterpolateUpwind(
      const FieldCell<Scal>& fcu, const MapEmbed<BCond<Scal>>& mebc, ConvSc sc,
      const FieldCell<Vect>& fcg, const FieldEmbed<Scal>& fev, const EB& eb);

  template <class T>
  static FieldFace<T> Interpolate(
      const FieldCell<T>& fcu, const MapEmbed<BCond<T>>& mebc, const M& m);
  template <class T>
  static FieldFace<T> Gradient(
      const FieldCell<T>& fcu, const MapEmbed<BCond<T>>& mebc, const M& m);
  static FieldFace<Scal> InterpolateUpwind(
      const FieldCell<Scal>& fcu, const MapEmbed<BCond<Scal>>& mebc, ConvSc sc,
      const FieldCell<Vect>& fcg, const FieldFace<Scal>& ffv, const M& m);

  // Gradient with bilinear interpolation in cut faces
  // and linear fit in embed faces.
  // fcu: field [a]
  // Returns:
  // grad dot GetNormal on faces [s]
  template <class T>
  static FieldFace<T> GradientBilinear(const FieldCell<T>& fcu, const EB& eb);

  // Gradient with bilinear interpolation in cut faces
  // and linear fit in embed faces.
  // fcu: field [a]
  // mcu: value (bc=0) or derivative (bc=1)
  // Returns:
  // grad dot GetNormal on embedded boundaries [s]
  template <class T>
  static FieldEmbed<T> GradientBilinear(
      const FieldCell<T>& fcu, size_t bc, const MapCell<T>& mcu, const EB& eb);

  template <class T>
  static FieldEmbed<T> GradientBilinear(
      const FieldCell<T>& fcu, const MapCondFace& mfc, size_t bc, T bcv,
      const EB& eb);

  template <class Expr>
  static Expr GetExprVal(IdxFace f, const BCond<Scal>& bc, const M& m) {
    Expr e;
    const auto nci = bc.nci;
    const IdxCell cm = m.GetCell(f, 0);
    const IdxCell cp = m.GetCell(f, 1);
    const IdxCell c = m.GetCell(f, nci);
    e.InsertTerm(0, cm);
    e.InsertTerm(0, cp);
    switch (bc.type) {
      case BCondType::dirichlet: {
        e.SetConstant(bc.val);
        break;
      }
      case BCondType::neumann: {
        const Scal g = (nci == 0 ? 1. : -1.);
        const Scal a = m.GetVolume(c) / m.GetArea(f) * 0.5 * g;
        e.SetConstant(a * bc.val);
        e.InsertTerm(1., c);
        break;
      }
      default:
        throw std::runtime_error(std::string() + __func__ + ": unknown");
    }
    return e;
  }
  template <class Expr>
  static Expr GetExprGrad(IdxFace f, const BCond<Scal>& bc, const M& m) {
    Expr e;
    const auto nci = bc.nci;
    const IdxCell c = m.GetCell(f, nci);
    const IdxCell cm = m.GetCell(f, 0);
    const IdxCell cp = m.GetCell(f, 1);
    e.InsertTerm(0, cm);
    e.InsertTerm(0, cp);
    switch (bc.type) {
      case BCondType::dirichlet: {
        const Scal g = (nci == 0 ? 1. : -1.);
        const Scal hr = m.GetArea(f) / m.GetVolume(c);
        const Scal a = hr * 2 * g;
        e.SetConstant(a * bc.val);
        e.InsertTerm(-a, c);
        break;
      }
      case BCondType::neumann: {
        e.SetConstant(bc.val);
        break;
      }
      default:
        throw std::runtime_error(std::string() + __func__ + ": unknown");
    }
    return e;
  }
};
