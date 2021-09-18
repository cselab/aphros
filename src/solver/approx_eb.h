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

namespace interp {

// Evaluates quadratic interpolant on points -1, 0, 1.
// x: target point
// um,u,up: values of function for x=-1, 0, 1.
template <class Scal, class T>
T Quad(Scal x, T um, T u, T up) {
  return (um * (x - 1) + up * (x + 1)) * x * 0.5 - u * (x - 1) * (x + 1);
}

// Evaluates linear interpolant on points 0, 1.
// x: target point
// u,up: values of function for x=0, 1.
template <class Scal, class T>
T Linear(Scal x, T u, T up) {
  return u * (1 - x) + up * x;
}

// Evaluates bilinear interpolant on points (0,0), (1,0), (0,1) and (1,1).
// x,y: target point
// u,ux,uy,uyx:  values of function for (x,y) = (0,0), (1,0), (0,1), (1,1)
template <class T, class Scal>
T Bilinear(Scal x, Scal y, T u, T ux, T uy, T uyx) {
  //                      //
  //   y                  //
  //   |                  //
  //   |*uy    *uyx       //
  //   |                  //
  //   |                  //
  //   |*u     *ux        //
  //   |-------------x    //
  //                      //
  const auto v = u * (1 - x) + ux * x;
  const auto vy = uy * (1 - x) + uyx * x;
  return v * (1 - y) + vy * y;
}

} // namespace interp

template <class Vect_>
struct ULinearFit {
  using Scal = typename Vect_::Scal;
  using Vect = Vect_;
  static constexpr size_t dim = Vect_::dim;

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
      if (aa(j, j) == 0) { // case of degenerate system
        aa(j, j) = 1;
      }
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
    T res = p.second;
    for (size_t d = 0; d < dim; ++d) {
      res += g[d] * x[d];
    }
    return res;
  }
};

template <class EB, class T>
auto FitLinear(IdxCell c, const FieldCell<T>& fcu, const EB& eb) {
  using Vect = typename EB::Vect;
  auto& m = eb.GetMesh();
  std::vector<Vect> xx;
  std::vector<T> uu;
  for (auto cn : eb.Stencil(c)) {
    xx.push_back(m.GetCenter(cn));
    uu.push_back(fcu[cn]);
  }
  return ULinearFit<Vect>::FitLinear(xx, uu);
}

template <class EB, class T>
T EvalLinearFit(
    typename EB::Vect x, IdxCell c, const FieldCell<T>& fcu, const EB& eb) {
  auto p = FitLinear(c, fcu, eb);
  return ULinearFit<typename EB::Vect>::EvalLinear(p, x);
}

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

  static void InitLevelSet(FieldNode<Scal>&, M& m, const Vars& var, bool verb);

  // feu: field on embedded boundaries [a]
  // Returns:
  // field on cells [a]
  template <class T>
  static FieldCell<T> Interpolate(const FieldEmbed<T>& feu, const EB& eb);
  template <class T>
  static FieldCell<T> Interpolate(const FieldFace<T>& feu, const M& m);

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
  template <class T>
  static FieldCell<T> RedistributeCutCells(const FieldCell<T>& fcu, const M& m);

  static FieldCell<Scal> RedistributeCutCellsAdvection(
      const FieldCell<Scal>& fcs, const FieldFace<Scal>& ffv, Scal cfl, Scal dt,
      const M& m);
  static FieldCell<Scal> RedistributeCutCellsAdvection(
      const FieldCell<Scal>& fcs, const FieldEmbed<Scal>& ffv, Scal cfl,
      Scal dt, const EB& eb);

  static void RedistributeConstTerms(
      FieldCell<typename M::Expr>& fce, const Embed<M>& eb, M& m);
  static void RedistributeConstTerms(
      FieldCell<typename M::Expr>& fce, const M& eb, M& m);

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
  template <class T>
  static FieldFace<T> InterpolateBilinearFaces(
      const FieldFace<T>& ffu, const M& m);

  template <class T>
  static FieldEmbed<T> Interpolate(
      const FieldCell<T>& fcu, const MapEmbed<BCond<T>>& mebc, const EB& eb);
  template <class T>
  static FieldEmbed<T> Gradient(
      const FieldCell<T>& fcu, const MapEmbed<BCond<T>>& mebc, const EB& eb);
  static FieldEmbed<Scal> InterpolateUpwind(
      const FieldCell<Scal>& fcu, const MapEmbed<BCond<Scal>>& mebc, ConvSc sc,
      const FieldCell<Vect>& fcg, const FieldEmbed<Scal>& fev, const EB& eb);
  // Bell-Colella-Glaz advection scheme (1989) without limiters
  static FieldEmbed<Scal> InterpolateBcg(
      const FieldCell<Scal>& fcu, const MapEmbed<BCond<Scal>>& mebc,
      const FieldEmbed<Scal>& fev, const FieldCell<Scal>& fc_src, const Scal dt,
      const EB& eb);

  template <class T>
  static FieldFace<T> Interpolate(
      const FieldCell<T>& fcu, const MapEmbed<BCond<T>>& mebc, const M& m);
  template <class T>
  static FieldFace<T> Gradient(
      const FieldCell<T>& fcu, const MapEmbed<BCond<T>>& mebc, const M& m);
  static FieldFace<Scal> InterpolateUpwind(
      const FieldCell<Scal>& fcu, const MapEmbed<BCond<Scal>>& mebc, ConvSc sc,
      const FieldCell<Vect>& fcg, const FieldFace<Scal>& ffv, const M& m);
  static FieldFace<Scal> InterpolateBcg(
      const FieldCell<Scal>& fcu, const MapEmbed<BCond<Scal>>& mebc,
      const FieldFace<Scal>& ffv, const FieldCell<Scal>& fc_src, const Scal dt,
      const M& m);

  using ExprFace = typename M::ExprFace;
  using Expr = typename M::Expr;

  // Implicit interpolation with deferred correction.
  // fcu: field cell from previous iteration [s]
  // mebc: boundary conditions
  // sc: interpolation scheme
  // deferred: deferred correction factor (1: fully deferred)
  // fcg: gradient [s]
  // ffv: flow direction [i]
  // Returns:
  // ffe: interpolation expressions [i]
  static FieldEmbed<ExprFace> InterpolateUpwindImplicit(
      const FieldCell<Scal>& fcu, const MapEmbed<BCond<Scal>>& mebc, ConvSc sc,
      Scal deferred, const FieldCell<Vect>& fcg, const FieldEmbed<Scal>& fev,
      const EB& eb);
  static FieldFace<ExprFace> InterpolateUpwindImplicit(
      const FieldCell<Scal>& fcu, const MapEmbed<BCond<Scal>>& mebc, ConvSc sc,
      Scal deferred, const FieldCell<Vect>& fcg, const FieldFace<Scal>& ffv,
      const M& m);

  // Implicit gradient.
  // mebc: boundary conditions
  // Returns:
  // ffe: interpolation expressions [i]
  static FieldEmbed<ExprFace> GradientImplicit(
      const MapEmbed<BCond<Scal>>& mebc, const EB& eb);
  static FieldFace<ExprFace> GradientImplicit(
      const MapEmbed<BCond<Scal>>& mebc, const M& m);

  // Implicit gradient with deferred correction.
  // fcu: field cell from previous iteration [s]
  // mebc: boundary conditions
  // deferred: deferred correction factor (1: fully deferred)
  // Returns:
  // ffe: interpolation expressions [i]
  static FieldEmbed<ExprFace> GradientImplicit(
      const FieldCell<Scal>& fcu, const MapEmbed<BCond<Scal>>& mebc,
      const EB& eb);
  static FieldFace<ExprFace> GradientImplicit(
      const FieldCell<Scal>& fcu, const MapEmbed<BCond<Scal>>& mebc,
      const M& m);

  static FieldCell<Vect> Gradient(const FieldEmbed<Scal>& feu, const EB& eb) {
    return GradientLinearFit(feu, eb);
  }
  static FieldCell<Vect> Gradient(const FieldFace<Scal>& ffu, const M& m);

  static FieldCell<Vect> AverageGradient(
      const FieldEmbed<Scal>& ffg, const EB& eb);
  static FieldCell<Vect> AverageGradient(
      const FieldFace<Scal>& ffg, const M& m);

  template <class MEB>
  static Scal Eval(
      const Expr& e, IdxCell c, const FieldCell<Scal>& fcu, const MEB& meb) {
    Scal r = e.back();
    r += fcu[c] * e[0];
    for (auto q : meb.Nci(c)) {
      r += fcu[meb.GetCell(c, q)] * e[1 + q.raw()];
    }
    return r;
  }
  template <class MEB>
  static Scal Eval(
      const ExprFace& e, IdxFace f, const FieldCell<Scal>& fcu,
      const MEB& meb) {
    const IdxCell cm = meb.GetCell(f, 0);
    const IdxCell cp = meb.GetCell(f, 1);
    return fcu[cm] * e[0] + fcu[cp] * e[1] + e[2];
  }
  static Scal Eval(
      const ExprFace& e, IdxCell c, const FieldCell<Scal>& fcu, const EB&) {
    return fcu[c] * e[0] + e[2];
  }

  template <class MEB>
  static auto InterpolateHarmonic(
      const FieldCell<Scal>& fcu, const MapEmbed<BCond<Scal>>& mebc,
      const MEB& eb) {
    auto inv = fcu;
    for (auto c : eb.AllCells()) {
      inv[c] = 1 / inv[c];
    }
    auto ff = Interpolate(inv, mebc, eb);
    eb.LoopSuFaces([&](auto cf) { //
      ff[cf] = 1 / ff[cf];
    });
    return ff;
  }

  // Computes vorticity of vector field.
  // fcv: vector field [s]
  // mf: boundary conditions for fcv
  // Returns:
  // fco: vorticity [i]
  template <class MEB>
  static FieldCell<Vect> GetVort(
      const FieldCell<Vect>& fcv, const MapEmbed<BCond<Vect>>& mebc,
      const MEB&);
  template <class MEB>
  static FieldCell<Scal> GetVortScal(
      const FieldCell<Vect>& fcv, const MapEmbed<BCond<Vect>>& mebc,
      const MEB&);
};

template <class T, class B>
MapEmbed<BCond<T>> GetBCondZeroGrad(const MapEmbed<B>& mebc) {
  MapEmbed<BCond<T>> r;
  mebc.LoopPairs([&](auto p) {
    auto& bc = r[p.first];
    bc.nci = p.second.nci;
    bc.type = BCondType::neumann;
  });
  return r;
}

// Smoothens a cell field.
// fc: fieldcell [s]
// mfc: boundary conditions
// rep: number of iterations
// Output:
// fc: smooth field [s]
template <class T, class MEB>
void Smoothen(
    FieldCell<T>& fc, const MapEmbed<BCond<T>>& mfc, MEB& eb, size_t iters);

// Converts vector conditions to scalar.
// mfv: vector velocity conditions
// d: direction, 0..2
template <class M>
MapEmbed<BCond<typename M::Scal>> GetScalarCond(
    const MapEmbed<BCond<typename M::Vect>>& mev, size_t d, const M& m);
