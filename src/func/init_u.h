// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <cmath>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "debug/isnan.h"
#include "dump/hdf.h"
#include "geom/block.h"
#include "geom/field.h"
#include "geom/vect.h"
#include "parse/vars.h"
#include "primlist.h"
#include "solver/multi.h"
#include "solver/reconst.h"

// Volume fraction cut by interface defined by level-set function
// ls: level-set function, interface ls=0, ls>0 for volume fraction 1
// xc: cell center
// h: cell size
template <class Scal, size_t dim>
static Scal GetLevelSetVolume(
    std::function<Scal(const generic::Vect<Scal, dim>&)> ls,
    const generic::Vect<Scal, dim>& xc, const generic::Vect<Scal, dim>& h) {
  using Vect = generic::Vect<Scal, dim>;
  using MIdx = generic::Vect<IntIdx, dim>;

  const Scal dx = 1e-3; // step to compute gradient relative to h

  GBlock<size_t, dim> b(MIdx(2));
  const Scal lsc = ls(xc);
  // traverse nodes of cell
  for (auto w : b) {
    const Vect x = xc + (Vect(w) - Vect(0.5)) * h;
    if ((lsc > 0.) != (ls(x) > 0.)) { // cell crossed by interface
      // linear approximation
      // ls(x) = ls(xc) + (x - xc).dot(grad(ls, xc))
      Vect n; // normal
      for (size_t i = 0; i < dim; ++i) {
        Vect xp(xc), xm(xc);
        const Scal dxh = dx * h[i];
        xp[i] += dxh * 0.5;
        xm[i] -= dxh * 0.5;
        n[i] = (ls(xp) - ls(xm)) / dxh;
      }
      return Reconst<Scal>::GetLineU(n, lsc, h);
    }
  }

  return lsc > 0. ? 1. : 0.;
}

// Returns point at which interpolant has value 0.
// x0,x1: points
// f0,f1: values
template <class Scal, class Vect>
Vect GetIso(Vect x0, Vect x1, Scal f0, Scal f1) {
  return (x0 * f1 - x1 * f0) / (f1 - f0);
}

template <class Scal>
Scal GetPositiveEdgeFraction(Scal l0, Scal l1) {
  if (l0 * l1 < 0) {
    return l0 < l1 ? l1 / (l1 - l0) : l0 / (l0 - l1);
  }
  return l0 + l1 > 0;
}

// Returns fraction of face area for which ls>0.
// ee: fraction of edge length for which ls>0
//       3
//   -------
//  0|     |
//   |     | 2
//   -------
//     1
template <class Scal>
Scal GetFaceAreaFraction(std::array<Scal, 4> ee) {
  using R = Reconst<Scal>;
  // face center is x=0 y=0
  // line equation
  //   n.dot(x) = a
  if (ee[0] < ee[2]) std::swap(ee[0], ee[2]);
  if (ee[1] < ee[3]) std::swap(ee[1], ee[3]);

  constexpr size_t ni = 4;

  const std::array<Scal, ni> xx{-0.5, ee[1] - 0.5, 0.5, ee[3] - 0.5};
  const std::array<Scal, ni> yy{ee[0] - 0.5, -0.5, ee[2] - 0.5, 0.5};

  // normal
  Scal nx = ee[0] - ee[2];
  Scal ny = ee[1] - ee[3];

  Scal a = 0;
  size_t aw = 0;
  for (size_t i = 0; i < ni; ++i) {
    const auto& e = ee[i];
    if (e > 0 && e < 1) {
      a += xx[i] * nx + yy[i] * ny;
      ++aw;
    }
  }
  if (aw == 0 || nx + ny == 0) {
    return (ee[0] + ee[1] + ee[2] + ee[3]) / ni;
  }
  a /= aw;

  if (nx > ny) {
    std::swap(nx, ny);
  }
  using Vect2 = generic::Vect<Scal, 2>;
  if (a < 0) {
    return R::GetLineU0(Vect2{nx, ny}, a);
  }
  return 1 - R::GetLineU0(Vect2{nx, ny}, -a);
}

// Computes face area for which ls > 0.
// fnl: level-set function ls on nodes [i]
template <class M>
FieldFace<typename M::Scal> GetPositiveAreaFraction(
    const FieldNode<typename M::Scal>& fnl, const M& m) {
  using Scal = typename M::Scal;
  FieldFace<Scal> ffs(m, 0);
  for (auto f : m.Faces()) {
    constexpr size_t em = 4;
    std::array<Scal, em> ll;
    for (size_t e = 0; e < em; ++e) {
      const size_t ep = (e + 1) % em;
      const IdxNode n = m.GetNode(f, e);
      const IdxNode np = m.GetNode(f, ep);
      ll[e] = GetPositiveEdgeFraction(fnl[n], fnl[np]);
      ffs[f] = GetFaceAreaFraction(ll);
    }
  }
  return ffs;
}

// fnl: level-set function ls computed on nodes [i]
template <class M>
FieldCell<typename M::Scal> GetPositiveVolumeFraction(
    const FieldNode<typename M::Scal>& fnl, const M& m) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using R = Reconst<Scal>;
  auto ffs = GetPositiveAreaFraction(fnl, m);
  FieldCell<Scal> fcu(m, 0);
  for (auto c : m.Cells()) {
    Vect nn(0); // normal
    for (auto q : m.Nci(c)) {
      const IdxFace f = m.GetFace(c, q);
      nn += m.GetNormal(f) * ffs[f] * m.GetOutwardFactor(c, q);
    }
    Scal u = 0;
    if (nn.norm1() == 0) {
      for (auto q : m.Nci(c)) {
        u += ffs[m.GetFace(c, q)];
      }
      u /= m.Nci(c).size();
    } else {
      nn /= -nn.norm1();
      Scal a = 0; // plane constant
      size_t aw = 0;
      for (auto q : m.Nci(c)) {
        const IdxFace f = m.GetFace(c, q);
        const size_t em = m.GetNumNodes(f);
        for (size_t e = 0; e < em; ++e) {
          const size_t ep = (e + 1) % em; // next node
          const IdxNode n = m.GetNode(f, e);
          const IdxNode np = m.GetNode(f, ep);
          const Scal l = fnl[n];
          const Scal lp = fnl[np];
          const Vect x = m.GetNode(n);
          const Vect xp = m.GetNode(np);

          if (l * lp < 0) {
            a += nn.dot(GetIso(x, xp, l, lp) - m.GetCenter(c));
            aw += 1;
          }
        }
      }
      assert(aw > 0);
      a /= aw;
      u = R::GetLineU(nn, a, m.GetCellSize());
    }
    fcu[c] = u;
  }
  return fcu;
}

// Fills volume fraction field from list of primitives.
// fc: field to fill
// primlist: list of primitives
// edim: effective dimension
// approx: 0: stepwise, 1: level-set, 2: overlap, 3: level-set on nodes
template <class M>
void InitVfList(
    FieldCell<typename M::Scal>& fc, std::istream& primlist, int approx,
    size_t edim, const M& m, bool verbose);

// Returns an initialized of the volume fraction field.
// par: parameters
// M: mesh
// Returns:
// std::function<void(GField<Cell>& fc,const M& m)>
// fc: field to fill [i]
// m: mesh
template <class M>
std::function<void(FieldCell<typename M::Scal>&, const M&)> CreateInitU(
    const Vars& par, bool verb);

// Reads list of primitives from either file `list_path`
// or inline if "inline" is the first word of `list_path`
std::stringstream ReadPrimList(std::string list_path, bool verbose);

template <class M>
void InitVf(
    FieldCell<typename M::Scal>& fcu, const Vars& var, M& m, bool verbose) {
  auto sem = m.GetSem("initvf");
  struct {
    std::vector<char> buf;
  } * ctx(sem);

  if (sem("zero")) {
    fcu.Reinit(m, 0);
  }
  const std::string init_vf = var.String["init_vf"];
  if (init_vf == "list") {
    if (sem("list-bcast")) {
      if (m.IsRoot()) {
        auto buf = ReadPrimList(var.String["list_path"], verbose);
        ctx->buf = std::vector<char>(
            std::istreambuf_iterator<char>(buf),
            std::istreambuf_iterator<char>());
      }
      m.Bcast(&ctx->buf);
    }
    if (sem("list-local")) {
      std::stringstream list;
      std::copy(
          ctx->buf.begin(), ctx->buf.end(), std::ostream_iterator<char>(list));
      InitVfList(fcu, list, var.Int["list_ls"], var.Int["dim"], m, verbose);
    }
  } else if (init_vf == "hdf") {
    if (sem.Nested()) {
      Hdf<M>::Read(fcu, var.String["init_vf_hdf_path"], m);
    }
  } else {
    if (sem("local")) {
      auto func = CreateInitU<M>(var, m.IsRoot());
      func(fcu, m);
    }
  }
  if (sem("comm")) {
    m.Comm(&fcu);
  }
}

// listpath: path to list of primitives b.dat
// fcu: volume fraction field
// fccl: color field
template <class M>
void InitOverlappingComponents(
    std::istream& primlist, const Multi<FieldCell<typename M::Scal>*>& fcu,
    const Multi<FieldCell<typename M::Scal>*>& fccl,
    const GRange<size_t>& layers, const M& m);
