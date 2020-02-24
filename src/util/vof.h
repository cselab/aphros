// Created by Petr Karnakov on 02.09.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <memory>

#include "solver/advection.h"
#include "solver/cond.h"
#include "solver/multi.h"

template <class M_>
class UVof {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  UVof();
  ~UVof();

  // Dumps PLIC polygons from multiple layers.
  // fn: filename
  // t: time
  // bin: binary vtk
  // merge: merge close points
  void DumpPoly(
      const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<const FieldCell<Scal>*>& fccl,
      const Multi<const FieldCell<Vect>*>& fcn,
      const Multi<const FieldCell<Scal>*>& fca,
      const Multi<const FieldCell<bool>*>& fci, std::string fn, Scal t,
      bool bin, bool merge, M& m);

  // Dumps marching cube triangles from multiple layers.
  // fn: filename
  // t: time
  // bin: binary vtk
  // merge: merge close points
  // iso: isovalue for surface fcu=iso
  void DumpPolyMarch(
      const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<const FieldCell<Scal>*>& fccl,
      const Multi<const FieldCell<Vect>*>& fcn,
      const Multi<const FieldCell<Scal>*>& fca,
      const Multi<const FieldCell<bool>*>& fci, std::string fn, Scal t,
      bool bin, bool merge, Scal iso, const FieldCell<Scal>*, M& m);

  // Dumps PLIC polygons from single layer.
  // fn: filename
  // t: time
  // bin: binary vtk
  // merge: merge close points
  void DumpPoly(
      const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcn,
      const FieldCell<Scal>& fca, const FieldCell<bool>& fci, std::string fn,
      Scal t, bool bin, bool merge, M& m);

  // Computes unique color for each connected component over all layers.
  // fcu: volume fraction
  // fccl: color to update
  // fccl0: known colors to keep (may be same as fccl)
  // clfixed: if >=0, override value for color in cell nearest to clfixed_x
  // coalth: merge two layers i,j if  u_i + u_j > coalth
  // mfcu: boundary conditions for u
  // verb: report color overflow (not enough layers)
  // unionfind: use union-find algorithm (otherwise iterative stencil updates)
  // reduce: reduce color space trying to keep the previous color
  static void Recolor(
      const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<FieldCell<Scal>*>& fccl,
      const Multi<const FieldCell<Scal>*>& fccl0, Scal clfixed, Vect clfixed_x,
      Scal coalth, const MapCondFace& mfcu, bool verb, bool unionfind,
      bool reduce, bool grid, M& m);

  static void GetAdvectionFaceCond(
      const M& m, const MapCondFaceAdvection<Scal>& mfc, MapCondFace& mfc_vf,
      MapCondFace& mfc_cl, MapCondFace& mfc_im, MapCondFace& mfc_n,
      MapCondFace& mfc_a);

  static std::tuple<
      MapEmbed<BCond<Scal>>, MapEmbed<BCond<Scal>>, MapEmbed<BCond<Scal>>,
      MapEmbed<BCond<Vect>>, MapEmbed<BCond<Scal>>>
  GetAdvectionBc(const M& m, const MapCondFaceAdvection<Scal>& mfc);

  // set volume fraction to 0 or 1 near wall
  static void BcClear(
      FieldCell<Scal>& uc, const MapCondFaceAdvection<Scal>& mfc, const M& m) {
    for (const auto& it : mfc) {
      auto& cfa = it.second;
      const IdxCell c = m.GetCell(it.first, cfa.GetNci());
      auto& u = uc[c];
      if (u < cfa.clear0) {
        u = 0;
      } else if (u > cfa.clear1) {
        u = 1;
      }
    }
  }

  // set volume fraction to 0 or 1 near wall
  static void BcClearOverrideColor(
      FieldCell<Scal>& uc, FieldCell<Scal>& clc, Scal cl0,
      const MapCondFaceAdvection<Scal>& mfc, const M& m) {
    static constexpr Scal kClNone = -1;
    for (const auto& it : mfc) {
      auto& cfa = it.second;
      const IdxCell c = m.GetCell(it.first, cfa.GetNci());
      auto& u = uc[c];
      if (u < cfa.clear0) {
        u = 0;
      } else if (u > cfa.clear1) {
        u = 1;
      }
      auto& cl = clc[c];
      if (cfa.clear1 < 1 && cl != cl0) {
        cl = kClNone;
        u = 0;
      }
    }
  }

 public:
  struct Imp;
  std::unique_ptr<Imp> imp;
};
