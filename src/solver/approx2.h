// Created by Petr Karnakov on 30.11.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <array>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "approx_eb.h"
#include "cond.h"
#include "embed.h"
#include "geom/mesh.h"
#include "solver.h"

// This file contains another set of numerical routines,
// complementary to `approx.h` and `approx_eb.h`.
template <class MEB>
struct Approx2 {
  using M = typename MEB::M;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  template <class T>
  using FieldFaceb = typename EmbedTraits<MEB>::template FieldFaceb<T>;

  // Returns divergence of face field
  // defined as sum of fluxes divided by volume of regular cell.
  static FieldCell<Scal> GetRegularDivergence(
      const FieldFaceb<Scal>& fev, const MEB& eb) {
    auto& m = eb.GetMesh();
    FieldCell<Scal> fcdiv(m, 0);
    for (auto c : eb.Cells()) {
      Scal div = 0;
      eb.LoopNci(c, [&](auto q) {
        const auto cf = eb.GetFace(c, q);
        div += fev[cf] * eb.GetOutwardFactor(c, q);
      });
      fcdiv[c] = div / m.GetVolume(c);
    }
    return fcdiv;
  }

  // Extrapolates values from inner faces to halo faces using boundary
  // conditions for the corresponding vector field.
  // The face field contains the normal components of the cell vector field.
  // ffv: scalar face field
  // mebc: boundary conditions
  static void ExtrapolateToHaloFaces(
      FieldFace<Scal>& ffv, const MapEmbed<BCond<Vect>>& mebc, M& m) {
    auto sem = m.GetSem(__func__);
    if (sem.Nested("comm")) {
      CommFieldFace(ffv, m);
    }
    if (sem("extrapolate")) {
      // Assign NaN in halo faces adjacent to faces with boundary conditions.
      for (const auto& p : mebc.GetMapFace()) {
        const auto f = m(p.first);
        const auto& bc = p.second;
        const auto c = f.cell(1 - bc.nci); // Outer (halo) adjacent cell.
        const auto d = f.direction();
        // FIXME: Specific for 3D.
        const auto du = d.next(1);
        const auto dv = d.next(2);
        ffv[c.face(du)] = GetNan<Scal>();
        ffv[c.face(-du)] = GetNan<Scal>();
        ffv[c.face(dv)] = GetNan<Scal>();
        ffv[c.face(-dv)] = GetNan<Scal>();
      }
      // Extrapolate to halo faces using boundary conditions.
      // Most faces are encountered twice.
      for (const auto& p : mebc.GetMapFace()) {
        const auto f = m(p.first);
        const auto& bc = p.second;
        const auto co = f.cell(1 - bc.nci); // Outer (halo) adjacent cell.
        const auto ci = f.cell(bc.nci); // Inner adjacent cell.
        const auto d = f.direction();
        // FIXME: Specific for 3D.
        const auto du = d.next(1);
        const auto dv = d.next(2);
        auto update = [&](IdxFace fo, Scal veln_bc, IdxFace fi) {
          auto& v = ffv[fo];
          if (IsNan(v)) {
            v = veln_bc - ffv[fi];
          } else {
            v = (v + veln_bc - ffv[fi]) * 0.5;
          }
        };
        switch (bc.type) {
          case BCondType::dirichlet:
            update(co.face(du), bc.val[du], ci.face(du));
            update(co.face(-du), bc.val[du], ci.face(-du));
            update(co.face(dv), bc.val[dv], ci.face(dv));
            update(co.face(-dv), bc.val[dv], ci.face(-dv));
            break;
          // TODO: implement other BCondType
          default:
            ffv[co.face(du)] = 0;
            ffv[co.face(-du)] = 0;
            ffv[co.face(dv)] = 0;
            ffv[co.face(-dv)] = 0;
            break;
        }
      }
    }
    if (sem()) {
      // XXX: empty stage
    }
  }

  // Evaluates a trilinear interpolant constructed from a scalar face field.
  // Assumes that face values adjacent to halo cells are filled.
  // The scalar face field stores normal components of a vector field.
  // The interpolant takes a point and returns a vector.
  // Points must be inside inner cells.
  // ffu: face scalar field
  // callback: function that takes the interpolant function as a parameter
  static void EvalTrilinearFromFaceField(
      const FieldFace<Scal>& ffu,
      const std::function<void(const std::function<Vect(Vect x)>&)>& callback,
      const M& m) {
    auto func = [&m, &ffu](Vect x) -> Vect {
      const Vect h = m.GetCellSize();
      const auto c = m(m.GetCellFromPoint(x));
      Vect res;
      for (auto di : m.dirs) {
        const Vect delta = (x - c.center) / h;
        const auto d = m.direction(di);
        const auto du = d.next(1).orient(delta);
        Scal vm, vp;
        if (M::dim == 2) {
          const IdxCellMesh<M> cc[] = {c, c + du};
          const IdxFace fm[] = {cc[0].face(-d), cc[1].face(-d)};
          const IdxFace fp[] = {cc[0].face(d), cc[1].face(d)};
          vm = interp::Linear(std::abs(delta[du]), ffu[fm[0]], ffu[fm[1]]);
          vp = interp::Linear(std::abs(delta[du]), ffu[fp[0]], ffu[fp[1]]);
        } else {
          const auto dv = d.next(2).orient(delta);
          const IdxCellMesh<M> cc[] = {c, c + du, c + dv, c + du + dv};
          const IdxFace fm[] = {
              cc[0].face(-d), cc[1].face(-d), cc[2].face(-d), cc[3].face(-d)};
          const IdxFace fp[] = {
              cc[0].face(d), cc[1].face(d), cc[2].face(d), cc[3].face(d)};
          vm = interp::Bilinear(
              std::abs(delta[du]), std::abs(delta[dv]), //
              ffu[fm[0]], ffu[fm[1]], ffu[fm[2]], ffu[fm[3]]);
          vp = interp::Bilinear(
              std::abs(delta[du]), std::abs(delta[dv]), //
              ffu[fp[0]], ffu[fp[1]], ffu[fp[2]], ffu[fp[3]]);
        }
        res[di] = interp::Linear(delta[d] + 0.5, vm, vp);
      }
      return res;
    };
    callback(func);
  }

};
