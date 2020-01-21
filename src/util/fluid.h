// Created by Petr Karnakov on 04.04.2019
// Copyright 2019 ETH Zurich

#pragma once

#include "geom/mesh.h"
#include "solver/cond.h"
#include "solver/fluid.h"

template <class M_>
class UFluid {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  // fcw: velocity
  // mfc: fluid face conditions
  // fcsv: volume source
  static void UpdateOutletBaseConditions(
      M& m, const FieldCell<Vect>& fcw, MapCondFaceFluid& mfc,
      const FieldCell<Scal>& fcsv) {
    // TODO: Consider seperate channels in one domain
    using namespace fluid_condition;

    auto sem = m.GetSem("outlet");

    struct {
      Scal fi; // total inlet volume flux
      Scal fo; // total outlet volume flux
      Scal ao; // total outlet area
    } * ctx(sem);
    auto& fi = ctx->fi;
    auto& fo = ctx->fo;
    auto& ao = ctx->ao;

    if (sem("local")) {
      fi = 0.;
      fo = 0.;
      ao = 0.;

      // Extrapolate velocity to outlet from neighbour cells,
      // and compute total fluxes
      for (auto& p : mfc) {
        IdxFace f = p.first;
        auto& cb = p.second; // cond base

        size_t nci = cb->GetNci();
        IdxCell c = m.GetCell(f, nci);
        if (m.IsInner(c)) {
          if (auto cd = cb.Get<Outlet<M>>()) {
            Scal q = (nci == 0 ? 1. : -1.);
            Vect vc = fcw[c];
            Vect s = m.GetSurface(f);
            // clip normal component, let only positive
            vc -= s * (q * std::min(0., vc.dot(s) * q) / s.dot(s));
            cd->SetVelocity(vc);
            fo += cd->GetVelocity().dot(s) * q;
            ao += m.GetArea(f);
          } else if (auto cd = cb.Get<Inlet<M>>()) {
            Scal q = (nci == 0 ? -1. : 1.);
            fi += cd->GetVelocity().dot(m.GetSurface(f)) * q;
          }
        }
      }

      // Append volume source to inlet flux
      for (auto c : m.Cells()) {
        fi += fcsv[c] * m.GetVolume(c);
      }

      m.Reduce(&fi, "sum");
      m.Reduce(&fo, "sum");
      m.Reduce(&ao, "sum");
    }

    if (sem("corr")) {
      Scal velcor = (fi - fo) / ao; // Additive correction for velocity

      // Apply correction on outlet faces
      for (auto& it : mfc) {
        IdxFace f = it.first;
        auto& cb = it.second; // cond base

        if (auto cd = cb.Get<Outlet<M>>()) {
          size_t nci = cd->GetNci();
          Scal q = (nci == 0 ? 1. : -1.);
          Vect n = m.GetNormal(f);
          cd->SetVelocity(cd->GetVelocity() + n * (velcor * q));
        }
      }
    }
  }

  // fcw: velocity
  // mfc: fluid face conditions
  // nid: maximum id of inlet flux for reduction
  static void UpdateInletFlux(
      M& m, const FieldCell<Vect>& fcw, MapCondFaceFluid& mfc, size_t nid) {
    // TODO: consider updating from predictor velocity
    using namespace fluid_condition;
    auto sem = m.GetSem("inletflux");

    struct {
      std::vector<Scal> ft; // target flux
      std::vector<Scal> fe; // extrapolated flux
      std::vector<Scal> ai; // inlet area
    } * ctx(sem);
    auto& ft = ctx->ft;
    auto& fe = ctx->fe;
    auto& ai = ctx->ai;

    if (sem("local")) {
      ft.resize(nid);
      fe.resize(nid);
      ai.resize(nid);

      for (size_t id = 0; id < nid; ++id) {
        ft[id] = 0;
        fe[id] = 0;
        ai[id] = 0;
      }

      // Extrapolate velocity to inlet from neighbour cells
      // and compute total fluxes
      for (auto& it : mfc) {
        IdxFace f = it.first;
        auto& cb = it.second;

        size_t nci = cb->GetNci();
        IdxCell c = m.GetCell(f, nci);
        if (m.IsInner(c)) {
          if (auto cd = cb.Get<InletFlux<M>>()) {
            size_t id = cd->GetId();
            assert(id < ft.size());
            Scal q = (nci == 0 ? -1. : 1.);
            // target flux
            ft[id] += cd->GetVelocity().dot(m.GetSurface(f)) * q;
            // extrapolate velocity
            cd->SetVelocity(fcw[c]);
            // extrapolated flux
            fe[id] += cd->GetVelocity().dot(m.GetSurface(f)) * q;
            // area
            ai[id] += m.GetArea(f);
          }
        }
      }

      for (size_t id = 0; id < nid; ++id) {
        m.Reduce(&ft[id], "sum");
        m.Reduce(&fe[id], "sum");
        m.Reduce(&ai[id], "sum");
      }
    }

    if (sem("corr")) {
      for (size_t id = 0; id < nid; ++id) {
        // additive correction of velocity
        Scal dv = (ft[id] - fe[id]) / ai[id];
        for (auto& it : mfc) {
          IdxFace f = it.first;
          auto& cb = it.second;

          if (auto cd = cb.Get<InletFlux<M>>()) {
            size_t nci = cd->GetNci();
            Scal q = (nci == 0 ? -1. : 1.);
            Vect n = m.GetNormal(f);
            cd->SetVelocity(cd->GetVelocity() + n * (dv * q));
          }
        }
      }
    }
  }
};
