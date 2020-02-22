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
          if (auto cd = cb.template Get<Outlet<M>>()) {
            Scal q = (nci == 0 ? 1. : -1.);
            Vect vc = fcw[c];
            Vect s = m.GetSurface(f);
            // clip normal component, let only positive
            vc -= s * (q * std::min(0., vc.dot(s) * q) / s.dot(s));
            cd->SetVelocity(vc);
            fo += cd->GetVelocity().dot(s) * q;
            ao += m.GetArea(f);
          } else if (auto cd = cb.template Get<Inlet<M>>()) {
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

        if (auto cd = cb.template Get<Outlet<M>>()) {
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
  // fcsv: volume source
  template <class EB>
  static void UpdateOutletVelocity(
      M& m, const EB& eb, const FieldCell<Vect>& fcvel,
      const MapEmbed<BCondFluid<Vect>>& mebc, const FieldCell<Scal>& fcsv,
      MapCondFace& mfvel) {
    auto sem = m.GetSem("outlet");
    struct {
      Scal fluxin; // total inlet volume flux
      Scal fluxout; // total outlet volume flux
      Scal areaout; // total outlet area
    } * ctx(sem);
    auto& fluxin = ctx->fluxin;
    auto& fluxout = ctx->fluxout;
    auto& areaout = ctx->areaout;

    if (sem("local")) {
      fluxin = 0.;
      fluxout = 0.;
      areaout = 0.;

      // Extrapolate velocity to outlet from neighbour cells,
      // and compute total fluxes
      for (auto& p : mebc.GetMapFace()) {
        const IdxFace f = p.first;
        const auto& bc = p.second;
        const auto nci = bc.nci;
        const IdxCell c = eb.GetCell(f, nci);
        if (m.IsInner(c)) {
          switch (bc.type) {
            case BCondFluidType::wall:
            case BCondFluidType::slipwall:
            case BCondFluidType::inlet:
            case BCondFluidType::inletflux:
            case BCondFluidType::symm: {
              const Scal q = (nci == 0 ? -1. : 1.);
              fluxin += bc.velocity.dot(eb.GetSurface(f)) * q;
              break;
            }
            case BCondFluidType::outlet: {
              const Scal q = (nci == 0 ? 1. : -1.);
              Vect vel = fcvel[c];
              // clip normal component, let only positive
              // (otherwise reversed flow leads to instability)
              const Vect n = eb.GetNormal(f);
              vel -= n * (q * std::min(0., vel.dot(n) * q));
              mfvel[f].Set<CondFaceValFixed<Vect>>(vel, nci);
              fluxout += vel.dot(eb.GetSurface(f)) * q;
              areaout += eb.GetArea(f);
              break;
            }
          }
        }
      }

      // Append volume source to inlet flux
      for (auto c : eb.Cells()) {
        fluxin += fcsv[c] * m.GetVolume(c);
      }

      m.Reduce(&fluxin, "sum");
      m.Reduce(&fluxout, "sum");
      m.Reduce(&areaout, "sum");
    }

    if (sem("corr")) {
      // additive correction of velocity
      const Scal velcor = (fluxin - fluxout) / areaout;

      // Apply correction on outlet faces
      for (auto& p : mebc.GetMapFace()) {
        const IdxFace f = p.first;
        const auto& bc = p.second;
        const auto nci = bc.nci;
        const IdxCell c = eb.GetCell(f, nci);
        if (m.IsInner(c)) {
          if (bc.type == BCondFluidType::outlet) {
            const Scal q = (nci == 0 ? 1. : -1.);
            const Vect n = eb.GetNormal(f);
            auto cd = mfvel[f].Get<CondFaceValFixed<Vect>>();
            cd->Set(cd->second() + n * (velcor * q));
          }
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
          if (auto cd = cb.template Get<InletFlux<M>>()) {
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

          if (auto cd = cb.template Get<InletFlux<M>>()) {
            size_t nci = cd->GetNci();
            Scal q = (nci == 0 ? -1. : 1.);
            Vect n = m.GetNormal(f);
            cd->SetVelocity(cd->GetVelocity() + n * (dv * q));
          }
        }
      }
    }
  }

  // fcw: velocity
  // mfc: fluid face conditions
  // nid: maximum id of inlet flux for reduction
  template <class EB>
  static void UpdateInletFlux(
      M & m, const EB& eb, const FieldCell<Vect>& fcvel,
      const MapEmbed<BCondFluid<Vect>>& mebc, size_t max_id,
      MapCondFace& mfvel) {
    using namespace fluid_condition;
    auto sem = m.GetSem("inletflux");

    struct {
      std::vector<Scal> flux_target;
      std::vector<Scal> flux_current;
      std::vector<Scal> area;
    } * ctx(sem);
    auto& flux_target = ctx->flux_target;
    auto& flux_current = ctx->flux_current;
    auto& area = ctx->area;

    if (sem("local")) {
      flux_target.resize(max_id);
      flux_current.resize(max_id);
      area.resize(max_id);

      for (size_t id = 0; id < max_id; ++id) {
        flux_target[id] = 0;
        flux_current[id] = 0;
        area[id] = 0;
      }

      // Extrapolate velocity to inlet from neighbor cells
      // and compute total fluxes
      for (auto& p : mebc.GetMapFace()) {
        const IdxFace f = p.first;
        const auto& bc = p.second;
        const auto nci = bc.nci;
        const IdxCell c = eb.GetCell(f, nci);
        if (m.IsInner(c)) {
          if (bc.type == BCondFluidType::inletflux) {
            const size_t id = 0; // TODO: implement other id's
            const Scal q = (nci == 0 ? -1. : 1.);
            flux_target[id] += bc.velocity.dot(eb.GetSurface(f)) * q;
            // extrapolate velocity
            const Vect vel = fcvel[c];
            mfvel[f].Set<CondFaceValFixed<Vect>>(vel, nci);
            flux_current[id] += vel.dot(eb.GetSurface(f)) * q;
            area[id] += eb.GetArea(f);
          }
        }
      }

      for (size_t id = 0; id < max_id; ++id) {
        m.Reduce(&flux_target[id], "sum");
        m.Reduce(&flux_current[id], "sum");
        m.Reduce(&area[id], "sum");
      }
    }

    if (sem("corr")) {
      for (size_t id = 0; id < max_id; ++id) {
        // additive correction of velocity
        const Scal velcor = (flux_target[id] - flux_current[id]) / area[id];
        for (auto& p : mebc.GetMapFace()) {
          const IdxFace f = p.first;
          const auto& bc = p.second;
          const auto nci = bc.nci;
          const IdxCell c = eb.GetCell(f, nci);
          if (m.IsInner(c)) {
            if (bc.type == BCondFluidType::inletflux) {
              const Scal q = (nci == 0 ? -1. : 1.);
              const Vect n = eb.GetNormal(f);
              auto cd = mfvel[f].Get<CondFaceValFixed<Vect>>();
              cd->Set(cd->second() + n * (velcor * q));
            }
          }
        }
      }
    }
  }
};
