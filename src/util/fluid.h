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
  template <class EB>
  static void UpdateOutletVelocity(
      M& m, const EB& eb, const FieldCell<Vect>& fcvel,
      const MapEmbed<BCondFluid<Vect>>& mebc, const FieldCell<Scal>& fcsv,
      MapEmbed<BCond<Vect>>& mebc_vel) {
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
      mebc.LoopBCond(eb, [&](auto cf, IdxCell c, auto& bc) {
        const auto nci = bc.nci;
        if (m.IsInner(c)) {
          switch (bc.type) {
            case BCondFluidType::wall:
            case BCondFluidType::slipwall:
            case BCondFluidType::inlet:
            case BCondFluidType::inletflux:
            case BCondFluidType::symm: {
              const Scal q = (nci == 0 ? -1. : 1.);
              fluxin += bc.velocity.dot(eb.GetSurface(cf)) * q;
              break;
            }
            case BCondFluidType::outlet: {
              const Scal q = (nci == 0 ? 1. : -1.);
              Vect vel = fcvel[c];
              const Vect n = eb.GetNormal(cf);
              // clip normal component, let only positive
              // (otherwise reversed flow leads to instability)
              vel -= n * (q * std::min(0., vel.dot(n) * q));
              fluxout += vel.dot(eb.GetSurface(cf)) * q;
              areaout += eb.GetArea(cf);
              mebc_vel.at(cf).val = vel;
              break;
            }
          }
        }
      });

      // Append volume source to inlet flux
      for (auto c : eb.Cells()) {
        fluxin += fcsv[c] * eb.GetVolume(c);
      }

      m.Reduce(&fluxin, "sum");
      m.Reduce(&fluxout, "sum");
      m.Reduce(&areaout, "sum");
    }

    if (sem("corr")) {
      // additive correction of velocity
      const Scal velcor = (fluxin - fluxout) / areaout;

      // Apply correction on outlet faces
      mebc.LoopBCond(eb, [&](auto cf, IdxCell c, auto& bc) {
        const auto nci = bc.nci;
        if (m.IsInner(c)) {
          if (bc.type == BCondFluidType::outlet) {
            const Scal q = (nci == 0 ? 1. : -1.);
            const Vect n = eb.GetNormal(cf);
            mebc_vel.at(cf).val += n * (velcor * q);
          }
        }
      });
    }
  }

  // fcw: velocity
  // mfc: fluid face conditions
  // nid: maximum id of inlet flux for reduction
  template <class EB>
  static void UpdateInletFlux(
      M & m, const EB& eb, const FieldCell<Vect>& fcvel,
      const MapEmbed<BCondFluid<Vect>>& mebc, size_t max_id,
      MapEmbed<BCond<Vect>>& mebc_vel) {
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
      mebc.LoopBCond(eb, [&](auto cf, IdxCell c, auto& bc) {
        const auto nci = bc.nci;
        if (m.IsInner(c)) {
          if (bc.type == BCondFluidType::inletflux) {
            const size_t id = 0; // TODO: implement other id's
            const Scal q = (nci == 0 ? -1. : 1.);
            flux_target[id] += bc.velocity.dot(eb.GetSurface(cf)) * q;
            // extrapolate velocity
            const Vect vel = fcvel[c];
            flux_current[id] += vel.dot(eb.GetSurface(cf)) * q;
            area[id] += eb.GetArea(cf);
            mebc_vel.at(cf).val = vel;
          }
        }
      });

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
        mebc.LoopBCond(eb, [&](auto cf, IdxCell c, auto& bc) {
          const auto nci = bc.nci;
          if (m.IsInner(c)) {
            if (bc.type == BCondFluidType::inletflux) {
              const Scal q = (nci == 0 ? -1. : 1.);
              const Vect n = eb.GetNormal(cf);
              mebc_vel.at(cf).val += n * (velcor * q);
            }
          }
        });
      }
    }
  }
};
