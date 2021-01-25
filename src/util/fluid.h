// Created by Petr Karnakov on 04.04.2019
// Copyright 2019 ETH Zurich

#pragma once

#include "geom/mesh.h"
#include "solver/approx_eb.h"
#include "solver/cond.h"
#include "solver/fluid.h"

template <class M_>
class UFluid {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using UEB = UEmbed<M>;

  // fcw: velocity
  // mfc: fluid face conditions
  // fcsv: volume source
  // relax: relaxation factor, 1: full extrapolation
  template <class EB>
  static void UpdateOutletVelocity(
      M& m, const EB& eb, const FieldCell<Vect>& fcvel,
      const MapEmbed<BCondFluid<Vect>>& mebc, const FieldCell<Scal>& fcsv,
      Scal relax, MapEmbed<BCond<Vect>>& mebc_vel) {
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
        switch (bc.type) {
          case BCondFluidType::wall:
          case BCondFluidType::slipwall:
          case BCondFluidType::inlet:
          case BCondFluidType::inletflux:
          case BCondFluidType::outletpressure:
          case BCondFluidType::symm: {
            const Scal q = (nci == 0 ? -1 : 1);
            if (m.IsInner(c)) {
              fluxin += mebc_vel.at(cf).val.dot(eb.GetSurface(cf)) * q;
            }
            break;
          }
          case BCondFluidType::inletpressure: {
            const Scal q = (nci == 0 ? -1 : 1);
            if (m.IsInner(c)) {
              fluxin += std::max(0., fcvel[c].dot(eb.GetSurface(cf)) * q);
            }
            break;
          }
          case BCondFluidType::outlet: {
            const Scal q = (nci == 0 ? 1 : -1);
            Vect vel = fcvel[c] * relax + mebc_vel.at(cf).val * (1 - relax);
            const Vect n = eb.GetNormal(cf);
            Scal vn = vel.dot(n);
            // clip normal component, let only positive
            // (otherwise reversed flow leads to instability)
            vn = (q > 0 ? std::max(0., vn) : std::min(0., vn));
            vel = n * vn + vel.orth(n);
            if (m.IsInner(c)) {
              fluxout += vel.dot(eb.GetSurface(cf)) * q;
              areaout += eb.GetArea(cf);
            }
            mebc_vel.at(cf).val = vel;
            break;
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
      mebc.LoopBCond(eb, [&](auto cf, IdxCell, auto& bc) {
        const auto nci = bc.nci;
        if (bc.type == BCondFluidType::outlet) {
          const Scal q = (nci == 0 ? 1. : -1.);
          const Vect n = eb.GetNormal(cf);
          mebc_vel.at(cf).val += n * (velcor * q);
        }
      });
    }
  }

  // fcw: velocity
  // mfc: fluid face conditions
  // nid: maximum id of inlet flux for reduction
  template <class EB>
  static void UpdateInletFlux(
      M& m, const EB& eb, const FieldCell<Vect>& fcvel,
      const MapEmbed<BCondFluid<Vect>>& mebc, size_t max_id,
      MapEmbed<BCond<Vect>>& mebc_vel) {
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
      flux_target.resize(max_id, 0);
      flux_current.resize(max_id, 0);
      area.resize(max_id, 0);

      // Extrapolate velocity to inlet from neighbor cells
      // and compute total fluxes
      mebc.LoopBCond(eb, [&](auto cf, IdxCell c, auto& bc) {
        const auto nci = bc.nci;
        if (bc.type == BCondFluidType::inletflux) {
          const size_t id = 0; // TODO: implement other id's
          const Scal q = (nci == 0 ? -1. : 1.);
          const Vect vel = fcvel[c].proj(eb.GetNormal(cf));
          if (m.IsInner(c)) {
            flux_target[id] += bc.velocity.dot(eb.GetSurface(cf)) * q;
            flux_current[id] += vel.dot(eb.GetSurface(cf)) * q;
            area[id] += eb.GetArea(cf);
          }
          mebc_vel.at(cf).val = vel;
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
        mebc.LoopBCond(eb, [&](auto cf, IdxCell, auto& bc) {
          const auto nci = bc.nci;
          if (bc.type == BCondFluidType::inletflux) {
            const Scal q = (nci == 0 ? -1. : 1.);
            const Vect n = eb.GetNormal(cf);
            mebc_vel.at(cf).val += n * (velcor * q);
          }
        });
      }
    }
  }

  // fcvel: velocity
  // fcp: pressure
  // mfc: fluid face conditions
  // factor: correction factor, velocity over pressure
  template <class EB>
  static void UpdateInletPressure(
      M& m, const EB& eb, const FieldCell<Scal>& fcp,
      MapEmbed<BCondFluid<Vect>>& mebc, Scal factor) {
    auto sem = m.GetSem("inletpressure");

    struct {
      Scal inlet_pressure = 0;
      Scal inlet_pressure_target = 0;
      Scal inlet_area = 0;
      Scal outlet_pressure = 0;
      Scal outlet_area = 0;
    } * ctx(sem);

    if (sem("local")) {
      auto& t = *ctx;
      mebc.LoopBCond(eb, [&](auto cf, IdxCell c, auto& bc) {
        if (bc.type == BCondFluidType::inletpressure) {
          if (m.IsInner(c)) {
            const auto area = eb.GetArea(cf);
            t.inlet_pressure += fcp[c] * area;
            t.inlet_pressure_target += bc.pressure * area;
            t.inlet_area += area;
          }
        } else if (bc.type == BCondFluidType::outlet) {
          if (m.IsInner(c)) {
            const auto area = eb.GetArea(cf);
            t.outlet_pressure += fcp[c] * area;
            t.outlet_area += area;
          }
        }
      });

      m.Reduce(&t.inlet_area, "sum");
      m.Reduce(&t.inlet_pressure, "sum");
      m.Reduce(&t.inlet_pressure_target, "sum");
      m.Reduce(&t.outlet_area, "sum");
      m.Reduce(&t.outlet_pressure, "sum");
    }

    if (sem("corr")) {
      auto& t = *ctx;
      t.inlet_pressure /= t.inlet_area;
      t.inlet_pressure_target /= t.inlet_area;
      t.outlet_pressure /= t.outlet_area;
      const Scal p_target = t.inlet_pressure_target;
      const Scal p_current = t.inlet_pressure - t.outlet_pressure;
      // additive correction of velocity
      const Scal velcor = (p_target - p_current) * factor;
      mebc.LoopBCond(eb, [&](auto cf, IdxCell, auto& bc) {
        if (bc.type == BCondFluidType::inletpressure) {
          const Scal q = (bc.nci == 0 ? -1. : 1.);
          const Vect n = eb.GetNormal(cf);
          mebc[cf].velocity += n * (velcor * q);
        }
      });
    }
  }

  // fcvel: velocity
  // fcp: pressure
  // mfc: fluid face conditions
  // factor: correction factor, velocity over pressure
  template <class EB>
  static void UpdateVelocityOnPressureBoundaries(
      MapEmbed<BCond<Vect>>& mebc_vel, const EB& eb,
      const FieldCell<Vect>& fcvel, const MapEmbed<BCondFluid<Vect>>& mebc) {
    mebc.LoopBCond(eb, [&](auto cf, IdxCell c, auto& bc) {
      if (bc.type == BCondFluidType::inletpressure) {
        const Vect vel = fcvel[c].proj(eb.GetNormal(cf));
        mebc_vel.at(cf).val = vel;
      } else if (bc.type == BCondFluidType::outletpressure) {
        const Scal q = (bc.nci == 0 ? 1 : -1);
        const Vect n = eb.GetNormal(cf);
        Vect vel = fcvel[c].proj(n);
        Scal vn = vel.dot(n);
        // clip normal component, let only positive
        // (otherwise reversed flow leads to instability)
        vn = (q > 0 ? std::max(0., vn) : std::min(0., vn));
        vel = n * vn + vel.orth(n);
        mebc_vel.at(cf).val = vel;
      }
    });
  }

  // fcvel: velocity
  // fcp: pressure
  // mfc: fluid face conditions
  // factor: correction factor, velocity over pressure
  template <class EB>
  static void UpdateVelocityOnPressureBoundaries(
      MapEmbed<BCond<Vect>>& mebc_vel, M& m, const EB& eb,
      const FieldEmbed<Scal>& fev, const MapEmbed<BCondFluid<Vect>>& mebc,
      Scal relax) {
    auto sem = m.GetSem("inletpressure");
    struct {
      Scal inlet_flux = 0;
      Scal inlet_area = 0;
      Scal outlet_flux = 0;
      Scal outlet_area = 0;
    } * ctx(sem);
    auto& t = *ctx;

    if (sem("reduce")) {
      mebc.LoopBCond(eb, [&](auto cf, IdxCell, auto& bc) {
        if (bc.type == BCondFluidType::inletpressure) {
          t.inlet_flux += fev[cf];
          t.inlet_area += eb.GetArea(cf);
        } else if (bc.type == BCondFluidType::outletpressure) {
          t.outlet_flux += fev[cf];
          t.outlet_area += eb.GetArea(cf);
        }
      });
      m.Reduce(&t.inlet_flux, "sum");
      m.Reduce(&t.inlet_area, "sum");
      m.Reduce(&t.outlet_flux, "sum");
      m.Reduce(&t.outlet_area, "sum");
    }
    if (sem("average")) {
      mebc.LoopBCond(eb, [&](auto cf, IdxCell, auto& bc) {
        if (bc.type == BCondFluidType::inletpressure) {
          const Vect vel = eb.GetNormal(cf) * (t.inlet_flux / t.inlet_area);
          mebc_vel.at(cf).val = vel;
        } else if (bc.type == BCondFluidType::outletpressure) {
          const Scal q = (bc.nci == 0 ? 1 : -1);
          const Vect n = eb.GetNormal(cf);
          const Vect velextrap = eb.GetNormal(cf) * (fev[cf] / eb.GetArea(cf));
          const Vect velavg =
              eb.GetNormal(cf) * (t.outlet_flux / t.outlet_area);
          Vect vel = velextrap * relax + velavg * (1 - relax);
          vel = vel.proj(n);
          Scal vn = vel.dot(n);
          // clip normal component, let only positive
          // (otherwise reversed flow leads to instability)
          vn = (q > 0 ? std::max(0., vn) : std::min(0., vn));
          vel = n * vn + vel.orth(n);
          mebc_vel.at(cf).val = vel;
        }
      });
    }
  }

  template <class MEB>
  static void AppendExplViscous(
      FieldCell<Vect>& fc_force, const FieldCell<Vect>& fc_vel,
      const MapEmbed<BCond<Vect>>& mebc_vel, const FieldFace<Scal>& ff_mu,
      const MEB& eb) {
    fc_force.Reinit(eb);
    auto& m = eb.GetMesh();
    for (auto d : GRange<size_t>(m.GetEdim())) {
      const auto fcu = GetComponent(fc_vel, d);
      const auto mebc = GetScalarCond(mebc_vel, d, m);
      const auto fcg = UEB::AverageGradient(UEB::Gradient(fcu, mebc, eb), eb);
      const auto ffg = UEB::Interpolate(fcg, GetBCondZeroGrad<Vect>(mebc), eb);
      for (auto c : eb.Cells()) {
        if (eb.IsRegular(c)) {
          Vect s(0);
          for (auto q : eb.Nci(c)) {
            const IdxFace f = eb.GetFace(c, q);
            s += ffg[f] * (ff_mu[f] * eb.GetOutwardSurface(c, q)[d]);
          }
          fc_force[c] += s / eb.GetVolume(c);
        }
      }
    }
  }

  template <class MEB>
  static void AppendExplViscousGradMu(
      FieldCell<Vect>& fc_force, const FieldCell<Vect>& fc_vel,
      const MapEmbed<BCond<Vect>>& mebc_vel, const FieldCell<Scal>& fc_mu,
      const MapEmbed<BCond<Scal>>& mebc_mu, const MEB& eb) {
    fc_force.Reinit(eb);
    auto& m = eb.GetMesh();
    const auto fcgm =
        UEB::AverageGradient(UEB::Gradient(fc_mu, mebc_mu, eb), eb);
    for (auto d : GRange<size_t>(m.GetEdim())) {
      const auto fcu = GetComponent(fc_vel, d);
      const auto mebc = GetScalarCond(mebc_vel, d, m);
      const auto fcg = UEB::AverageGradient(UEB::Gradient(fcu, mebc, eb), eb);
      for (auto c : eb.Cells()) {
        if (eb.IsRegular(c)) {
          fc_force[c] += fcg[c] * fcgm[c][d];
        }
      }
    }
  }
};
