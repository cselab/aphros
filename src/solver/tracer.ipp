// Created by Petr Karnakov on 24.06.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <array>
#include <exception>
#include <fstream>
#include <limits>
#include <memory>

#include "approx.h"
#include "approx_eb.h"
#include "util/vof.h"

#include "tracer.h"

template <class EB_>
struct Tracer<EB_>::Imp {
  using Owner = Tracer<EB_>;
  using UEB = UEmbed<M>;

  static Scal Clip(Scal a, Scal low, Scal high) {
    return a < low ? low : a > high ? high : a;
  }
  Imp(Owner* owner, M& m_, const EB& eb_,
      const Multi<const FieldCell<Scal>*>& vfcu,
      const Multi<const MapEmbed<BCond<Scal>>*>& vmebc, Scal time, Conf conf_)
      : owner_(owner)
      , m(m_)
      , eb(eb_)
      , conf(conf_)
      , time_(time)
      , layers(conf.layers)
      , vmebc_(vmebc)
      , vfcu_(vfcu)
      , fc_rho_(m, 0)
      , fc_mu_(m, 0)
      , fc_vf_(m, 0) {}
  Scal GetMixtureViscosity(IdxCell c) const {
    Scal sum = 0;
    for (auto l : layers) {
      sum += vfcu_[l][c] * conf.viscosity[l];
    }
    return sum;
  }
  Scal GetMixtureDensity(IdxCell c) const {
    Scal sum = 0;
    for (auto l : layers) {
      sum += vfcu_[l][c] * conf.density[l];
    }
    return sum;
  }
  Scal GetMixtureVolumeFraction(IdxCell c) const {
    Scal sum = 0;
    for (auto l : layers) {
      sum += vfcu_[l][c];
    }
    return Clip(sum, 0, 1);
  }
  // rho: mixture density
  // mu: mixture viscosity
  Vect GetSlipVelocity(size_t l, Scal rho, Scal mu) const {
    switch (conf.slip[l].type) {
      case SlipType::none:
        return Vect(0);
      case SlipType::stokes:
        return conf.gravity *
               ((conf.density[l] - rho) * sqr(conf.diameter[l]) / (18 * mu));
      case SlipType::constant:
        return conf.slip[l].velocity;
      default:
        fassert(false, "not implemented");
    }
    return Vect(0);
  }
  void Step(Scal dt, const FieldEmbed<Scal>& fev) {
    auto sem = m.GetSem("step");
    struct {
      Multi<FieldCell<Scal>> vfct; // change of conserved quantity
    } * ctx(sem);
    auto& t = *ctx;
    if (sem("local")) {
      t.vfct.Reinit(layers, m, 0);
      for (auto c : eb.AllCells()) {
        fc_rho_[c] = GetMixtureDensity(c);
        fc_mu_[c] = GetMixtureViscosity(c);
        fc_vf_[c] = GetMixtureVolumeFraction(c);
      }
      MapEmbed<BCond<Scal>> me_neumann;
      vmebc_[0].LoopBCond(eb, [&](auto cf, IdxCell, auto bc) { //
        me_neumann[cf] = BCond<Scal>(BCondType::neumann, bc.nci);
      });
      const auto fe_rho = UEmbed<M>::Interpolate(fc_rho_, me_neumann, eb);
      const auto fe_mu = UEmbed<M>::Interpolate(fc_mu_, me_neumann, eb);
      FieldFaceb<Scal> fev_carrier(m, 0);
      eb.LoopFaces([&](auto cf) { //
        fev_carrier[cf] = fev[cf];
      });
      for (auto l : layers) {
        const auto feu = UEmbed<M>::Interpolate(vfcu_[l], vmebc_[l], eb);
        eb.LoopFaces([&, this](auto cf) { //
          auto vel = this->GetSlipVelocity(l, fe_rho[cf], fe_mu[cf]);
          if (vmebc_[l].find(cf)) { // zero flux on boundaries
            vel = Vect(0);
          }
          fev_carrier[cf] -= feu[cf] * vel.dot(eb.GetSurface(cf));
        });
      }

      for (auto l : layers) {
        auto& fcu = vfcu_[l];
        const auto ffg = UEB::Gradient(fcu, vmebc_[l], eb);
        const auto fcg = UEB::AverageGradient(ffg, eb);
        FieldFaceb<Scal> fevl(m); // phase l flux with slip
        eb.LoopFaces([&, this](auto cf) { //
          auto vel = this->GetSlipVelocity(l, fe_rho[cf], fe_mu[cf]);
          if (vmebc_[l].find(cf)) { // zero flux on boundaries
            vel = Vect(0);
          }
          fevl[cf] = fev_carrier[cf] + vel.dot(eb.GetSurface(cf));
        });
        auto feu = UEmbed<M>::InterpolateUpwind(
            fcu, vmebc_[l], conf.scheme, fcg, fevl, eb);

        FieldEmbed<Scal> fe_flux(m, 0);
        eb.LoopFaces([&](auto cf) {
          // advection
          fe_flux[cf] = -feu[cf] * fevl[cf];
          // diffusion
          fe_flux[cf] += conf.diffusion[l] * ffg[cf] * eb.GetArea(cf);
        });
        for (auto c : eb.Cells()) {
          Scal sum = 0;
          eb.LoopNci(c, [&](auto q) {
            const auto cf = eb.GetFace(c, q);
            sum += fe_flux[cf] * eb.GetOutwardFactor(c, q);
          });
          t.vfct[l][c] = dt * sum;
        }
        m.Comm(&t.vfct[l]);
      }
    }
    if (sem("local-redistr")) {
      for (auto l : layers) {
        t.vfct[l] = UEmbed<M>::RedistributeCutCells(t.vfct[l], eb);
        for (auto c : eb.Cells()) {
          vfcu_[l][c] += t.vfct[l][c] / eb.GetVolume(c);
        }
        if (auto* fc_source = conf.fc_source[l]) {
          for (auto c : eb.Cells()) {
            vfcu_[l][c] += dt * (*fc_source)[c];
          }
        }
      }
      if (conf.clip) {
        for (auto c : eb.Cells()) {
          for (auto l : layers) {
            auto& u = vfcu_[l][c];
            u = Clip(u, 0, 1);
          }
        }
      }
      for (auto l : layers) {
        m.Comm(&vfcu_[l]);
      }
    }
    if (!m.flags.fc_innermask.empty() && sem("clear-excluded")) {
      // Clear fields in excluded cells.
      for (auto c : m.AllCells()) {
        if (m.IsExcluded(c)) {
          for (auto l : layers) {
            vfcu_[l][c] = 0;
          }
        }
      }
    }
    if (sem("stat")) {
      time_ += dt;
    }
  }

  Owner* owner_;
  M& m;
  const EB& eb;
  Conf conf;
  Scal time_;
  GRange<size_t> layers;
  Multi<MapEmbed<BCond<Scal>>> vmebc_;
  Multi<FieldCell<Scal>> vfcu_;
  FieldCell<Scal> fc_rho_;
  FieldCell<Scal> fc_mu_;
  FieldCell<Scal> fc_vf_;
};

template <class EB_>
Tracer<EB_>::Tracer(
    M& m, const EB& eb, const Multi<const FieldCell<Scal>*>& vfcu,
    const Multi<const MapEmbed<BCond<Scal>>*>& vmebc, Scal time, Conf conf)
    : imp(new Imp(this, m, eb, vfcu, vmebc, time, conf)) {}

template <class EB_>
Tracer<EB_>::~Tracer() = default;

template <class EB_>
auto Tracer<EB_>::GetConf() const -> const Conf& {
  return imp->conf;
}

template <class EB_>
void Tracer<EB_>::SetConf(Conf conf) {
  imp->conf = conf;
}

template <class EB_>
void Tracer<EB_>::Step(Scal dt, const FieldEmbed<Scal>& fe_flux) {
  imp->Step(dt, fe_flux);
}

template <class EB_>
auto Tracer<EB_>::GetVolumeFraction() const -> const Multi<FieldCell<Scal>>& {
  return imp->vfcu_;
}

template <class EB_>
void Tracer<EB_>::SetVolumeFraction(const Multi<FieldCell<Scal>>& vfcu) {
  imp->vfcu_ = vfcu;
}

template <class EB_>
auto Tracer<EB_>::GetMixtureDensity() const -> const FieldCell<Scal>& {
  return imp->fc_rho_;
}

template <class EB_>
auto Tracer<EB_>::GetMixtureViscosity() const -> const FieldCell<Scal>& {
  return imp->fc_mu_;
}

template <class EB_>
auto Tracer<EB_>::GetMixtureVolumeFraction() const -> const FieldCell<Scal>& {
  return imp->fc_vf_;
}

template <class EB_>
auto Tracer<EB_>::GetView() const -> TracerView {
  return {imp->layers, imp->vfcu_, imp->vmebc_};
}

template <class EB_>
auto Tracer<EB_>::GetTime() const -> Scal {
  return imp->time_;
}

template <class EB_>
auto Tracer<EB_>::GetBCondMutable() -> Multi<MapEmbed<BCond<Scal>>>& {
  return imp->vmebc_;
}
