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
  static Scal Clip(Scal a) {
    return Clip(a, 0, 1);
  }
  Imp(Owner* owner, M& m, const EB& eb_,
      const Multi<const FieldCell<Scal>*>& vfcu,
      const MapEmbed<BCond<Scal>>& mebc, Scal time, Conf conf_)
      : owner_(owner)
      , m(m)
      , eb(eb_)
      , conf(conf_)
      , time_(time)
      , layers(conf.layers)
      , mebc_(mebc)
      , vfcu_(vfcu)
      , fc_density_(m, conf.density[0])
      , fc_viscosity_(m, conf.viscosity[0]) {
    for (auto c : m.AllCells()) {
      Scal sum = 0;
      for (auto l : layers) {
        if (l > 0) {
          sum += vfcu_[l][c];
        }
      }
      sum = Clip(sum);
      vfcu_[0][c] = 1 - sum;
    }
  }
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
        throw std::runtime_error(FILELINE + ": not implemented");
    }
    return Vect(0);
  }
  void Step(Scal dt, const FieldEmbed<Scal>& fev) {
    auto sem = m.GetSem("step");
    if (sem("local")) {
      auto& fc_rho = fc_density_;
      auto& fc_mu = fc_viscosity_;
      for (auto c : eb.AllCells()) {
        fc_rho[c] = GetMixtureDensity(c);
        fc_mu[c] = GetMixtureViscosity(c);
      }
      const auto fe_rho = UEmbed<M>::Interpolate(fc_rho, mebc_, eb);
      const auto fe_mu = UEmbed<M>::Interpolate(fc_mu, mebc_, eb);
      FieldFaceb<Scal> fev_carrier(m, 0);
      eb.LoopFaces([&](auto cf) { //
        fev_carrier[cf] = fev[cf];
      });
      for (auto l : layers) {
        const auto feu = UEmbed<M>::Interpolate(vfcu_[l], mebc_, eb);
        eb.LoopFaces([&](auto cf) { //
          auto vel = GetSlipVelocity(l, fe_rho[cf], fe_mu[cf]);
          if (mebc_.find(cf)) { // zero flux on boundaries
            vel = Vect(0);
          }
          fev_carrier[cf] -= feu[cf] * vel.dot(eb.GetSurface(cf));
        });
      }

      for (auto l : layers) {
        auto& fcu = vfcu_[l];
        const auto fcg =
            UEB::AverageGradient(UEB::Gradient(fcu, mebc_, eb), eb);
        FieldFaceb<Scal> fevl(m); // phase l flux with slip
        eb.LoopFaces([&](auto cf) { //
          auto vel = GetSlipVelocity(l, fe_rho[cf], fe_mu[cf]);
          if (mebc_.find(cf)) { // zero flux on boundaries
            vel = Vect(0);
          }
          fevl[cf] =
              fev_carrier[cf] + vel.dot(eb.GetSurface(cf));
        });
        auto feu = UEmbed<M>::InterpolateUpwind(
            fcu, mebc_, ConvSc::sou, fcg, fevl, eb);
        // auto feu =
        //    InterpolateSuperbee(fcu, fcg, {}, fev.GetFieldFace(), m);

        FieldEmbed<Scal> fevu(m, 0);
        eb.LoopFaces([&](auto cf) { //
          fevu[cf] = feu[cf] * fevl[cf];
        });
        FieldCell<Scal> fct(eb, 0);
        for (auto c : eb.Cells()) {
          Scal sum = 0.;
          for (auto q : eb.Nci(c)) {
            const auto f = eb.GetFace(c, q);
            sum += fevu[f] * eb.GetOutwardFactor(c, q);
          }
          fct[c] = -dt * sum;
        }
        fct = UEmbed<M>::RedistributeCutCells(fct, eb);
        for (auto c : eb.Cells()) {
          fcu[c] += fct[c] / eb.GetVolume(c);
        }
      }
      // clip to [0, 1] and normalize to sum 1
      for (auto c : eb.Cells()) {
        Scal sum = 0;
        for (auto l : layers) {
          auto& u = vfcu_[l][c];
          u = Clip(u);
          sum += u;
        }
        if (sum > 0) {
          for (auto l : layers) {
            auto& u = vfcu_[l][c];
            u /= sum;
          }
        }
      }
      for (auto l : layers) {
        m.Comm(&vfcu_[l]);
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
  MapEmbed<BCond<Scal>> mebc_;
  Multi<FieldCell<Scal>> vfcu_;
  FieldCell<Scal> fc_density_;
  FieldCell<Scal> fc_viscosity_;
};

template <class EB_>
Tracer<EB_>::Tracer(
    M& m, const EB& eb, const Multi<const FieldCell<Scal>*>& vfcu,
    const MapEmbed<BCond<Scal>>& mebc, Scal time, Conf conf)
    : imp(new Imp(this, m, eb, vfcu, mebc, time, conf)) {}

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
auto Tracer<EB_>::GetMixtureDensity() const -> const FieldCell<Scal>& {
  return imp->fc_density_;
}

template <class EB_>
auto Tracer<EB_>::GetMixtureViscosity() const -> const FieldCell<Scal>& {
  return imp->fc_viscosity_;
}

template <class EB_>
auto Tracer<EB_>::GetView() const -> TracerView {
  return {imp->layers, imp->vfcu_, imp->mebc_};
}

template <class EB_>
auto Tracer<EB_>::GetTime() const -> Scal {
  return imp->time_;
}
