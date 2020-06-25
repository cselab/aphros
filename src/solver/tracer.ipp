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
      , fc_density_(m, GetNan<Scal>())
      , fc_viscosity_(m, GetNan<Scal>()) {}
  void Step(Scal dt, const FieldEmbed<Scal>& fe_flux) {
    auto sem = m.GetSem("step");
    if (sem("init")) {
      for (auto l : layers) {
        for (auto c : m.Cells()) {
          vfcu_[l][c] += 0.01;
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
