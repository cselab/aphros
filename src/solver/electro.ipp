// Created by Petr Karnakov on 27.09.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <array>
#include <exception>
#include <fstream>
#include <limits>
#include <memory>

#include "approx.h"
#include "approx_eb.h"
#include "solver/pois.h"
#include "util/vof.h"

#include "electro.h"

template <class EB_>
struct Electro<EB_>::Imp {
  using Owner = Electro<EB_>;
  using UEB = UEmbed<M>;

  Imp(Owner* owner, M& m, const EB& eb, const MapEmbed<BCond<Scal>>& mebc_pot,
      Scal time, Conf conf)
      : owner_(owner)
      , m(m)
      , eb(eb)
      , conf(conf)
      , time_(time)
      , fc_pot_(m, 0)
      , mebc_pot_(mebc_pot) {}
  void Step(
      Scal dt, const FieldCell<Scal>& fc_permit,
      const FieldCell<Scal>& fc_charge) {
    auto sem = m.GetSem("step");
    struct {
      FieldCell<Scal> fc_rhs;
    } * ctx(sem);
    auto& t = *ctx;
    if (sem("local")) {
      t.fc_rhs.Reinit(m, 0);
      for (auto c : eb.CellsM()) {
        t.fc_rhs[c] = std::sin(c.center[0] * M_PI * 2);
      }
    }
    if (sem.Nested("solve")) {
      SolvePoisson(fc_pot_, t.fc_rhs, mebc_pot_, m);
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
  FieldCell<Scal> fc_pot_;
  const MapEmbed<BCond<Scal>>& mebc_pot_;
};

template <class EB_>
Electro<EB_>::Electro(
    M& m, const EB& eb, const MapEmbed<BCond<Scal>>& mebc_pot, Scal time,
    Conf conf)
    : imp(new Imp(this, m, eb, mebc_pot, time, conf)) {}

template <class EB_>
Electro<EB_>::~Electro() = default;

template <class EB_>
auto Electro<EB_>::GetConf() const -> const Conf& {
  return imp->conf;
}

template <class EB_>
void Electro<EB_>::SetConf(Conf conf) {
  imp->conf = conf;
}

template <class EB_>
void Electro<EB_>::Step(
    Scal dt, const FieldCell<Scal>& fc_permit,
    const FieldCell<Scal>& fc_charge) {
  imp->Step(dt, fc_permit, fc_charge);
}

template <class EB_>
auto Electro<EB_>::GetPotential() const -> const FieldCell<Scal>& {
  return imp->fc_pot_;
}

template <class EB_>
auto Electro<EB_>::GetTime() const -> Scal {
  return imp->time_;
}
