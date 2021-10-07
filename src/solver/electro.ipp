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
#include "linear/linear.h"
#include "solver/pois.h"
#include "util/vof.h"

#include "electro.h"

template <class EB_>
struct Electro<EB_>::Imp {
  using Owner = Electro<EB_>;
  using UEB = UEmbed<M>;

  Imp(Owner* owner, M& m_, const EB& eb_, const MapEmbed<BCond<Scal>>& mebc_pot,
      Scal time, const Conf& conf_)
      : owner_(owner)
      , m(m_)
      , eb(eb_)
      , conf(conf_)
      , time_(time)
      , fc_pot_(m, 0)
      , fc_current_(m, Vect(0))
      , ff_current_(m, 0)
      , mebc_pot_(mebc_pot) {}
  void Step(Scal dt, const FieldCell<Scal>& fc_vf) {
    auto sem = m.GetSem("step");
    using Expr = typename M::Expr;
    using ExprFace = typename M::ExprFace;
    struct {
      FieldCell<Scal> fc_rhs;
      FieldCell<Expr> fcl;
      FieldFaceb<Scal> ff_resist;
    } * ctx(sem);
    auto& t = *ctx;
    if (sem("local")) {
      t.fc_rhs.Reinit(m, 0);
    }
    if (sem("solve")) {
      const FieldFaceb<Scal> ff_vf =
          UEB::Interpolate(fc_vf, GetBCondZeroGrad<Scal>(mebc_pot_), eb);
      auto r1 = conf.var.Double["resist1"];
      auto r2 = conf.var.Double["resist2"];
      t.ff_resist.Reinit(m);
      eb.LoopFaces([&](auto cf) { //
        t.ff_resist[cf] = 1 / (1 / r2 * ff_vf[cf] + 1 / r1 * (1 - ff_vf[cf]));
      });

      // FIXME: Don't rely on `time_`.
      //        Zero initial guess causes wrong dirichlet boundary conditions
      //        at first step if second order.
      FieldFaceb<ExprFace> ffg;
      if (time_ == 0) {
        ffg = UEB::GradientImplicit(mebc_pot_, eb);
      } else {
        ffg = UEB::GradientImplicit(fc_pot_, mebc_pot_, eb);
      }
      t.fcl.Reinit(m, Expr::GetUnit(0));
      for (auto c : eb.Cells()) {
        if (eb.IsExcluded(c)) continue;
        Expr sum(0);
        eb.LoopNci(c, [&](auto q) {
          const auto cf = eb.GetFace(c, q);
          const ExprFace flux = ffg[cf] / t.ff_resist[cf] * eb.GetArea(cf);
          eb.AppendExpr(sum, flux * eb.GetOutwardFactor(c, q), q);
        });
        t.fcl[c] = sum;
      }
      t.fcl.SetName("electro");
    }
    if (sem.Nested("solve")) {
      conf.linsolver->Solve(t.fcl, &fc_pot_, fc_pot_, m);
    }
    if (!m.flags.fc_innermask.empty() && sem("clear-excluded")) {
      // Clear field in excluded cells.
      for (auto c : m.AllCells()) {
        if (m.IsExcluded(c)) {
          fc_pot_[c] = 0;
        }
      }
    }
    if (sem("post")) {
      time_ += dt;

      const auto ffg = UEB::Gradient(fc_pot_, mebc_pot_, eb);
      ff_current_.Reinit(m, 0);
      eb.LoopFaces([&](auto cf) { //
        ff_current_[cf] = ffg[cf] / t.ff_resist[cf];
      });
      fc_current_ = UEB::AverageGradient(ff_current_, eb);

      stat_.current = 0;
      stat_.potential_min = std::numeric_limits<Scal>::max();
      stat_.potential_max = -std::numeric_limits<Scal>::max();
      mebc_pot_.LoopBCond(eb, [&](auto cf, IdxCell c, auto& bc) {
        const auto nci = bc.nci;
        if (m.IsInner(c)) {
          const Scal current = ff_current_[cf] * (nci == 0 ? -1 : 1);
          if (current < 0) {
            stat_.current -= current * eb.GetArea(cf);
          }
          const Scal pot = fc_pot_[c];
          if (pot < stat_.potential_min) {
            stat_.potential_min = std::min(stat_.potential_min, pot);
          }
          if (pot > stat_.potential_max) {
            stat_.potential_max = std::max(stat_.potential_max, pot);
          }
        }
      });
      m.Reduce(&stat_.current, "sum");
      m.Reduce(&stat_.potential_min, "min");
      m.Reduce(&stat_.potential_max, "max");
    }
    if (sem()) {
      stat_.potential = stat_.potential_max - stat_.potential_min;
    }
  }

  Owner* owner_;
  M& m;
  const EB& eb;
  Conf conf;
  Scal time_;
  FieldCell<Scal> fc_pot_;
  FieldCell<Vect> fc_current_;
  FieldEmbed<Scal> ff_current_;
  const MapEmbed<BCond<Scal>>& mebc_pot_;
  Stat stat_;
};

template <class EB_>
Electro<EB_>::Electro(
    M& m, const EB& eb, const MapEmbed<BCond<Scal>>& mebc_pot, Scal time,
    const Conf& conf)
    : imp(new Imp(this, m, eb, mebc_pot, time, conf)) {}

template <class EB_>
Electro<EB_>::~Electro() = default;

template <class EB_>
auto Electro<EB_>::GetConf() const -> Conf& {
  return imp->conf;
}

template <class EB_>
void Electro<EB_>::Step(Scal dt, const FieldCell<Scal>& fc_vf) {
  imp->Step(dt, fc_vf);
}

template <class EB_>
auto Electro<EB_>::GetPotential() const -> const FieldCell<Scal>& {
  return imp->fc_pot_;
}

template <class EB_>
auto Electro<EB_>::GetCurrent() const -> const FieldCell<Vect>& {
  return imp->fc_current_;
}

template <class EB_>
auto Electro<EB_>::GetFaceCurrent() const -> const FieldEmbed<Scal>& {
  return imp->ff_current_;
}

template <class EB_>
auto Electro<EB_>::GetStat() const -> const Stat& {
  return imp->stat_;
}
template <class EB_>
auto Electro<EB_>::GetTime() const -> Scal {
  return imp->time_;
}
