// Created by Petr Karnakov on 26.06.2020
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

#include "particles.h"

template <class EB_>
struct Particles<EB_>::Imp {
  using Owner = Particles<EB_>;
  using UEB = UEmbed<M>;

  struct State {
    std::vector<Vect> x;
    std::vector<Vect> v;
    std::vector<Scal> r;
    std::vector<Scal> rho;
  };

  Imp(Owner* owner, M& m, const EB& eb_, const ParticlesView& init, Scal time,
      Conf conf_)
      : owner_(owner)
      , m(m)
      , eb(eb_)
      , conf(conf_)
      , time_(time)
      , state_({init.x, init.v, init.r, init.rho}) {}
  void Step(Scal dt, const FieldEmbed<Scal>& fev) {
    auto sem = m.GetSem("step");
    auto& s = state_;
    if (sem("local")) {
      for (size_t i = 0; i < s.x.size(); ++i) {
        s.v[i] += conf.gravity * dt;
        s.x[i] += s.v[i] * dt;
      }
    }
    if (sem.Nested()) {
      Comm(s.x, {}, {}, m);
    }
    if (sem("stat")) {
      time_ += dt;
    }
  }
  // Exchanges particle positions and data between blocks
  // such that the positions are inside the owning block.
  // x: particle positions
  // attr_scal: scalar attributes
  // attr_vect: vector attributes
  static void Comm(
      std::vector<Vect>& x, const std::vector<std::vector<Scal>*>& attr_scal,
      const std::vector<std::vector<Vect>*>& attr_vect, M& m) {
    auto sem = m.GetSem("particles-comm");
    struct {
      std::vector<Vect> send_x;
      std::vector<std::vector<Scal>> send_attr_scal;
      std::vector<std::vector<Vect>> send_attr_vect;
    } * ctx(sem);
    auto& send_x = ctx->send_x;
    auto& send_attr_scal = ctx->send_attr_scal;
    auto& send_attr_vect = ctx->send_attr_vect;
    if (sem("local")) {
      std::vector<size_t> outside;
      const auto box = m.GetBoundingBox();
      for (size_t i = 0; i < x.size(); ++i) {
        if (!box.IsInside(x[i])) {
          outside.push_back(i);
        }
      }
      for (auto i : outside) {
        send_x.push_back(x[i]);
      }
      for (auto& attr : attr_scal) {
        send_attr_scal.emplace_back();
        for (auto i : outside) {
          send_attr_scal.back().push_back((*attr)[i]);
        }
      }
      for (auto& attr : attr_vect) {
        send_attr_vect.emplace_back();
        for (auto i : outside) {
          send_attr_vect.back().push_back((*attr)[i]);
        }
      }
    }
    if (sem()) {
    }
  }
  void DumpCsv(std::string path) const {
    auto sem = m.GetSem("dumpcsv");
    struct {
      State s;
    } * ctx(sem);
    if (sem("gather")) {
      auto& s = state_;
      fassert(s.x.size() == s.x.size());
      fassert(s.r.size() == s.x.size());
      fassert(s.rho.size() == s.x.size());
      ctx->s.x = s.x;
      ctx->s.v = s.v;
      ctx->s.r = s.r;
      ctx->s.rho = s.rho;
      using TV = typename M::template OpCatT<Vect>;
      m.Reduce(std::make_shared<TV>(&ctx->s.x));
      m.Reduce(std::make_shared<TV>(&ctx->s.v));
      using TS = typename M::template OpCatT<Scal>;
      m.Reduce(std::make_shared<TS>(&ctx->s.r));
      m.Reduce(std::make_shared<TS>(&ctx->s.rho));
    }
    if (sem("init")) {
      if (m.IsRoot()) {
        std::ofstream o(path);
        o.precision(16);
        // header
        o << "x,y,z,vx,vy,vz,r";
        o << std::endl;
        // content
        auto& s = ctx->s;
        for (size_t i = 0; i < s.x.size(); ++i) {
          o << s.x[i][0] << ',' << s.x[i][1] << ',' << s.x[i][2];
          o << ',' << s.v[i][0] << ',' << s.v[i][1] << ',' << s.v[i][2];
          o << ',' << s.r[i];
          o << "\n";
        }
      }
    }
    if (sem()) {
    }
  }

  Owner* owner_;
  M& m;
  const EB& eb;
  Conf conf;
  Scal time_;
  State state_;
};

template <class EB_>
Particles<EB_>::Particles(
    M& m, const EB& eb, const ParticlesView& init, Scal time, Conf conf)
    : imp(new Imp(this, m, eb, init, time, conf)) {}

template <class EB_>
Particles<EB_>::~Particles() = default;

template <class EB_>
auto Particles<EB_>::GetConf() const -> const Conf& {
  return imp->conf;
}

template <class EB_>
void Particles<EB_>::SetConf(Conf conf) {
  imp->conf = conf;
}

template <class EB_>
void Particles<EB_>::Step(Scal dt, const FieldEmbed<Scal>& fe_flux) {
  imp->Step(dt, fe_flux);
}

template <class EB_>
auto Particles<EB_>::GetView() const -> ParticlesView {
  return {imp->state_.x, imp->state_.v, imp->state_.r, imp->state_.rho};
}

template <class EB_>
void Particles<EB_>::DumpCsv(std::string path) const {
  imp->DumpCsv(path);
}

template <class EB_>
auto Particles<EB_>::GetTime() const -> Scal {
  return imp->time_;
}
