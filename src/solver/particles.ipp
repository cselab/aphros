// Created by Petr Karnakov on 26.06.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <array>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <type_traits>

#include "approx.h"
#include "approx_eb.h"
#include "util/format.h"
#include "util/vof.h"

#include "particles.h"

template <class EB_>
struct Particles<EB_>::Imp {
  using Owner = Particles<EB_>;
  using UEB = UEmbed<M>;

  struct State {
    // See description of attributes in ParticlesView.
    std::vector<Vect> x;
    std::vector<bool> is_inner;
    std::vector<Vect> v;
    std::vector<Scal> r;
    std::vector<Scal> source;
    std::vector<Scal> rho;
    std::vector<Scal> termvel;
    std::vector<Scal> removed;
  };

  Imp(Owner* owner, M& m_, const EB& eb_, const ParticlesView& init, Scal time,
      Conf conf_)
      : owner_(owner)
      , m(m_)
      , eb(eb_)
      , conf(conf_)
      , time_(time)
      , state_({init.x, init.is_inner, init.v, init.r, init.source, init.rho,
                init.termvel, init.removed}) {}
  static ParticlesView GetView(State& s) {
    return {s.x, s.is_inner, s.v, s.r, s.source, s.rho, s.termvel, s.removed};
  }
  static void CheckSize(const State& s) {
    const size_t n = s.x.size();
    fassert_equal(s.is_inner.size(), n);
    fassert_equal(s.v.size(), n);
    fassert_equal(s.r.size(), n);
    fassert_equal(s.source.size(), n);
    fassert_equal(s.rho.size(), n);
    fassert_equal(s.termvel.size(), n);
    fassert_equal(s.removed.size(), n);
  }
  static void CheckSize(const ParticlesView& s) {
    const size_t n = s.x.size();
    fassert_equal(s.is_inner.size(), n);
    fassert_equal(s.v.size(), n);
    fassert_equal(s.r.size(), n);
    fassert_equal(s.source.size(), n);
    fassert_equal(s.rho.size(), n);
    fassert_equal(s.termvel.size(), n);
    fassert_equal(s.removed.size(), n);
  }
  template <class F>
  static void ForEachAttribute(const ParticlesView& view, F func) {
    func(view.x);
    func(view.is_inner);
    func(view.v);
    func(view.r);
    func(view.source);
    func(view.rho);
    func(view.termvel);
    func(view.removed);
  }
  static void SwapParticles(const ParticlesView& view, size_t i, size_t j) {
    ForEachAttribute(view, [&](auto& v) { //
      std::swap(v[i], v[j]);
    });
  }
  static void ResizeParticles(const ParticlesView& view, size_t size) {
    ForEachAttribute(view, [&](auto& v) { //
      v.resize(size);
    });
  }
  // Clears particles marked as removed.
  static void ClearRemoved(const ParticlesView& view) {
    const size_t n = view.x.size();
    if (!n) {
      return;
    }
    size_t i = 0;
    for (size_t j = n - 1; j >= i;) { // j: index of last non-removed particle
      if (view.removed[i]) {
        while (j > i && view.removed[j]) {
          --j;
        }
        if (i == j) {
          break;
        }
        SwapParticles(view, i, j);
      }
      ++i;
    }
    { // check
      size_t j = 0;
      while (j < n && !view.removed[j]) {
        ++j;
      }
      fassert_equal(j, i);
      while (j < n && view.removed[j]) {
        ++j;
      }
      fassert_equal(j, n);
    }
    ResizeParticles(view, i);
  }
  void Step(Scal dt, const FieldEmbed<Scal>& fev) {
    auto sem = m.GetSem("step");
    auto& s = state_;
    if (sem("local")) {
      // convert flux to normal velocity component
      FieldFaceb<Scal> ff_vel = fev.template Get<FieldFaceb<Scal>>();
      eb.LoopFaces([&](auto cf) { //
        ff_vel[cf] /= eb.GetArea(cf);
      });
      // restore vector velocity field
      const FieldCell<Vect> fc_vel = UEB::AverageGradient(ff_vel, eb);

      for (size_t i = 0; i < s.x.size(); ++i) {
        const auto c = m.GetCellFromPoint(s.x[i]);
        Scal tau = 0;
        switch (conf.mode) {
          case ParticlesMode::stokes:
            tau = (s.rho[i] - conf.mixture_density) * sqr(2 * s.r[i]) /
                  (18 * conf.mixture_viscosity);
            break;
          case ParticlesMode::termvel:
            tau = s.termvel[i] / conf.gravity.norm();
            break;
          case ParticlesMode::tracer:
            tau = 0;
            break;
        }

        // Update velocity from viscous drag, implicit in time.
        //   dv/dt = (u - v) / tau + g
        // particle velocity `v`, liquid velocity `u`,
        // relaxation time `tau`, gravity `g`
        s.v[i] = (fc_vel[c] + s.v[i] * (tau / dt) + conf.gravity * tau) /
                 (1 + tau / dt);
        s.x[i] += s.v[i] * dt;

        if (s.source[i] != 0) {
          // Apply volume source
          const Scal pi = M_PI;
          const Scal k = 4. / 3 * pi;
          Scal vol = k * std::pow(s.r[i], 3);
          vol = std::max<Scal>(0, vol + s.source[i] * dt);
          s.r[i] = std::pow(vol / k, 1. / 3);
        }

        const auto gl = m.GetGlobalLength();
        for (size_t d = 0; d < m.GetEdim(); ++d) {
          if (m.flags.is_periodic[d]) {
            if (s.x[i][d] < 0) {
              s.x[i][d] += gl[d];
            }
            if (s.x[i][d] > gl[d]) {
              s.x[i][d] -= gl[d];
            }
          }
        }
        // Remove particles that have left the domain
        if (!m.GetGlobalBoundingBox().IsInside(s.x[i]) || eb.IsCut(c) ||
            eb.IsExcluded(c)) {
          s.removed[i] = 1;
        }
      }
      ClearRemoved(GetView(s));
      typename M::CommPartRequest req;
      req.x = &s.x;
      req.is_inner = &s.is_inner;
      req.attr_scal = {&s.r, &s.source, &s.rho, &s.termvel, &s.removed};
      req.attr_vect = {&s.v};
      m.CommPart(req);
    }
    if (sem("stat")) {
      time_ += dt;
    }
  }
  static void Append(State& s, ParticlesView& other) {
    auto append = [](auto& v, const auto& a) {
      v.insert(v.end(), a.begin(), a.end());
    };
    CheckSize(other);
    append(s.x, other.x);
    append(s.is_inner, other.is_inner);
    append(s.v, other.v);
    append(s.r, other.r);
    append(s.source, other.source);
    append(s.rho, other.rho);
    append(s.termvel, other.termvel);
    append(s.removed, other.removed);
    CheckSize(s);
  }
  void DumpCsv(std::string path) const {
    auto sem = m.GetSem("dumpcsv");
    struct {
      State s;
      std::vector<int> block;
    } * ctx(sem);
    auto& t = *ctx;
    if (sem("gather")) {
      CheckSize(state_);
      auto& s = state_;
      t.s = s;
      for (size_t i = 0; i < s.x.size(); ++i) {
        t.block.push_back(m.GetId());
      }
      m.Reduce(&t.s.x, Reduction::concat);
      m.Reduce(&t.s.v, Reduction::concat);
      m.Reduce(&t.s.r, Reduction::concat);
      m.Reduce(&t.s.source, Reduction::concat);
      m.Reduce(&t.s.rho, Reduction::concat);
      m.Reduce(&t.s.termvel, Reduction::concat);
      m.Reduce(&t.s.removed, Reduction::concat);
      m.Reduce(&t.block, Reduction::concat);
    }
    if (sem("write") && m.IsRoot()) {
      std::ofstream o(path);
      o.precision(16);
      // header
      o << "x,y,z,vx,vy,vz,r,source,rho,termvel,removed,block";
      o << std::endl;
      // content
      auto& s = t.s;
      for (size_t i = 0; i < s.x.size(); ++i) {
        o << s.x[i][0] << ',' << s.x[i][1] << ',' << s.x[i][2];
        o << ',' << s.v[i][0] << ',' << s.v[i][1] << ',' << s.v[i][2];
        o << ',' << s.r[i];
        o << ',' << s.source[i];
        o << ',' << s.rho[i];
        o << ',' << s.termvel[i];
        o << ',' << s.removed[i];
        o << ',' << t.block[i];
        o << "\n";
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
  size_t nrecv_;
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
  return imp->GetView(imp->state_);
}

template <class EB_>
void Particles<EB_>::Append(ParticlesView& app) {
  imp->Append(imp->state_, app);
}

template <class EB_>
void Particles<EB_>::DumpCsv(std::string path) const {
  imp->DumpCsv(path);
}

template <class EB_>
auto Particles<EB_>::GetNumRecv() const -> size_t {
  return imp->nrecv_;
}

template <class EB_>
auto Particles<EB_>::GetTime() const -> Scal {
  return imp->time_;
}
