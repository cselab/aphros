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
#include <unordered_map>

#include "approx.h"
#include "approx2.h"
#include "approx_eb.h"
#include "util/format.h"

#include "dump/dump.h"
#include "dump/raw.h"

#include "particles.h"

template <class EB_>
struct Particles<EB_>::Imp {
  using Owner = Particles<EB_>;
  using UEB = UEmbed<M>;

  struct State {
    // See description of attributes in ParticlesView.
    std::vector<Vect> x;
    std::vector<bool> inner;
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
      , state_(
            {init.x, init.inner, init.v, init.r, init.source, init.rho,
             init.termvel, init.removed}) {
    CheckSize(init);
    auto& s = state_;
    const auto req = GetCommPartRequest<M>(GetView(s));
    m.CommPart(req);
  }
  static ParticlesView GetView(State& s) {
    return {s.x, s.inner, s.v, s.r, s.source, s.rho, s.termvel, s.removed};
  }
  static void CheckSize(const State& s) {
    const size_t n = s.x.size();
    fassert_equal(s.inner.size(), n);
    fassert_equal(s.v.size(), n);
    fassert_equal(s.r.size(), n);
    fassert_equal(s.source.size(), n);
    fassert_equal(s.rho.size(), n);
    fassert_equal(s.termvel.size(), n);
    fassert_equal(s.removed.size(), n);
  }
  static void CheckSize(const ParticlesView& s) {
    const size_t n = s.x.size();
    fassert_equal(s.inner.size(), n);
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
    func(view.inner);
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
  void Step(
      Scal dt, const FieldEmbed<Scal>& fev,
      const MapEmbed<BCond<Vect>>& mebc_velocity,
      std::function<void(const ParticlesView&)> velocity_hook) {
    auto sem = m.GetSem("step");
    struct {
      FieldFace<Scal> ff_veln;
    } * ctx(sem);
    auto& t = *ctx;
    auto& s = state_;
    if (sem()) {
      t.ff_veln = fev.template Get<FieldFaceb<Scal>>();
      for (auto f : m.Faces()) {
        t.ff_veln[f] /= eb.GetArea(f);
      }
    }
    if (sem.Nested()) {
      Approx2<EB>::ExtrapolateToHaloFaces(t.ff_veln, mebc_velocity, m);
    }
    if (sem("local")) {
      /*
      if (m.IsRoot() && time_ == 0) {
        using MIdx = typename M::MIdx;
        const MIdx size(64);
        GBlock<IdxCell, dim> block(size);
        GIndex<IdxCell, dim> index(size);
        FieldCell<Scal> fc(index);
        const IdxCell c0 = m.GetIndexCells().GetIdx(MIdx(0));
        auto callback = [&](const std::function<Vect(Vect x)>& func) {
          const Vect h = m.GetCellSize();
          for (MIdx w : block) {
            Vect x((Vect(w) + Vect(0.5)) / Vect(size));
            x = m.GetCenter(c0) - h * 0.5 + x * (2 * h);
            fc[index.GetIdx(w)] = func(x)[0];
          }
        };
        Approx2<EB>::EvalTrilinearFromFaceField(t.ff_veln, callback, m);
        dump::Raw<M>::WritePlainArrayWithXmf("test.raw", "u", fc.data(), size);
      }
      */

      // Project liquid velocity to particles.
      std::vector<Vect> v_liquid(s.x.size());
      auto callback = [&](const std::function<Vect(Vect x)>& func) {
        for (size_t i = 0; i < s.x.size(); ++i) {
          v_liquid[i] = func(s.x[i]);
        }
      };
      Approx2<EB>::EvalTrilinearFromFaceField(t.ff_veln, callback, m);

      // Compute velocity on particles and advance positions.
      for (size_t i = 0; i < s.x.size(); ++i) {
        const auto c = m.GetCellFromPoint(s.x[i]);
        // Update velocity.
        {
          Vect& v = s.v[i];
          const Vect g = conf.gravity;
          const Vect u = v_liquid[i];
          switch (conf.mode) {
            case ParticlesMode::stokes: {
              // Acceleration from Stokes drag and gravity, implicit in time
              //   rho_p V dv/dt = 6 pi mu R (u - v) + (rho_p - rho) g V
              // with particle volume `V`, particle density `rho_p`,
              // particle velocity `v`, liquid velocity `u`,
              // fluid viscosity `mu` and density `rho`,
              // particle volume `V`, gravity `g`
              // With notation for coefficients,
              //   kt dv/dt = km (u - v) + kg g
              const Scal r = s.r[i];
              const Scal rho_p = s.rho[i];
              const Scal rho = conf.mixture_density;
              const Scal mu = conf.mixture_viscosity;
              // XXX Generated in `gen/particles_stokes.py`.
              v = (9 * dt * mu * u +
                   g * (-2 * dt * sqr(r) * rho + 2 * dt * sqr(r) * rho_p) +
                   2 * sqr(r) * rho_p * v) /
                  (9 * dt * mu + 2 * sqr(r) * rho_p);
              break;
            }
            case ParticlesMode::termvel:
              // Given terminal velocity.
              v = s.termvel[i] * g.normalized();
              break;
            case ParticlesMode::tracer:
              v = u;
              break;
          }
        }

        // Update velocity from viscous drag, implicit in time.
        //   dv/dt = (u - v) / tau + g
        // particle velocity `v`, liquid velocity `u`,
        // relaxation time `tau`, gravity `g`
        /*
        s.v[i] = (v_liquid[i] + s.v[i] * (tau / dt) + conf.gravity * tau) /
                 (1 + tau / dt);
                 */
        // s.v[i] += dt * (s.v[i] / tau + conf.gravity);
        s.v[i] = 2. / 9 * (s.rho[i] - conf.mixture_density) /
                 conf.mixture_viscosity * conf.gravity * sqr(s.r[i]);
        velocity_hook(GetView(s));
        s.x[i] += s.v[i] * dt;

        if (s.source[i] != 0) {
          // Apply volume source.
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
        // Remove particles that have left the domain.
        if (!m.GetGlobalBoundingBox().IsInside(s.x[i]) || eb.IsCut(c) ||
            eb.IsExcluded(c)) {
          s.removed[i] = 1;
        }
      }
      const auto& view = GetView(s);
      ClearRemoved(view);
      const auto req = GetCommPartRequest<M>(view);
      m.CommPart(req);
    }
    if (sem("stat")) {
      time_ += dt;
    }
  }
  static void Append(State& s, const ParticlesView& other) {
    auto append = [](auto& v, const auto& a) {
      v.insert(v.end(), a.begin(), a.end());
    };
    CheckSize(other);
    append(s.x, other.x);
    append(s.inner, other.inner);
    append(s.v, other.v);
    append(s.r, other.r);
    append(s.source, other.source);
    append(s.rho, other.rho);
    append(s.termvel, other.termvel);
    append(s.removed, other.removed);
    CheckSize(s);
  }
  void DumpCsv(
      const std::string& path,
      const std::unordered_set<std::string>& sel) const {
    auto sem = m.GetSem("dumpcsv");
    struct {
      State s;
      std::vector<int> block;
      std::vector<Scal> inner;
    } * ctx(sem);
    auto& t = *ctx;
    if (sem("gather")) {
      CheckSize(state_);
      auto& s = state_;
      t.s = s;
      for (size_t i = 0; i < s.x.size(); ++i) {
        t.block.push_back(m.GetId());
        t.inner.push_back(s.inner[i]);
      }
      m.Reduce(&t.s.x, Reduction::concat);
      m.Reduce(&t.s.v, Reduction::concat);
      m.Reduce(&t.s.r, Reduction::concat);
      m.Reduce(&t.s.source, Reduction::concat);
      m.Reduce(&t.s.rho, Reduction::concat);
      m.Reduce(&t.s.termvel, Reduction::concat);
      m.Reduce(&t.s.removed, Reduction::concat);
      m.Reduce(&t.block, Reduction::concat);
      m.Reduce(&t.inner, Reduction::concat);
    }
    if (sem("write") && m.IsRoot()) {
      const auto& s = t.s;
      std::vector<std::pair<std::string, std::vector<Scal>>> data;
      // Dump halo particles if field "inner" is selected.
      const bool dump_halo = sel.count("inner");
      auto append_scal = [&](std::string name, const auto& v) {
        if (sel.count(name)) {
          std::vector<Scal> res;
          for (size_t i = 0; i < v.size(); ++i) {
            if (t.inner[i] || dump_halo) {
              res.push_back(v[i]);
            }
          }
          data.emplace_back(name, std::move(res));
        }
      };
      auto append_vect = [&](std::string prefix, const std::vector<Vect>& v) {
        for (auto d : m.dirs) {
          const std::string dletter(1, GDir<dim>(d).letter());
          const auto name = prefix + dletter;
          if (sel.count(name)) {
            std::vector<Scal> res;
            for (size_t i = 0; i < v.size(); ++i) {
              if (t.inner[i] || dump_halo) {
                res.push_back(v[i][d]);
              }
            }
            data.emplace_back(name, std::move(res));
          }
        }
      };
      append_vect("", s.x); // Position.
      append_vect("v", s.v); // Velocity.
      append_scal("r", s.r);
      append_scal("source", s.source);
      append_scal("rho", s.rho);
      append_scal("termvel", s.termvel);
      append_scal("removed", s.removed);
      append_scal("block", t.block);
      append_scal("inner", t.inner);
      dump::DumpCsv(data, path);
    }
    if (sem()) {
    }
  }
  static ReadCsvStatus ReadCsv(
      std::istream& in, const ParticlesView& view, char delim) {
    ReadCsvStatus status;
    const auto data = dump::ReadCsv<Scal>(in, delim);
    auto set_component = [](std::vector<Vect>& dst, size_t d,
                            const std::vector<Scal>& src) {
      fassert(d < Vect::dim);
      dst.resize(src.size());
      for (size_t i = 0; i < src.size(); ++i) {
        dst[i][d] = src[i];
      }
    };
    for (const auto& p : data) {
      const auto name = p.first;
      for (size_t d = 0; d < Vect::dim; ++d) {
        const std::string dletter(1, GDir<dim>(d).letter());
        // Position.
        if (name == dletter) {
          set_component(view.x, d, p.second);
          status.x = true;
        }
        // Velocity.
        if (name == "v" + dletter) {
          set_component(view.v, d, p.second);
          status.v = true;
        }
      }
      if (name == "r") {
        view.r = p.second;
        status.r = true;
      }
      if (name == "source") {
        view.source = p.second;
        status.source = true;
      }
      if (name == "rho") {
        view.rho = p.second;
        status.rho = true;
      }
      if (name == "termvel") {
        view.termvel = p.second;
        status.termvel = true;
      }
    }
    // Fill the remaining fields with default values.
    const size_t n = view.x.size();
    view.inner.resize(n, true);
    view.v.resize(n, Vect(0));
    view.r.resize(n, 0);
    view.source.resize(n, 0);
    view.rho.resize(n, 0);
    view.termvel.resize(n, 0);
    view.removed.resize(n, false);
    return status;
  }
  static ReadCsvStatus ReadCsv(
      const std::string& path, const ParticlesView& view, char delim) {
    std::ifstream fin(path);
    fassert(fin.good(), "Can't open file '" + path + "' for reading");
    return ReadCsv(fin, view, delim);
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
void Particles<EB_>::Step(
    Scal dt, const FieldEmbed<Scal>& fe_flux,
    const MapEmbed<BCond<Vect>>& mebc_velocity,
    std::function<void(const ParticlesView&)> velocity_hook) {
  imp->Step(dt, fe_flux, mebc_velocity, velocity_hook);
}

template <class EB_>
auto Particles<EB_>::GetView() const -> ParticlesView {
  return imp->GetView(imp->state_);
}

template <class EB_>
void Particles<EB_>::Append(const ParticlesView& app) {
  imp->Append(imp->state_, app);
}

template <class EB_>
void Particles<EB_>::DumpCsv(
    const std::string& path, const std::unordered_set<std::string>& sel) const {
  imp->DumpCsv(path, sel);
}

template <class EB_>
auto Particles<EB_>::ReadCsv(
    std::istream& fin, const ParticlesView& view, char delim) -> ReadCsvStatus {
  return Imp::ReadCsv(fin, view, delim);
}

template <class EB_>
auto Particles<EB_>::ReadCsv(
    const std::string& path, const ParticlesView& view, char delim)
    -> ReadCsvStatus {
  return Imp::ReadCsv(path, view, delim);
}

template <class EB_>
auto Particles<EB_>::GetNumRecv() const -> size_t {
  return imp->nrecv_;
}

template <class EB_>
auto Particles<EB_>::GetTime() const -> Scal {
  return imp->time_;
}
