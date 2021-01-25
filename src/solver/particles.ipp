// Created by Petr Karnakov on 26.06.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <array>
#include <exception>
#include <fstream>
#include <limits>
#include <memory>
#include <type_traits>

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
    std::vector<Scal> termvel;
  };

  Imp(Owner* owner, M& m_, const EB& eb_, const ParticlesView& init, Scal time,
      Conf conf_)
      : owner_(owner)
      , m(m_)
      , eb(eb_)
      , conf(conf_)
      , time_(time)
      , state_({init.x, init.v, init.r, init.rho, init.termvel}) {}
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
        const Scal tau_stokes = (s.rho[i] - conf.mixture_density) *
                                sqr(2 * s.r[i]) / (18 * conf.mixture_viscosity);
        const Scal tau =
            (conf.use_termvel ? s.termvel[i] / conf.gravity.norm()
                              : tau_stokes);
        const auto gl = m.GetGlobalLength();
        // implicit drag
        // dv/dt = (u - v) / tau + g
        if (!m.GetGlobalBoundingBox().IsInside(s.x[i]) || eb.IsCut(c) ||
            eb.IsExcluded(c)) {
          s.v[i] = Vect(0);
        } else {
          s.v[i] = //
              (fc_vel[c] + s.v[i] * (tau / dt) + conf.gravity * tau) /
              (1 + tau / dt);
        }
        const Vect x_old = s.x[i];
        s.x[i] += s.v[i] * dt;
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
        if (!m.GetGlobalBoundingBox().IsInside(s.x[i]) || eb.IsCut(c) ||
            eb.IsExcluded(c)) {
          s.x[i] = x_old;
        }
      }
    }
    if (sem.Nested()) {
      Comm(s.x, {&s.r, &s.rho, &s.termvel}, {&s.v}, m, nrecv_);
    }
    if (sem("stat")) {
      time_ += dt;
    }
  }
  static void CheckSize(const State& s) {
    fassert_equal(s.x.size(), s.x.size());
    fassert_equal(s.r.size(), s.x.size());
    fassert_equal(s.rho.size(), s.x.size());
    fassert_equal(s.termvel.size(), s.x.size());
  }
  static void Append(State& s, ParticlesView& app) {
    auto append = [](auto& v, const auto& a) {
      v.insert(v.end(), a.begin(), a.end());
    };
    append(s.x, app.x);
    append(s.v, app.v);
    append(s.r, app.r);
    append(s.rho, app.rho);
    append(s.termvel, app.termvel);
    CheckSize(s);
  }
  // Exchanges particle positions and data between blocks
  // such that the positions are inside the owning block.
  // x: particle positions
  // attr_scal: scalar attributes
  // attr_vect: vector attributes
  // Output:
  // x,attr_scal,attr_vect: updated
  // ncomm: number of recevied particles
  static void Comm(
      std::vector<Vect>& x, const std::vector<std::vector<Scal>*>& attr_scal,
      const std::vector<std::vector<Vect>*>& attr_vect, M& m, size_t& ncomm) {
    auto sem = m.GetSem("particles-comm");
    struct {
      std::vector<Vect> gather_x;
      std::vector<std::vector<Scal>> gather_scal;
      std::vector<std::vector<Vect>> gather_vect;
      std::vector<std::vector<Scal>> scatter_serial;
      std::vector<Scal> recv_serial;
      std::vector<int> id;
    } * ctx(sem);
    if (sem("local")) {
      const auto box = m.GetBoundingBox();
      std::vector<size_t> inside; // indices of particles inside box
      std::vector<size_t> outside; // indices of particles outside box
      for (size_t i = 0; i < x.size(); ++i) {
        if (box.IsInside(x[i])) {
          inside.push_back(i);
        } else {
          outside.push_back(i);
        }
      }

      // Returns elements with indices from `outside`,
      // and overwrites `attr` with elements from `inside`.
      auto select = [](const auto& attr, const std::vector<size_t>& indices) {
        using T = typename std::decay<decltype(attr[0])>::type;
        std::vector<T> res;
        for (auto i : indices) {
          res.push_back(attr[i]);
        }
        return res;
      };

      ctx->gather_x = select(x, outside);
      x = select(x, inside);
      for (auto& attr : attr_scal) {
        ctx->gather_scal.push_back(select(*attr, outside));
        (*attr) = select(*attr, inside);
      }
      for (auto& attr : attr_vect) {
        ctx->gather_vect.push_back(select(*attr, outside));
        (*attr) = select(*attr, inside);
      }
      // gather
      m.Reduce(&ctx->gather_x, Reduction::concat);
      for (auto& attr : ctx->gather_scal) {
        m.Reduce(&attr, Reduction::concat);
      }
      for (auto& attr : ctx->gather_vect) {
        m.Reduce(&attr, Reduction::concat);
      }
      ctx->id.push_back(m.GetId());
      m.Reduce(&ctx->id, Reduction::concat);
    }
    if (sem()) {
      if (m.IsRoot()) {
        const size_t gather_size = ctx->gather_x.size();
        const size_t maxid = ctx->id.size();

        // block id to index in `ctx->id`
        std::vector<size_t> id_to_index(ctx->id.size());
        for (size_t i = 0; i < ctx->id.size(); ++i) {
          id_to_index[ctx->id[i]] = i;
        }
        // index in `gather_x` to index in `ctx->id`.
        std::vector<size_t> remote_index;
        for (size_t k = 0; k < ctx->gather_x.size(); ++k) {
          const size_t id = m.GetIdFromPoint(ctx->gather_x[k]);
          fassert(id < maxid);
          remote_index.push_back(id_to_index[id]);
        }

        // FIXME: serialization, revise with m.Scatter() that takes Vect
        auto scatter = [&remote_index, gather_size, maxid](
                           std::vector<std::vector<Scal>>& vscatter,
                           const std::vector<Scal>& vgather) {
          fassert_equal(vgather.size(), gather_size);
          vscatter.resize(maxid);
          for (size_t i = 0; i < vgather.size(); ++i) {
            vscatter[remote_index[i]].push_back(vgather[i]);
          }
        };
        auto scatterv = [&remote_index, gather_size, maxid](
                            std::vector<std::vector<Scal>>& vscatter,
                            const std::vector<Vect>& vgather) {
          fassert_equal(vgather.size(), gather_size);
          vscatter.resize(maxid);
          for (size_t i = 0; i < vgather.size(); ++i) {
            vscatter[remote_index[i]].push_back(vgather[i][0]);
            vscatter[remote_index[i]].push_back(vgather[i][1]);
            vscatter[remote_index[i]].push_back(vgather[i][2]);
          }
        };

        scatterv(ctx->scatter_serial, ctx->gather_x);
        for (size_t i = 0; i < attr_scal.size(); ++i) {
          scatter(ctx->scatter_serial, ctx->gather_scal[i]);
        }
        for (size_t i = 0; i < attr_vect.size(); ++i) {
          scatterv(ctx->scatter_serial, ctx->gather_vect[i]);
        }
      }
      m.Scatter({&ctx->scatter_serial, &ctx->recv_serial});
    }
    if (sem()) {
      const size_t nscal = 3 + attr_scal.size() + attr_vect.size() * 3;
      fassert_equal(ctx->recv_serial.size() % nscal, 0);
      // number of particles received
      const size_t recv_size = ctx->recv_serial.size() / nscal;
      ncomm = recv_size;
      auto deserial = [recv_size](auto& attr, auto& serial) {
        for (size_t i = 0; i < recv_size; ++i) {
          attr.push_back(serial.back());
          serial.pop_back();
        }
      };
      auto deserialv = [recv_size](auto& attr, auto& serial) {
        for (size_t i = 0; i < recv_size; ++i) {
          Vect v;
          v[2] = serial.back();
          serial.pop_back();
          v[1] = serial.back();
          serial.pop_back();
          v[0] = serial.back();
          serial.pop_back();
          attr.push_back(v);
        }
      };
      // ctx->recv_serial contains serialized particles from root
      // deserialize in reversed order
      for (size_t i = attr_vect.size(); i > 0;) {
        --i;
        deserialv(*attr_vect[i], ctx->recv_serial);
      }
      for (size_t i = attr_scal.size(); i > 0;) {
        --i;
        deserial(*attr_scal[i], ctx->recv_serial);
      }
      deserialv(x, ctx->recv_serial);
    }
    if (sem()) { // FIXME: empty stage
    }
  }
  void DumpCsv(std::string path) const {
    auto sem = m.GetSem("dumpcsv");
    struct {
      State s;
      std::vector<int> block;
    } * ctx(sem);
    if (sem("gather")) {
      CheckSize(state_);
      auto& s = state_;
      ctx->s.x = s.x;
      ctx->s.v = s.v;
      ctx->s.r = s.r;
      ctx->s.rho = s.rho;
      ctx->s.termvel = s.termvel;
      for (size_t i = 0; i < s.x.size(); ++i) {
        ctx->block.push_back(m.GetId());
      }
      m.Reduce(&ctx->s.x, Reduction::concat);
      m.Reduce(&ctx->s.v, Reduction::concat);
      m.Reduce(&ctx->s.r, Reduction::concat);
      m.Reduce(&ctx->s.rho, Reduction::concat);
      m.Reduce(&ctx->s.termvel, Reduction::concat);
      m.Reduce(&ctx->block, Reduction::concat);
    }
    if (sem("write") && m.IsRoot()) {
      std::ofstream o(path);
      o.precision(16);
      // header
      o << "x,y,z,vx,vy,vz,r,rho,termvel,block";
      o << std::endl;
      // content
      auto& s = ctx->s;
      for (size_t i = 0; i < s.x.size(); ++i) {
        o << s.x[i][0] << ',' << s.x[i][1] << ',' << s.x[i][2];
        o << ',' << s.v[i][0] << ',' << s.v[i][1] << ',' << s.v[i][2];
        o << ',' << s.r[i];
        o << ',' << s.rho[i];
        o << ',' << s.termvel[i];
        o << ',' << ctx->block[i];
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
  return {
      imp->state_.x, imp->state_.v, imp->state_.r, imp->state_.rho,
      imp->state_.termvel};
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
