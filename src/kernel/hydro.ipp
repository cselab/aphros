// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <array>
#include <cassert>
#include <chrono>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>

#include "debug/isnan.h"
#include "dump/dump.h"
#include "dump/dumper.h"
#include "dump/hdf.h"
#include "dump/raw.h"
#include "func/init.h"
#include "func/init_contang.h"
#include "func/init_u.h"
#include "geom/mesh.h"
#include "kernelmeshpar.h"
#include "linear/linear.h"
#include "parse/curv.h"
#include "parse/parser.h"
#include "parse/proj.h"
#include "parse/simple.h"
#include "parse/util.h"
#include "parse/vars.h"
#include "parse/vof.h"
#include "solver/advection.h"
#include "solver/approx.h"
#include "solver/approx2.h"
#include "solver/approx_eb.h"
#include "solver/curv.h"
#include "solver/electro.h"
#include "solver/embed.h"
#include "solver/fluid_dummy.h"
#include "solver/multi.h"
#include "solver/normal.h"
#include "solver/particles.h"
#include "solver/proj.h"
#include "solver/reconst.h"
#include "solver/simple.h"
#include "solver/solver.h"
#include "solver/tracer.h"
#include "solver/vof.h"
#include "solver/vofm.h"
#include "util/convdiff.h"
#include "util/events.h"
#include "util/filesystem.h"
#include "util/format.h"
#include "util/hydro.h"
#include "util/hydro_post.h"
#include "util/linear.h"
#include "util/metrics.h"
#include "util/posthook.h"
#include "util/stat.h"
#include "util/sysinfo.h"
#include "util/timer.h"
#include "util/vof.h"

#include "hydro.h"

template <class Vect>
Vect Sqrt(Vect v) {
  for (size_t d = 0; d < Vect::dim; ++d) {
    v[d] = std::sqrt(v[d]);
  }
  return v;
}

template <class M>
void Hydro<M>::UpdateAdvectionPar() {
  VofPar<M> par = ParseVofPar<M>(var);
  if (auto as = dynamic_cast<ASV*>(as_.get())) {
    par.labeling = as->GetPar().labeling;
    as->SetPar(par);
  }
  if (auto as = dynamic_cast<ASVM*>(as_.get())) {
    par.labeling = as->GetPar().labeling;
    as->SetPar(par);
  }
  if (auto as = dynamic_cast<ASVEB*>(as_.get())) {
    par.labeling = as->GetPar().labeling;
    as->SetPar(par);
  }
  if (auto as = dynamic_cast<ASVMEB*>(as_.get())) {
    par.labeling = as->GetPar().labeling;
    as->SetPar(par);
  }
}
template <class M>
typename M::Scal Hydro<M>::GetSurfaceTensionDt() const {
  const Scal sigma = var.Double["sigma"];
  const Scal cflst = var.Double("cflst", 0);
  if (cflst > 0 && sigma != 0) {
    const Scal pi = M_PI;
    const Scal h3 = std::pow(m.GetCellSize()[0], 3);
    const Scal rho1 = var.Double["rho1"];
    const Scal rho2 = var.Double["rho2"];
    return cflst * std::sqrt(h3 * (rho1 + rho2) / (4. * pi * std::abs(sigma)));
  }
  return 0;
}
template <class M>
typename M::Scal Hydro<M>::GetViscosityDt() const {
  const Scal rho1 = var.Double["rho1"];
  const Scal rho2 = var.Double["rho2"];
  const Scal mu1 = var.Double["mu1"];
  const Scal mu2 = var.Double["mu2"];
  const Scal nu1 = mu1 / rho1;
  const Scal nu2 = mu2 / rho2;
  const Scal nu_max = std::max(nu1, nu2);
  const Scal cflvis = var.Double("cflvis", 0);
  if (cflvis > 0 && nu_max != 0) {
    const Scal h2 = sqr(m.GetCellSize()[0]);
    return cflvis * h2 / nu_max;
  }
  return 0;
}
template <class M>
void Hydro<M>::CalcVort() {
  auto& fcv = fs_->GetVelocity();
  if (eb_) {
    fcom_ = UEmbed<M>::GetVort(fcv, fs_->GetVelocityCond(), *eb_);
  } else {
    fcom_ = UEmbed<M>::GetVort(fcv, fs_->GetVelocityCond(), m);
  }
  fcomm_.Reinit(m);
  for (auto c : m.Cells()) {
    fcomm_[c] = fcom_[c].norm();
  }
}
template <class M>
auto Hydro<M>::CalcStrain(const FieldCell<Vect> fcvel) const
    -> FieldCell<Scal> {
  auto& fcv = fcvel;
  auto ffv = UEB::Interpolate(fcv, fs_->GetVelocityCond(), m);

  std::array<FieldCell<Vect>, dim> g; // g[i][c][j] is derivative du_i/dx_j
  for (size_t i = 0; i < dim; ++i) {
    g[i] = UEB::AverageGradient(GetComponent(ffv, i), m);
  }

  FieldCell<Scal> fcs(m, 0);
  int edim = var.Int["dim"];
  for (auto c : m.Cells()) {
    for (int i = 0; i < edim; ++i) {
      for (int j = 0; j < edim; ++j) {
        fcs[c] += sqr(g[i][c][j]) + g[i][c][j] * g[j][c][i];
      }
    }
    fcs[c] *= 0.5;
  }
  return fcs;
}
template <class M>
auto Hydro<M>::GetDiv() const -> FieldCell<Scal> {
  if (eb_) {
    return Approx2<EB>::GetRegularDivergence(fs_->GetVolumeFlux(), *eb_);
  } else {
    return Approx2<M>::GetRegularDivergence(
        fs_->GetVolumeFlux().GetFieldFace(), m);
  }
}

template <class M>
void Hydro<M>::InitEmbed() {
  if (var.Int["enable_embed"]) {
    auto sem = m.GetSem("embed");
    struct {
      FieldNode<Scal> fnl;
    } * ctx(sem);
    if (sem("ctor")) {
      eb_.reset(new Embed<M>(m));
    }
    if (sem.Nested("levelset")) {
      UEB::InitLevelSet(ctx->fnl, m, var, m.IsRoot() && !silent_);
    }
    if (sem("hook")) {
      InitEmbedHook(ctx->fnl, var, m);
    }
    if (sem.Nested("smoothen")) {
      SmoothenNode(ctx->fnl, m, var.Int["embed_smoothen_iters"]);
    }
    if (sem.Nested("init")) {
      eb_->Init(ctx->fnl);
    }
  }
}
template <class M>
void Hydro<M>::InitStepwiseBody(
    FieldCell<bool>& fc_innermask, FieldFace<bool>& ff_innermask) {
  auto sem = m.GetSem("stepwise");
  using Xmf = dump::Xmf<Vect>;
  struct {
    FieldCell<Scal> fc_inner; // 1 inner cells, 0 inside solid body
    Vars varbody;
    typename Xmf::Meta meta;
  } * ctx(sem);
  auto& t = *ctx;
  if (sem()) {
    auto& varbody = t.varbody;
    varbody.String.Set("init_vf", var.String["body_init"]);
    varbody.String.Set("list_path", var.String["body_list_path"]);
    varbody.Int.Set("dim", var.Int["dim"]);
    varbody.Int.Set("list_ls", 3);
  }
  if (sem.Nested("body-mask")) {
    // Fill volume fraction of body in `t.fc_inner`.
    // Halos cells are filled from periodic boundary conditions.
    InitVf(t.fc_inner, t.varbody, m, !silent_);
  }
  if (sem("body-bc")) {
    if (var.Int["body_init_inverse"]) {
      for (auto c : m.AllCells()) {
        t.fc_inner[c] = 1 - t.fc_inner[c];
      }
    }
    // Set volume fraction of body outside the domain to 1.
    const MIdx gs = m.GetGlobalSize();
    for (auto c : m.AllCellsM()) {
      for (auto d : m.dirs) {
        if (!m.flags.is_periodic[d] && (c[d] < 0 || c[d] >= gs[d])) {
          t.fc_inner[c] = 0;
        }
      }
    }

    fc_innermask.Reinit(m);
    const Scal thres = var.Double["stepwise_threshold"];
    for (auto c : m.AllCells()) {
      t.fc_inner[c] = (t.fc_inner[c] > thres ? 1 : 0);
      fc_innermask[c] = t.fc_inner[c];
    }
    // Inner faces are those adjacent to an inner cell.
    ff_innermask.Reinit(m, false);
    for (auto c : m.AllCellsM()) {
      if (fc_innermask[c]) {
        for (auto q : m.Nci(c)) {
          ff_innermask[c.face(q)] = true;
        }
      }
    }
  }
  if (var.Int["dump_stepwise_body"]) {
    if (sem.Nested()) {
      // Dump scalar field which is 1 in inner cells and 0 in excluded cells.
      dump::Raw<M>::WriteWithXmf(t.fc_inner, "stepwise", "stepwise.raw", m);
    }
  }
}
template <class M>
void Hydro<M>::SpawnTracer() {
  const Scal sphere_r = var.Double["tracer_spawn_sphere_r"];
  if (sphere_r > 0) {
    auto& conf = tracer_->GetConf();
    const auto tracer_layers = GRange<size_t>(conf.layers);
    const Vect sphere_c(var.Vect["tracer_spawn_sphere_c"]);
    const size_t edim = m.GetEdim();
    auto vfcu = tracer_->GetVolumeFraction();
    for (auto l : tracer_layers) {
      const std::string prefix = "tracer" + std::to_string(l);
      auto k = var.Double[prefix + "_factor"];
      for (auto c : m.AllCells()) {
        const auto xc = m.GetCenter(c);
        Vect dx = xc - sphere_c;
        if (edim == 2 && M::dim > 2) {
          dx[2] = 0;
        }
        if (dx.sqrnorm() < sqr(sphere_r)) {
          vfcu[l][c] = k;
        }
      }
    }
    tracer_->SetVolumeFraction(vfcu);
  }
}

template <class M>
void Hydro<M>::OverwriteBc() {
  // piecewise-linear function
  auto piecewise = [&](Scal t, const std::vector<Scal>& times,
                       const std::vector<Scal>& values) {
    fassert_equal(values.size(), times.size());
    if (times.size() == 0) {
      return GetNan<Scal>();
    }
    size_t i = 0;
    while (i < times.size() && times[i] <= t) {
      ++i;
    }
    if (i == 0) { // t < times[0]
      return GetNan<Scal>();
    }
    if (i < times.size()) { // times[i - 1] <= t < times[i]
      const Scal t0 = times[i - 1];
      const Scal t1 = times[i];
      const Scal v0 = values[i - 1];
      const Scal v1 = values[i];
      return t0 < t1 ? v0 + (v1 - v0) * (t - t0) / (t1 - t0) : v0;
    } else {
      return values.back();
    }
  };

  // restore from original
  mebc_fluid_ = mebc_fluid_orig_;

  { // correct inlet velocity by factor
    const auto factor = piecewise(
        fs_->GetTime(), var.Vect["overwrite_inlet_times"],
        var.Vect["overwrite_inlet_factors"]);
    if (!IsNan(factor)) {
      mebc_fluid_.LoopPairs([&](auto cf_bc) {
        auto& curr = mebc_fluid_[cf_bc.first];
        if (curr.type == BCondFluidType::inlet) {
          curr.velocity *= factor;
        }
      });
    }
  }

  { // apply body velocity to meshvel and subtract from boundary conditions
    Vect bodyvel;
    bodyvel[0] = piecewise(
        fs_->GetTime(), var.Vect["bodyvel_times"], var.Vect["bodyvel_x"]);
    bodyvel[1] = piecewise(
        fs_->GetTime(), var.Vect["bodyvel_times"], var.Vect["bodyvel_y"]);
    if (M::dim > 2) {
      bodyvel[2] = piecewise(
          fs_->GetTime(), var.Vect["bodyvel_times"], var.Vect["bodyvel_z"]);
    }
    if (!IsNan(bodyvel)) {
      st_.meshvel = bodyvel;
      mebc_fluid_.LoopPairs([&](auto cf_bc) {
        auto& curr = mebc_fluid_[cf_bc.first];
        curr.velocity = bodyvel;
      });
    } else {
      st_.meshvel = Vect(0);
    }
  }

  { // overwrite inlet pressure
    const auto p_new = piecewise(
        fs_->GetTime(), var.Vect["overwrite_inletpressure_times"],
        var.Vect["overwrite_inletpressure_pressure"]);
    if (!IsNan(p_new)) {
      mebc_fluid_.LoopPairs([&](auto cf_bc) {
        auto& curr = mebc_fluid_[cf_bc.first];
        if (curr.type == BCondFluidType::inletpressure) {
          curr.pressure = p_new;
        }
      });
    }
  }

  if (var.Int["enable_inlet_periodic"]) {
    const Scal t0 = var.Double["inlet_periodic_t0"];
    const Scal dur0 = var.Double["inlet_periodic_dur0"];
    const Scal dur1 = var.Double["inlet_periodic_dur1"];
    const Scal vf0 = var.Double["inlet_periodic_vf0"];
    const Scal vf1 = var.Double["inlet_periodic_vf1"];
    const Scal t = fs_->GetTime();
    if (t >= t0) {
      Scal intpart;
      const Scal fractpart = std::modf((t - t0) / (dur0 + dur1), &intpart);
      const Scal vf = (fractpart <= dur0 / (dur0 + dur1) ? vf0 : vf1);
      mebc_fluid_.LoopPairs([&](auto cf_bc) {
        auto cf = cf_bc.first;
        auto& curr = mebc_fluid_[cf];
        if (curr.type == BCondFluidType::inlet ||
            curr.type == BCondFluidType::inletpressure) {
          mebc_adv_[cf].fill_vf = vf;
        }
      });
    }
  }
}

template <class M>
void Hydro<M>::SpawnParticles(ParticlesView& view) {
  const Vect sphere_c(var.Vect["particles_spawn_sphere_c"]);
  Vect sphere_r;
  if (auto* v = var.Vect.Find("particles_spawn_sphere_r")) {
    sphere_r = Vect(*v);
  } else {
    sphere_r = Vect(var.Double["particles_spawn_sphere_r"]);
  }
  const size_t edim = m.GetEdim();
  const Scal sphere_vol =
      (edim == 3 ? 4. / 3. * M_PI * sphere_r.prod()
                 : M_PI * sphere_r[0] * sphere_r[1]);
  const Vect h = m.GetCellSize();
  const Scal cell_vol = (edim == 3 ? h.prod() : h[0] * h[1]);
  const Vect velocity(var.Vect["particles_spawn_velocity"]);
  const Scal density = var.Double["particles_density"];
  const auto spawn_rate = var.Vect["particles_spawn_rate"];
  const auto radius = var.Vect["particles_radius"];
  const auto termvel = var.Vect["particles_termvel"];
  const size_t num_rates = spawn_rate.size();
  fassert_equal(radius.size(), num_rates);
  fassert_equal(termvel.size(), num_rates);
  std::uniform_real_distribution<Scal> u(0, 1);
  std::uniform_real_distribution<Scal> um(-0.5, 0.5);
  auto& g = randgen_;

  for (auto c : m.Cells()) {
    const auto xc = m.GetCenter(c);
    Vect dx = xc - sphere_c;
    if (edim == 2 && M::dim > 2) {
      dx[2] = 0;
    }
    if ((dx / sphere_r).sqrnorm() < 1) {
      for (size_t i = 0; i < num_rates; ++i) {
        const Scal prob = particles_dt_ * cell_vol * spawn_rate[i] / sphere_vol;
        if (u(g) < prob) {
          Vect xrand;
          for (auto d : M::dirs) {
            xrand[d] = um(g);
          }
          view.x.push_back(m.GetCenter(c) + xrand * h);
          view.inner.push_back(true);
          view.v.push_back(velocity);
          view.id.push_back(0);
          view.r.push_back(radius[i]);
          view.source.push_back(0);
          view.rho.push_back(density);
          view.termvel.push_back(termvel[i]);
          view.removed.push_back(0);
        }
      }
    }
  }
}

template <class M>
void Hydro<M>::InitParticles() {
  typename ParticlesInterface<M>::Conf conf;
  conf.mixture_density = var.Double["rho1"];
  conf.mixture_viscosity = var.Double["mu1"];
  conf.gravity = Vect(var.Vect["gravity"]);
  const auto mode = var.String["particles_mode"];
  if (mode == "tracer") {
    conf.mode = ParticlesMode::tracer;
  } else if (mode == "stokes") {
    conf.mode = ParticlesMode::stokes;
  } else if (mode == "termvel") {
    conf.mode = ParticlesMode::termvel;
  } else {
    fassert(
        false,
        "Unknown mode=" + mode + ". Known modes are tracer, stokes, termvel");
  }
  std::vector<Vect> p_x;
  std::vector<bool> p_inner;
  std::vector<Vect> p_v;
  std::vector<Scal> p_id;
  std::vector<Scal> p_r;
  std::vector<Scal> p_source;
  std::vector<Scal> p_rho;
  std::vector<Scal> p_termvel;
  std::vector<Scal> p_removed;
  ParticlesView view{p_x,      p_inner, p_v,       p_id,     p_r,
                     p_source, p_rho,   p_termvel, p_removed};
  const auto init_csv = var.String["particles_init_csv"];
  if (init_csv.length() && m.IsRoot()) {
    auto status = Particles<M>::ReadCsv(init_csv, view);
    if (!status.termvel) {
      std::fill(
          view.termvel.begin(), view.termvel.end(),
          var.Vect["particles_termvel"][0]);
    }
    if (!status.r) {
      std::fill(view.r.begin(), view.r.end(), var.Vect["particles_radius"][0]);
    }
  }
  SpawnParticles(view);
  if (eb_) {
    particles_.reset(new Particles<EB>(m, *eb_, view, fs_->GetTime(), conf));
  } else {
    particles_.reset(new Particles<M>(m, m, view, fs_->GetTime(), conf));
  }

  fc_wall_dist_.Reinit(m, GetNan<Vect>());
  for (const auto& p : mebc_fluid_.GetMapFace()) {
    const IdxFace f = p.first;
    const auto& bc = p.second;
    const auto nci = bc.nci;
    auto cc = m.GetCellColumn(f, nci);
    if (bc.type == BCondFluidType::wall ||
        bc.type == BCondFluidType::slipwall ||
        bc.type == BCondFluidType::symm) {
      for (IdxCell c : {cc[2], cc[3]}) {
        const Vect dist = m.GetCenter(f) - m.GetCenter(c);
        auto& mindist = fc_wall_dist_[c];
        if (IsNan(mindist) || dist.sqrnorm() < mindist.sqrnorm()) {
          mindist = dist;
        }
      }
    }
    // TODO: consider mebc_fluid_.GetMapCell()
  }
}

// Nucleation points are randomly positioned on walls
// with a given areal number density.
template <class M>
void Hydro<M>::InitNucleationPoints() {
    // Areal number density of nucleation points.
    const Scal numdens = var.Double("nucleation_number_density", 0);
    if (numdens) {
      std::uniform_real_distribution<Scal> uniform(0, 1);
      // Traverse boundary faces and seed at most one point per face.
      for (const auto& p : mebc_fluid_.GetMapFace()) {
        const auto f = m(p.first);
        const auto& bc = p.second;
        if (m.IsInner(f.cell(bc.nci))) {
          if (const size_t* group = me_group_.find(f)) {
            const auto& custom = bc_group_custom_[*group];
            const auto it = custom.find("nucleation");
            if (it != custom.end()) {
              const Scal factor = it->second;
              if (numdens * factor * f.area > uniform(randgen_)) {
                // Uniform displacement.
                Vect dx(0);
                for (size_t d : m.dirs) {
                  dx[d] = uniform(randgen_) - 0.5;
                }
                dx *= m.GetCellSize();
                // Normal towards inner cell.
                const Vect normal = (bc.nci == 0 ? -1 : 1) * f.normal();
                // Cancel out the component normal to the face.
                dx = dx.orth(normal);
                // Move inside the inner cell.
                dx += normal * (m.GetCellSize()[0] * 1e-3);
                nucl_points_.insert(f.center() + dx);
              }
            }
          }
        }
      }
    }
}

template <class M>
void Hydro<M>::InitElectro() {
  if (var.Int["enable_electro"]) {
    if (!linsolver_electro_) {
      if (var.String.Find("linsolver_electro")) {
        linsolver_electro_ = ULinear<M>::MakeLinearSolver(var, "electro", m);
      } else if (linsolver_symm_) {
        linsolver_electro_ = linsolver_symm_;
      } else {
        linsolver_electro_ = ULinear<M>::MakeLinearSolver(var, "symm", m);
      }
    }
    typename ElectroInterface<M>::Conf conf{var, linsolver_electro_};
    if (eb_) {
      electro_.reset(
          new Electro<EB>(m, *eb_, mebc_electro_, fs_->GetTime(), conf));
    } else {
      electro_.reset(new Electro<M>(m, m, mebc_electro_, fs_->GetTime(), conf));
    }
  }
}

template <class M>
void Hydro<M>::InitTracer(Multi<FieldCell<Scal>>& vfcu) {
  if (var.Int["enable_tracer"]) {
    auto multi = [](const std::vector<Scal>& v) {
      Multi<Scal> w(v.size());
      w.data() = v;
      return w;
    };
    typename TracerInterface<M>::Conf conf;
    conf.layers = var.Int["tracer_layers"];
    const auto tracer_layers = GRange<size_t>(conf.layers);

    conf.density = multi(var.Vect["tracer_density"]);
    fassert(conf.density.size() >= conf.layers);

    conf.viscosity = multi(var.Vect["tracer_viscosity"]);
    fassert(conf.viscosity.size() >= conf.layers);
    conf.gravity = Vect(var.Vect["gravity"]);

    auto termvel = multi(var.Vect["tracer_termvel"]);
    conf.diameter.resize(conf.layers);
    if (var.Int["tracer_use_termvel"]) {
      const Scal mixture_density = var.Double["rho1"];
      for (auto l : tracer_layers) {
        conf.diameter[l] = std::sqrt(std::abs(
            18 * conf.viscosity[l] * termvel[l] /
            (conf.gravity.norm() * (conf.density[l] - mixture_density))));
      }
    } else {
      conf.diameter = multi(var.Vect["tracer_diameter"]);
    }
    fassert(conf.diameter.size() >= conf.layers);

    conf.diffusion = multi(var.Vect["tracer_diffusion"]);
    fassert(conf.diffusion.size() >= conf.layers);

    conf.scheme = GetConvSc(var.String["tracer_scheme"]);

    using SlipType = typename TracerInterface<M>::SlipType;
    conf.slip.resize(conf.layers);
    for (auto l : tracer_layers) {
      auto& slip = conf.slip[l];
      std::stringstream arg(var.String["tracer" + std::to_string(l) + "_slip"]);
      std::string type;
      arg >> type;
      if (type == "none") {
        slip.type = SlipType::none;
      } else if (type == "stokes") {
        slip.type = SlipType::stokes;
      } else if (type == "termvel") {
        slip.type = SlipType::constant;
        slip.velocity = conf.gravity * (termvel[l] / conf.gravity.norm());
      } else if (type == "constant") {
        slip.type = SlipType::constant;
        arg >> slip.velocity;
      } else {
        fassert(false, "Unknown slip='" + type + "'");
      }
    }
    Multi<MapEmbed<BCond<Scal>>> vmebc(conf.layers); // boundary conditions
    me_group_.LoopPairs([&](auto cf_group) {
      auto cf = cf_group.first;
      size_t group = cf_group.second;
      auto nci = mebc_fluid_[cf].nci;

      const auto& custom = bc_group_custom_[group];
      auto getptr = [&](std::string key) -> const Scal* {
        auto it = custom.find(key);
        if (it != custom.end()) {
          return &it->second;
        }
        return nullptr;
      };

      for (auto l : tracer_layers) {
        auto sl = std::to_string(l);
        if (auto* dirichlet = getptr("tracer" + sl + "_dirichlet")) {
          vmebc[l][cf] = BCond<Scal>(BCondType::dirichlet, nci, *dirichlet);
        }
        if (auto* neumann = getptr("tracer" + sl + "_neumann")) {
          vmebc[l][cf] = BCond<Scal>(BCondType::neumann, nci, *neumann);
        }
      }
    });

    fc_tracer_source_.Reinit(tracer_layers, m, 0);
    conf.fc_source = fc_tracer_source_;

    if (eb_) {
      tracer_.reset(new Tracer<EB>(m, *eb_, vfcu, vmebc, fs_->GetTime(), conf));
    } else {
      tracer_.reset(new Tracer<M>(m, m, vfcu, vmebc, fs_->GetTime(), conf));
    }
  }
}

template <class M>
void Hydro<M>::InitFluid(const FieldCell<Vect>& fc_vel) {
  fcvm_ = fc_vel;

  std::string fs = var.String["fluid_solver"];
  if (fs == "simple") {
    auto par = ParsePar<Simple<M>>()(var);
    par.meshvel = st_.meshvel;
    std::shared_ptr<linear::Solver<M>> linsolver_symm(
        ULinear<M>::MakeLinearSolver(var, "symm", m));
    std::shared_ptr<linear::Solver<M>> linsolver_gen(
        ULinear<M>::MakeLinearSolver(var, "gen", m));
    const SimpleArgs<M> args{
        fc_vel,     mebc_fluid_,    mc_velcond_,   &fc_rho_,  &fc_mu_,
        &fc_force_, &febp_,         &fc_src_,      &fc_srcm_, 0.,
        st_.dt,     linsolver_symm, linsolver_gen, par};
    if (eb_) {
      fs_.reset(new Simple<Embed<M>>(m, *eb_, args));
    } else {
      fs_.reset(new Simple<M>(m, m, args));
    }
  } else if (fs == "proj") {
    auto par = ParsePar<Proj<M>>()(var);
    par.meshvel = st_.meshvel;
    std::shared_ptr<linear::Solver<M>> linsolver(
        ULinear<M>::MakeLinearSolver(var, "symm", m));
    linsolver_symm_ = linsolver;
    const ProjArgs<M> args{fc_vel,    mebc_fluid_, mc_velcond_, &fc_rho_,
                           &fc_mu_,   &fc_force_,  &febp_,      &fc_src_,
                           &fc_srcm_, 0.,          st_.dt,      linsolver,
                           par};
    if (eb_) {
      fs_.reset(new Proj<Embed<M>>(m, *eb_, args));
    } else {
      fs_.reset(new Proj<M>(m, m, args));
    }
  } else if (fs == "dummy") {
    if (eb_) {
      fs_.reset(new FluidDummy<Embed<M>>(
          m, *eb_, fc_vel, mebc_fluid_, &fc_rho_, &fc_mu_, &fc_force_, &febp_,
          &fc_src_, &fc_srcm_, 0., st_.dt, var));
    } else {
      fs_.reset(new FluidDummy<M>(
          m, m, fc_vel, mebc_fluid_, &fc_rho_, &fc_mu_, &fc_force_, &febp_,
          &fc_src_, &fc_srcm_, 0., st_.dt, var));
    }
  } else {
    fassert(false, "Unknown fluid_solver=" + fs);
  }
}

template <class M>
void Hydro<M>::InitAdvection(
    const FieldCell<Scal>& fcvf, const FieldCell<Scal>& fccl) {
  const std::string solver = var.String["advection_solver"];

  const auto modname = var.String["labeling"];
  if (auto* mod = ModuleLabeling<M>::GetInstance(modname)) {
    auto conf = ModuleLabeling<M>::ParseConf(var);
    labeling_ = mod->Make(conf, m);
  } else {
    std::string known;
    for (auto p : ModuleLabeling<M>::GetInstances()) {
      if (!known.empty()) {
        known += ", ";
      }
      known += p.first;
    }
    fassert(
        false,
        util::Format(
            "Unknown module labeling='{}' for connected component labeling. "
            "Known modules: {}",
            modname, known));
  }

  auto par = ParseVofPar<M>(var);
  par.labeling = labeling_.get();

  if (solver == "vof") {
    if (eb_) {
      as_.reset(new ASVEB(
          m, *eb_, fcvf, fccl, mebc_adv_, &fs_->GetVolumeFlux(Step::time_curr),
          &fc_src2_, 0., st_.dta, par));
    } else {
      as_.reset(new ASV(
          m, m, fcvf, fccl, mebc_adv_, &fs_->GetVolumeFlux(Step::time_curr),
          &fc_src2_, 0., st_.dta, par));
    }
    layers = GRange<size_t>(1);
  } else if (solver == "vofm") {
    if (eb_) {
      auto as = new ASVMEB(
          m, *eb_, fcvf, fccl, mebc_adv_, &fs_->GetVolumeFlux(Step::time_curr),
          &fc_src2_, 0., st_.dta, par);
      as_.reset(as);
      layers = GRange<size_t>(as->GetNumLayers());
    } else {
      auto as = new ASVM(
          m, m, fcvf, fccl, mebc_adv_, &fs_->GetVolumeFlux(Step::time_curr),
          &fc_src2_, 0., st_.dta, par);
      as_.reset(as);
      layers = GRange<size_t>(as->GetNumLayers());
    }
  } else {
    fassert(false, "Unknown advection_solver=" + solver);
  }

  fck_.resize(layers);
  fck_.InitAll(FieldCell<Scal>(m, GetNan<Scal>()));
  curv_estimator_ = curvature::MakeEstimator(var, m, layers);
}

template <class M>
template <class MEB>
void Hydro<M>::InitStat(const MEB& eb) {
  stat_.reset(new Stat<M>(m, eb_.get()));
  auto& stat = *stat_;
  auto& vf = as_->GetField();
  auto& vel = fs_->GetVelocity();
  auto& rho = fc_rho_;
  auto& p = fs_->GetPressure();
  // mm needed to prevent Wshadow
  stat.AddSum(
      "vol", "volume of domain", //
      [](IdxCell c, const M& mm) { return mm.GetVolume(c); });
  stat.AddSum(
      "vol1", "volume of phase 1", //
      [&vf](IdxCell c, const M& mm) { return (1 - vf[c]) * mm.GetVolume(c); });
  if (eb_) {
    stat.AddSum(
        "vol_eb", "volume of domain with embed", //
        [](IdxCell c, const EB& eeb) { return eeb.GetVolume(c); });
    stat.AddSum(
        "vol1_eb", "volume of phase 1 with embed", //
        [&vf](IdxCell c, const EB& eeb) {
          return eeb.GetVolume(c) - vf[c] * eeb.GetMesh().GetVolume(c);
        });
  }
  stat.AddSum(
      "vol2", "volume of phase 2", //
      [&vf](IdxCell c, const M& mm) { return vf[c] * mm.GetVolume(c); });
  stat.AddSum(
      "ekin", "kinetic energy of mixture", //
      [&vel, &rho](IdxCell c, const MEB& meb) {
        return 0.5 * rho[c] * vel[c].sqrnorm() * meb.GetVolume(c);
      });
  stat.AddSum(
      "ekin1", "kinetic energy of phase 1", //
      [&vf, &vel, &rho](IdxCell c, const MEB& meb) {
        return 0.5 * rho[c] * vel[c].sqrnorm() * meb.GetVolume(c) * (1 - vf[c]);
      });
  stat.AddSum(
      "ekin2", "kinetic energy of phase 2", //
      [&vf, &vel, &rho](IdxCell c, const MEB& meb) {
        return 0.5 * rho[c] * vel[c].sqrnorm() * meb.GetVolume(c) * vf[c];
      });
  stat.AddMax(
      "pmax", "maximum pressure", //
      [&p](IdxCell c, const M&) { return p[c]; });
  stat.AddMin(
      "pmin", "minimum pressure", //
      [&p](IdxCell c, const M&) { return p[c]; });
  stat.AddNone("t", "time", [&fs_ = fs_]() { return fs_->GetTime(); });
  stat.AddNone(
      "dt", "fluid time step", //
      [&fs_ = fs_]() { return fs_->GetTimeStep(); });
  stat.AddNone(
      "dta", "advection time step", //
      [&as_ = as_]() { return as_->GetTimeStep(); });
  stat.AddNone(
      "wt", "wall-clock time", //
      [&timer_ = timer_]() { return timer_.GetSeconds(); });
  stat.AddNone("iter", "total fluid iteration number", [&st_ = st_]() {
    return st_.iter;
  });
  stat.AddNone(
      "step", "fluid step number", [&st_ = st_]() { return st_.step; });
  stat.AddNone("meshvel", "moving mesh velocity", [&st_ = st_]() {
    return st_.meshvel;
  });
  stat.AddNone(
      "diff", "velocity difference between last two iterations", //
      [&fs_ = fs_]() { return fs_->GetError(); });

  stat.AddSumHidden("x*vol1", "", [&vf](IdxCell c, const M& mm) { //
    return mm.GetCenter(c) * mm.GetVolume(c) * (1 - vf[c]);
  });
  stat.AddSumHidden("x*vol2", "", [&vf](IdxCell c, const M& mm) { //
    return mm.GetCenter(c) * mm.GetVolume(c) * vf[c];
  });
  stat.AddSumHidden("v*vol1", "", [&vf, &vel](IdxCell c, const M& mm) { //
    return vel[c] * mm.GetVolume(c) * (1 - vf[c]);
  });
  stat.AddSumHidden("v*vol2", "", [&vf, &vel](IdxCell c, const M& mm) { //
    return vel[c] * mm.GetVolume(c) * vf[c];
  });
  stat.AddSumHidden("p*vol1", "", [&p, &vf](IdxCell c, const M& mm) { //
    return p[c] * (1 - vf[c]) * mm.GetVolume(c);
  });
  stat.AddSumHidden("p*vol2", "", [&p, &vf](IdxCell c, const M& mm) { //
    return p[c] * vf[c] * mm.GetVolume(c);
  });

  auto div = [](auto v, Scal d) {
    if (d == 0) {
      return v * 0;
    }
    return v / d;
  };

  auto& ffv = fs_->GetVolumeFlux();
  stat.AddSum(
      "q_inlet", "inlet volume rate", //
      [&ffv, this, &mm = m, &meb = eb]() {
        Scal sum = 0;
        mebc_fluid_.LoopBCond(meb, [&](auto cf, IdxCell c, auto bc) { //
          if (mm.IsInner(c)) {
            if (bc.type == BCondFluidType::inlet ||
                bc.type == BCondFluidType::inletflux) {
              sum += ffv[cf] * (bc.nci == 0 ? -1 : 1);
            }
          }
        });
        return sum;
      });
  stat.AddSum(
      "q_inletpressure", "inletpressure volume rate", //
      [&ffv, this, &mm = m, &meb = eb]() {
        Scal sum = 0;
        mebc_fluid_.LoopBCond(meb, [&](auto cf, IdxCell c, auto bc) { //
          if (mm.IsInner(c)) {
            if (bc.type == BCondFluidType::inletpressure) {
              sum += ffv[cf] * (bc.nci == 0 ? -1 : 1);
            }
          }
        });
        return sum;
      });
  stat.AddSum(
      "q_outlet", "outlet volume rate", //
      [&ffv, this, &mm = m, &meb = eb]() {
        Scal sum = 0;
        mebc_fluid_.LoopBCond(meb, [&](auto cf, IdxCell c, auto bc) { //
          if (mm.IsInner(c)) {
            if (bc.type == BCondFluidType::outlet ||
                bc.type == BCondFluidType::outletpressure) {
              sum += ffv[cf] * (bc.nci == 0 ? 1 : -1);
            }
          }
        });
        return sum;
      });
  stat.AddSum(
      "area_inletpressure", "inletpressure area", //
      [this, &mm = m, &meb = eb]() {
        Scal sum = 0;
        mebc_fluid_.LoopBCond(meb, [&](auto cf, IdxCell c, auto bc) { //
          if (mm.IsInner(c)) {
            if (bc.type == BCondFluidType::inletpressure) {
              sum += meb.GetArea(cf);
            }
          }
        });
        return sum;
      });
  stat.AddSumHidden(
      "p*area_inletpressure", "inletpressure pressure * area", //
      [&p, this, &mm = m, &meb = eb]() {
        Scal sum = 0;
        mebc_fluid_.LoopBCond(meb, [&](auto cf, IdxCell c, auto bc) { //
          if (mm.IsInner(c)) {
            if (bc.type == BCondFluidType::inletpressure) {
              sum += p[c] * meb.GetArea(cf);
            }
          }
        });
        return sum;
      });
  stat.AddDerived(
      "p_inletpressure", "inletpressure average pressure", //
      [div](const Stat<M>& s) {
        return div(s["p*area_inletpressure"], s["area_inletpressure"]);
      });

  stat.AddDerived(
      "c1", "centeroid of phase 1", //
      [div](const Stat<M>& s) { return div(s.vect["x*vol1"], s["vol1"]); });
  stat.AddDerived(
      "c2", "centeroid of phase 2", //
      [div](const Stat<M>& s) { return div(s.vect["x*vol2"], s["vol2"]); });
  stat.AddDerived(
      "v1", "velocity of phase 1", //
      [div](const Stat<M>& s) { return div(s.vect["v*vol1"], s["vol1"]); });
  stat.AddDerived(
      "v2", "velocity of phase 2", //
      [div](const Stat<M>& s) { return div(s.vect["v*vol2"], s["vol2"]); });
  stat.AddDerived(
      "pd", "pressure max-min", //
      [](const Stat<M>& s) { return s["pmax"] - s["pmin"]; });
  stat.AddDerived(
      "pavg1", "average pressure of phase 1", //
      [div](const Stat<M>& s) { return div(s["p*vol1"], s["vol1"]); });
  stat.AddDerived(
      "pavg2", "average pressure of phase 2", //
      [div](const Stat<M>& s) { return div(s["p*vol2"], s["vol2"]); });
  // XXX: relies on lexicographical order, "vol2_0" < "vol2_diff"
  //      to have "vol2_0" defined before "vol2_diff"
  stat.AddDerived(
      "vol2_0", "initial volume of phase 2", //
      [](const Stat<M>& s) {
        return s["vol2_0"] == 0 ? s["vol2"] : s["vol2_0"];
      });
  stat.AddDerived(
      "vol2_diff", "relative difference between vol2 and vol2_0", //
      [div](const Stat<M>& s) {
        return div(s["vol2"] - s["vol2_0"], s["vol2_0"]);
      });

  if (var.Int["stat_dissip"]) {
    stat.AddSum(
        "dissip", "dissipation rate of mixture", //
        [&str = fc_strain_, &mu = fc_mu_](IdxCell c, const M& mm) {
          return 2 * mu[c] * str[c] * mm.GetVolume(c);
        });
    stat.AddSum(
        "dissip1", "dissipation rate of phase 1", //
        [&str = fc_strain_, &mu = fc_mu_, &vf](IdxCell c, const M& mm) {
          return 2 * mu[c] * str[c] * (1 - vf[c]) * mm.GetVolume(c);
        });
    stat.AddSum(
        "dissip2", "dissipation rate of phase 2", //
        [&str = fc_strain_, &mu = fc_mu_, &vf](IdxCell c, const M& mm) {
          return 2 * mu[c] * str[c] * vf[c] * mm.GetVolume(c);
        });
    stat.AddDerived(
        "edis", "dissipated energy of mixture", //
        [&fs_ = fs_](const Stat<M>& s) {
          return s["edis"] + fs_->GetTimeStep() * s["dissip"];
        });
    stat.AddDerived(
        "edis1", "dissipated energy of phase 1", //
        [&fs_ = fs_](const Stat<M>& s) {
          return s["edis1"] + fs_->GetTimeStep() * s["dissip1"];
        });
    stat.AddDerived(
        "edis2", "dissipated energy of phase 2", //
        [&fs_ = fs_](const Stat<M>& s) {
          return s["edis2"] + fs_->GetTimeStep() * s["dissip2"];
        });
  }
  if (var.Int["enstrophy"]) {
    stat.AddSum(
        "enstr", "enstrophy", //
        [&omm = fcomm_, &rho = fc_rho_](IdxCell c, const M& mm) {
          return 0.5 * sqr(omm[c]) * rho[c] * mm.GetVolume(c);
        });
  }
  if (var.Int["statvel"]) {
    const Vect vel0(var.Vect["vel"]);
    stat.AddSumHidden("dv*dv*vol", "", [&vel0, &vel](IdxCell c, const M& mm) {
      auto dv = vel[c] - vel0;
      return dv * dv * mm.GetVolume(c);
    });
    stat.AddSumHidden(
        "dv*dv*vol2", "", [&vel0, &vel, &vf](IdxCell c, const M& mm) {
          auto dv = vel[c] - vel0;
          return dv * dv * vf[c] * mm.GetVolume(c);
        });
    stat.AddMax(
        "dvel_max", "maximum velocity difference with `vel`",
        [&vel0, &vel](IdxCell c, const M&) {
          auto dv = vel[c] - vel0;
          return dv.abs();
        });
    stat.AddDerived(
        "dvel_l2", "L2 norm of velocity difference with `vel`", //
        [](const Stat<M>& s) { return Sqrt(s.vect["dv*dv*vol"] / s["vol"]); });
    stat.AddMax(
        "dvel2_max", "maximum phase 2 velocity difference with `vel`",
        [&vel0, &vel, &vf](IdxCell c, const M&) {
          auto dv = (vel[c] - vel0) * vf[c];
          return dv.abs();
        });
    stat.AddDerived(
        "dvel2_l2", "L2 norm of phase 2 velocity difference with `vel`", //
        [](const Stat<M>& s) {
          return Sqrt(s.vect["dv*dv*vol2"] / s["vol2"]);
        });
  }
  if (var.Int["stat_vofm"]) {
    auto add_vofm = [this, &stat](auto as) {
      for (auto l : layers) {
        auto sl = std::to_string(l);
        const auto& vfl = *as->GetFieldM()[l];
        stat.AddSum(
            "vofm_cells_vf" + sl, "cells with positive volume fraction", //
            [&vfl](IdxCell c, const M&) -> Scal { return vfl[c] > 0 ? 1 : 0; });
        const auto& cl = *as->GetColor()[l];
        stat.AddSum(
            "vofm_cells_cl" + sl, "cells with defined color", //
            [&cl](IdxCell c, const M&) -> Scal { return cl[c] > 0 ? 1 : 0; });
        stat.AddSum(
            "vofm_vol" + sl, "integral of volume fraction", //
            [&vfl](IdxCell c, const M& mm) {
              return vfl[c] * mm.GetVolume(c);
            });
        auto slp = std::to_string(l + 1);
        // TODO: revise with `++hist[cnt]`, now takes L^2 operations.
        stat.AddSum(
            "vofm_hist" + slp,
            "number of cells with " + slp + " non-empty layers", //
            [vfm = as->GetFieldM(), &layers = layers, cnt_target = l + 1](
                IdxCell c, const M&) -> Scal {
              size_t cnt = 0;
              for (auto ll : layers) {
                if ((*vfm[ll])[c] > 0) {
                  ++cnt;
                }
              }
              return cnt == cnt_target ? 1 : 0;
            });
      }
    };
    if (auto as = dynamic_cast<ASVM*>(as_.get())) {
      add_vofm(as);
    }
    if (auto as = dynamic_cast<ASVMEB*>(as_.get())) {
      add_vofm(as);
    }
  }
  if (var.Int["stat_area"]) {
    if (auto as = dynamic_cast<ASV*>(as_.get())) {
      auto& fcn = as->GetNormal();
      auto& fca = as->GetAlpha();
      auto& fcvf = as->GetField();
      const Vect h = m.GetCellSize();
      stat.AddSum(
          "area", "area of the interface", //
          [&fcn, &fca, &fcvf, &h](IdxCell c, const M&) -> Scal {
            if (fcvf[c] > 0 && fcvf[c] < 1 && !IsNan(fca[c])) {
              using R = Reconst<Scal>;
              auto s = std::abs(
                  R::GetArea(R::GetCutPoly2(fcn[c], fca[c], h), fcn[c]));
              if (!IsNan(s)) {
                return s;
              }
            }
            return 0;
          });
    }
  }
  if (eb_) {
    stat.AddSum("pdrag", "pressure drag", [this]() {
      auto d = CalcPressureDrag(fs_->GetPressure(), *eb_);
      if (m.GetEdim() == 2 && M::dim > 2) {
        d /= m.GetCellSize()[2];
      }
      return d;
    });
    stat.AddSum("vdrag", "viscous drag", [this]() {
      auto d = CalcViscousDrag(fs_->GetVelocity(), fc_mu_, *eb_);
      if (m.GetEdim() == 2 && M::dim > 2) {
        d /= m.GetCellSize()[2];
      }
      return d;
    });
    stat.AddDerived("drag", "total drag", [](const Stat<M>& s) {
      return s.vect["pdrag"] + s.vect["vdrag"];
    });
  }
  if (particles_) {
    stat.AddSum(
        "particles_n", "number of particles",
        [&p = particles_]() -> Scal { return p->GetView().x.size(); });
    stat.AddSum(
        "particles_nrecv", "number of particles transferred at last step",
        [&p = particles_]() -> Scal { return p->GetNumRecv(); });
  }
  if (electro_) {
    stat.AddNone(
        "el_current", "total electric current", //
        [&electro_ = electro_]() { return electro_->GetStat().current; });
    stat.AddNone(
        "el_potential", "potential difference between min and max", //
        [&electro_ = electro_]() { return electro_->GetStat().potential; });
    if (var.Int["enable_electro_control"]) {
      stat.AddNone(
          "el_potential_control", "electric potential from controller", //
          [&st_ = st_]() { return st_.electro_control_potential; });
    }
  }
  if (tracer_) {
    for (auto l : GRange<size_t>(tracer_->GetConf().layers)) {
      const auto sl = std::to_string(l);
      stat.AddSum(
          "tu" + sl + "_vol", "volume of tracer " + sl, //
          [&tr = tracer_, l](IdxCell c, const MEB& meb) {
            return tr->GetVolumeFraction()[l][c] * meb.GetVolume(c);
          });
      stat.AddMax(
          "tu" + sl + "_max", "maximum of tracer " + sl, //
          [&tr = tracer_, l](IdxCell c, const MEB&) {
            return tr->GetVolumeFraction()[l][c];
          });
      stat.AddMin(
          "tu" + sl + "_min", "minimum of tracer " + sl, //
          [&tr = tracer_, l](IdxCell c, const MEB&) {
            return tr->GetVolumeFraction()[l][c];
          });
      stat.AddSum(
          "tu" + sl + "_outrate", "outlet volumetric rate of tracer " + sl, //
          [&ffv, &fcu = tracer_->GetVolumeFraction()[l], this, &meb = eb, l]() {
            Scal sum = 0;
            mebc_fluid_.LoopBCond(meb, [&, this](auto cf, IdxCell c, auto bc) {
              if (m.IsInner(c)) {
                if (bc.type == BCondFluidType::outlet ||
                    bc.type == BCondFluidType::outletpressure) {
                  sum += fcu[c] * ffv[cf] * (bc.nci == 0 ? 1 : -1);
                }
              }
            });
            return sum;
          });
      stat.AddDerived(
          "tu" + sl + "_outvol", "outlet volume of tracer " + sl, //
          [this, sl](const Stat<M>& s) {
            return s["tu" + sl + "_outvol"] +
                   s["tu" + sl + "_outrate"] * st_.dt;
          });
    }
    stat.AddNone(
        "total_src2", "total source of component 2", //
        [& st_ = st_]() { return st_.total_src2; });
  }

  stat_->SortNames();

  auto blacklist = GetWords(var.String("stat_blacklist", ""));
  for (auto name : stat_->GetNames()) {
    stat_->SetEnabled(name, !blacklist.count(name));
  }

  if (m.IsRoot() && dumpstat_) {
    {
      const std::string path = "stat.dat";
      fstat_.open(path);
      fassert(fstat_.good(), "Can't open file '" + path + "' for writing");
      fstat_.precision(16);
      stat_->WriteHeader(fstat_);
    }

    {
      const std::string path = "stat_summary";
      std::ofstream fsum(path);
      fassert(fsum.good(), "Can't open file '" + path + "' for writing");
      stat_->WriteSummary(fsum, true);
    }
  }
}

template <class M>
void Hydro<M>::Init() {
  using namespace fluid_condition;
  using Xmf = dump::Xmf<Vect>;
  auto sem = m.GetSem("init");
  struct {
    FieldCell<Vect> fcvel; // initial velocity
    FieldCell<Scal> fcvf; // initial volume fraction
    FieldCell<Scal> fccl; // initial color
    FieldCell<Vect> fcpot; // velocity potential (stream function)
    FieldCell<Scal> fcpot_scal; // first componenet of velocity potential
    typename Xmf::Meta meta; // metadata for velocity potential dump
    Multi<FieldCell<Scal>> tracer_vfcu;
    std::shared_ptr<linear::Solver<M>> linsolver_vort;
    FieldCell<bool> fc_innermask;
    FieldFace<bool> ff_innermask;
    std::vector<std::pair<std::string, std::vector<Scal>>> nucl_data;
  } * ctx(sem);
  auto& t = *ctx;
  auto& fcvel = t.fcvel;
  auto& fcvf = t.fcvf;
  auto& fccl = t.fccl;
  if (sem("flags")) {
    silent_ = var.Int("silent", 0);
    dumpstat_ = var.Int("dumpstat", 1);
    if (m.IsRoot() && var.Int("dumpconfig", 1)) {
      std::ofstream out("out.conf");
      Parser::PrintVars(var, out);
    }

    m.flags.linreport = var.Int["linreport"];
    m.flags.check_symmetry = var.Int["check_symmetry"];
    m.flags.check_symmetry_dump_threshold =
        var.Double["check_symmetry_dump_threshold"];
    randgen_.seed(m.GetId() + 1);

    if (auto* name = var.String.Find("module_post_step")) {
      module_post_step_ = ModulePostStep<M>::GetInstance(*name);
      fassert(
          module_post_step_,
          "ModulePostStep instance '" + *name + "' not found");
    }
  }
  if (sem.Nested("sysinfo") && !silent_) {
    ReportSysinfo(std::cerr);
  }
  if (sem.Nested("embed")) {
    InitEmbed();
  }
  if (sem.Nested()) {
    InitVf(fcvf, var, m, !silent_);
  }
  if (sem.Nested("stepwise") && var.Int["enable_stepwise_body"]) {
    InitStepwiseBody(t.fc_innermask, t.ff_innermask);
  }
  if (sem.Nested()) {
    if (var.Int["enable_tracer"]) {
      InitTracerFields(t.tracer_vfcu);
    }
  }
  if (sem("fields")) {
    if (!t.fc_innermask.empty()) {
      // XXX: Must be set before InitBc().
      m.flags.SetInnerMask(t.fc_innermask, t.ff_innermask);
    }
    if (eb_) {
      auto& eb = *eb_;
      for (auto c : m.AllCells()) {
        if (eb.GetType(c) == M::Type::excluded) {
          fcvf[c] = 0;
        }
      }
    }

    fc_src_.Reinit(m, 0.);
    fc_src2_.Reinit(m, 0.);
    fc_srcm_.Reinit(m, 0.);

    // initial surface tension sigma
    fc_sig_.Reinit(m, 0);
    CreateInitSig<M>(var)(fc_sig_, m);
    m.Comm(&fc_sig_);

    // initial contact angle
    {
      fc_contang_.Reinit(m, -1);
      const std::string name = var.String["init_contang"];
      if (auto ptr = ModuleInitContang<M>::GetInstance(name)) {
        (*ptr)(fc_contang_, var, m);
      } else {
        if (m.IsRoot()) {
          std::cerr << "Known values of 'init_contang': ";
          for (auto& p : ModuleInitContang<M>::GetInstances()) {
            std::cerr << p.first << " ";
          }
          std::cerr << std::endl;
        }
        fassert(false, "Unknown init_contang=" + name);
      }
      m.Comm(&fc_contang_);
    }

    // initial velocity
    fcvel.Reinit(m, Vect(0));
    InitVel(fcvel, var, m);
    if (eb_) {
      InitVelHook(fcvel, var, m, *eb_);
    } else {
      InitVelHook(fcvel, var, m);
    }
    m.Comm(&fcvel);
    fcvel.SetHalo(2);
    fcvel.SetName("fcvel");

    if (m.IsRoot() && !silent_) {
      std::cerr << "global mesh=" << m.GetGlobalSize() << '\n';
      std::cerr << "surface tension dt=" << GetSurfaceTensionDt() << '\n';
      std::cerr << "viscosity dt=" << GetViscosityDt() << '\n';
    }

    // Known names of keys for custom boundary conditions.
    std::set<std::string> known_keys = {
        "electro_dirichlet", "electro_neumann",   "tracer0_dirichlet",
        "tracer0_neumann",   "tracer1_dirichlet", "tracer1_neumann",
        "nucleation", // factor to `nucleation_number_density`
        "electro_control", // factor to potential from controller
    };
    if (eb_) {
      auto initbc = InitBc(var, *eb_, known_keys, t.fc_innermask);
      mebc_fluid_ = std::get<0>(initbc);
      mebc_adv_ = std::get<1>(initbc);
      me_group_ = std::get<2>(initbc);
      bc_group_desc_ = std::get<3>(initbc);
      bc_group_custom_ = std::get<4>(initbc);
      mebc_adv_.LoopBCond(*eb_, [&](auto cf, auto c, auto) { //
        me_contang_[cf] = fc_contang_[c];
        mebc_adv_[cf].contang = fc_contang_[c];
      });
    } else {
      auto initbc = InitBc(var, m, known_keys, t.fc_innermask);
      mebc_fluid_ = std::get<0>(initbc);
      mebc_adv_ = std::get<1>(initbc);
      me_group_ = std::get<2>(initbc);
      bc_group_desc_ = std::get<3>(initbc);
      bc_group_custom_ = std::get<4>(initbc);
      for (auto& p : mebc_adv_.GetMapFace()) {
        const auto& f = p.first;
        auto& bc = p.second;
        const auto c = m.GetCell(f, bc.nci);
        me_contang_[f] = fc_contang_[c];
        bc.contang = fc_contang_[c];
      }
    }
    mebc_fluid_orig_ = mebc_fluid_;

    if (var.Int["enable_electro"]) {
      me_group_.LoopPairs([&](auto cf_group) {
        auto cf = cf_group.first;
        size_t group = cf_group.second;
        auto nci = mebc_fluid_[cf].nci;

        const auto& custom = bc_group_custom_[group];
        auto getptr = [&](std::string key) -> const Scal* {
          auto it = custom.find(key);
          if (it != custom.end()) {
            return &it->second;
          }
          return nullptr;
        };

        if (auto* dirichlet = getptr("electro_dirichlet")) {
          mebc_electro_[cf] =
              BCond<Scal>(BCondType::dirichlet, nci, *dirichlet);
        } else if (auto* neumann = getptr("electro_neumann")) {
          mebc_electro_[cf] = BCond<Scal>(BCondType::neumann, nci, *neumann);
        } else {
          fassert(
              false,
              "Unknown electro conditions for group" + std::to_string(group));
        }
      });
      if (var.Int["enable_electro_control"]) {
        me_group_.LoopPairs([&](auto cf_group) {
          auto cf = cf_group.first;
          size_t group = cf_group.second;
          const auto& custom = bc_group_custom_[group];
          const auto it = custom.find("electro_control");
          if (it != custom.end()) {
            electro_control_factor_[cf] = it->second;
          }
        });
        st_.electro_control_potential =
            var.Double["electro_control_initial_potential"];
      }
    }

    // boundary conditions for smoothing of volume fraction
    mebc_vfsm_ = GetBCondZeroGrad<Scal>(mebc_fluid_);
  }

  if (var.Int["bc_wall_init_vel"] && sem("bc_wall_init_vel")) {
    // velocity on walls from initial conditions in neighbor cells
    for (auto& p : mebc_fluid_.GetMapFace()) {
      const IdxFace f = p.first;
      auto& bc = p.second;
      if (bc.type == BCondFluidType::wall) {
        const IdxCell c = m.GetCell(f, bc.nci);
        bc.velocity = fcvel[c];
      }
    }
  }

  if (var.Int["initvort"]) {
    if (sem()) {
      t.linsolver_vort = ULinear<M>::MakeLinearSolver(var, "vort", m);
    }
    if (sem.Nested("initvort")) {
      InitVort(
          fcvel, fcvel, &t.fcpot, mebc_fluid_, t.linsolver_vort, m,
          var.Int("initvort_zero_dirichlet", 0));
    }
    if (var.Int("initvort_dumppot", 0)) {
      if (sem.Nested()) {
        dump::Raw<M>::WriteWithXmf(t.fcpot_scal, "velpot", "velpot.raw", m);
      }
    }
  }

  if (var.Int["vel_init_random"] && sem("random")) {
    Scal amp = var.Double["random_amp"];
    Vect vel(var.Vect["random_vel"]);
    std::default_random_engine g(m.GetId());
    std::uniform_real_distribution<Scal> u(-amp, amp);
    for (auto c : m.Cells()) {
      Vect v = vel * u(g) / 7;
      for (auto q : m.Nci(c)) {
        IdxCell cn = m.GetCell(c, q);
        fcvel[cn] += v;
      }
      fcvel[c] += v;
    }
    m.Comm(&fcvel);
  }

  if (sem.Nested("smooth")) {
    Smoothen(fcvf, mebc_vfsm_, m, var.Int["vf_init_sm"]);
  }

  if (sem.Nested("mixture")) {
    CalcMixture(fcvf);
  }

  if (sem("color-ini")) {
    if (var.Int["enable_color"]) {
      // initial color
      // TODO revise with bcast
      auto icl = CreateInitCl<M>(var, m.IsRoot());
      icl(fccl, fcvf, m);
      m.Comm(&fccl);
    } else {
      fccl.Reinit(m, 0.);
    }
  }

  if (sem.Nested("cellcond")) {
    GetFluidCellCond(var, m, mc_velcond_);
  }
  if (sem("calcdt0")) {
    const Scal dt = var.Double["dt0"];
    st_.dt = dt;
    st_.dta = dt;
    tracer_dt_ = dt;
    particles_dt_ = dt;
  }
  if (sem("solv")) {
    InitFluid(fcvel);

    InitAdvection(fcvf, fccl);

    InitTracer(t.tracer_vfcu);

    InitElectro();

    if (var.Int["enable_particles"]) {
      InitParticles();
    }

    if (tracer_ && var.Int["enable_nucleation_points"]) {
      InitNucleationPoints();
    }

    st_.iter = 0;
    st_.step = 0;

    if (m.IsLead()) {
      this->var_mutable.Int.Set("iter", st_.iter);
    }

    if (eb_) {
      InitStat(*eb_);
    } else {
      InitStat(m);
    }

    if (var.Int["fill_halo_nan"]) {
      std::vector<std::pair<IdxFace, size_t>> vf;
      for (auto& p : mebc_adv_.GetMapFace()) {
        vf.emplace_back(p.first, p.second.GetNci());
      }
      m.SetNanFaces(vf);
      m.flags.nan_faces_value = var.Double["fill_halo_nan_value"];
    }

    events_ = std::unique_ptr<Events>(
        new Events(this->var_mutable, m.IsRoot(), m.IsLead(), !silent_));
    events_->AddHandler(
        "vf_save_state",
        [&path = vf_save_state_path_](std::string arg) { //
          path = arg;
        });
    events_->Parse();
  }
  if (sem.Nested("inithook")) {
    InitHook(this);
  }
  if (sem.Nested("vofm-load")) {
    if (var.Int["init_vf_load_state"]) {
      if (auto as = dynamic_cast<ASVMEB*>(as_.get())) {
        as->LoadState(var.String["init_vf_state_dir"]);
      }
    }
  }
  if (var.Int["dumpbc"]) {
    if (bc_group_desc_.size()) {
      if (sem("dump-bcgroups")) {
        if (m.IsRoot()) {
          std::ofstream fdesc("bc_groups.dat");
          for (size_t i = 0; i < bc_group_desc_.size(); ++i) {
            fdesc << i << " " << bc_group_desc_[i] << std::endl;
          }
        }
      }
      if (sem.Nested("bcdump")) {
        if (eb_) {
          DumpBcPoly("bc.vtk", me_group_, me_contang_, *eb_, m);
        } else {
          DumpBcPoly("bc.vtk", me_group_, me_contang_, m, m);
        }
      }
    }
  }
  if (var.Int["enable_nucleation_points"] && var.Int("dump_nucl", 0)) {
    const std::string path = "nucl.csv";
    if (sem("dump-nucl-local")) {
      if (m.IsRoot()) {
        std::cerr << util::Format("dump nucleation points to {}\n", path);
      }
      t.nucl_data.push_back({"x", {}});
      t.nucl_data.push_back({"y", {}});
      t.nucl_data.push_back({"z", {}});
      for (auto xd : nucl_points_) {
        using Vect3 = generic::Vect<Scal, 3>;
        Vect3 x(xd);
        t.nucl_data[0].second.push_back(x[0]);
        t.nucl_data[1].second.push_back(x[1]);
        t.nucl_data[2].second.push_back(x[2]);
      }
    }
    if (sem.Nested()) {
      dump::DumpCsv(t.nucl_data, path, m);
    }
  }
  if (eb_ && sem.Nested() && var.Int("dump_eb", 1)) {
    eb_->DumpPoly("eb.vtk", var.Int["vtkbin"], var.Int["vtkmerge"]);
  }
  if (sem.Nested("stat") && dumpstat_) {
    stat_->Update();
  }
  if (var.Int["dumpinit"]) {
    if (sem.Nested()) {
      Dump(true);
    }
  }
  if (sem()) {
    initialized_ = true;
  }
}

template <class M>
Hydro<M>::Hydro(Vars& var_, const BlockInfoProxy& bi, Par& par)
    : KernelMeshPar<M, Par>(var_, bi, par)
    , dumper_(var, "dump_field_")
    , dmptraj_(var, "dump_traj_")
    , dmptrep_(var, "dump_trep_")
    , bubgen_(var, "bubgen_") {}

template <class M>
void Hydro<M>::CalcStat() {
  auto sem = m.GetSem("stat");
  if (sem("local")) {
    if (var.Int["stat_dissip"]) {
      fc_strain_ = CalcStrain(fs_->GetVelocity());
    }
    if (var.Int["enstrophy"]) {
      CalcVort();
    }
  }
}

template <class M>
void Hydro<M>::CalcDt() {
  auto sem = m.GetSem("calcdt");
  struct {
    Scal dtmin;
    Scal tmin;
  } * ctx(sem);
  auto& t = *ctx;

  if (sem("local")) {
    st_.t = fs_->GetTime();
    t.dtmin = fs_->GetAutoTimeStep();
    m.Reduce(&t.dtmin, Reduction::min);
  }
  if (sem("reduce")) {
    // Determine fluid time step.
    // Set from cfl if defined.
    if (auto* cfl = var.Double.Find("cfl")) {
      st_.dt = t.dtmin * (*cfl);
      st_.dt = std::min<Scal>(st_.dt, var.Double["dtmax"]);
    }

    { // Constraint from surface tension.
      const auto dt = GetSurfaceTensionDt();
      if (dt > 0) {
        st_.dt = std::min<Scal>(st_.dt, dt);
      }
    }

    { // Constraint from viscosity.
      const auto dt = GetViscosityDt();
      if (dt > 0) {
        st_.dt = std::min<Scal>(st_.dt, dt);
      }
    }
    // FIXME: Without this reduction the time step turns out slighly
    //        different on various blocks. Which shouldn't happen.
    m.Reduce(&st_.dt, Reduction::min);
  }
  if (sem("reduce2")) {
    fs_->SetTimeStep(st_.dt);

    // Determine advection time step.
    // Set from cfla if defined.
    if (auto* cfla = var.Double.Find("cfla")) {
      st_.dta = t.dtmin * (*cfla);
      st_.dta = std::min<Scal>(st_.dta, var.Double["dtmax"]);
    }
    // Round up dta to such that dt / dta is integer.
    const Scal dt = fs_->GetTime() + fs_->GetTimeStep() - as_->GetTime();
    st_.dta = dt / std::max(1, int(dt / st_.dta + 0.5));
    as_->SetTimeStep(st_.dta);

    // Update tracer time step.
    if (tracer_) {
      tracer_dt_ = fs_->GetTimeStep();
      if (auto* cfl = var.Double.Find("cfl_tracer_advection")) {
        tracer_dt_ = std::min(tracer_dt_, t.dtmin * (*cfl));
      }
      if (auto* cfl = var.Double.Find("cfl_tracer_diffusion")) {
        Scal diffusion = 0;
        for (auto d : var.Vect["tracer_diffusion"]) {
          diffusion = std::max(diffusion, d);
        }
        if (diffusion > 0) {
          tracer_dt_ = std::min<Scal>(
              tracer_dt_, (*cfl) * sqr(m.GetCellSize()[0]) / diffusion);
        }
      }
      // Round up such that fluid time step is a multiple of `tracer_dt_`.
      const Scal dtwhole =
          fs_->GetTime() - tracer_->GetTime() + fs_->GetTimeStep();
      tracer_dt_ = dtwhole / std::max<Scal>(1, int(dtwhole / tracer_dt_ + 0.5));
    }

    if (particles_) {
      if (auto* cflp = var.Double.Find("cflp")) {
        particles_dt_ = t.dtmin * (*cflp);
        // round up dta to such that dt / dta is integer
        const Scal dtwhole =
            fs_->GetTime() + fs_->GetTimeStep() - particles_->GetTime();
        particles_dt_ =
            dtwhole / std::max(1, int(dtwhole / particles_dt_ + 0.5));
      } else {
        particles_dt_ = fs_->GetTimeStep();
      }
    }
    t.tmin = fs_->GetTime();
    m.Reduce(&t.tmin, Reduction::min);
  }
  if (sem()) {
    fassert_equal(fs_->GetTime(), t.tmin);
  }
}

template <class M>
IdxCell GetCell(IdxFace f, Side nci, const M& m) {
  return m.GetCell(f, nci);
}

template <class M>
IdxCell GetCell(IdxCell c, Side, const M&) {
  return c;
}

template <class M>
void Hydro<M>::CalcMixture(const FieldCell<Scal>& fc_vf0) {
  auto sem = m.GetSem("mixture");

  if (sem("init")) {
    fc_mu_.Reinit(m);
    fc_rho_.Reinit(m);
    fc_force_.Reinit(m, Vect(0));
    febp_.Reinit(m, 0);
    fc_src_.Reinit(m, 0);
    fc_src2_.Reinit(m, 0);
    fc_smvf_ = fc_vf0;
    if (eb_ && var.Int["vfsmooth_extrapolate_cut"]) {
      auto& eb = *eb_;
      for (auto c : eb.CFaces()) {
        fc_smvf_[c] = fc_vf0[eb.GetRegularNeighbor(c)];
      }
      m.Comm(&fc_smvf_);
    }

    // update surface tension sigma
    CreateInitSig<M>(var)(fc_sig_, m);
    m.Comm(&fc_sig_);
  }

  if (sem.Nested("smooth")) {
    if (eb_ && var.Int["vfsmooth_extrapolate_cut"]) {
      Smoothen(fc_smvf_, mebc_vfsm_, *eb_, var.Int["vfsmooth"]);
    } else {
      Smoothen(fc_smvf_, mebc_vfsm_, m, var.Int["vfsmooth"]);
    }
  }

  if (sem("calc")) {
    FieldCell<Scal>& fcvfsm = fc_smvf_;
    FieldFace<Scal>& ffvf = ff_smvf_;
    if (eb_) {
      ffvf = UEB::Interpolate(fcvfsm, mebc_vfsm_, *eb_).GetFieldFace();
    } else {
      ffvf = UEB::Interpolate(fcvfsm, mebc_vfsm_, m);
    }

    const Vect force(var.Vect["force"]);
    const Vect grav(var.Vect["gravity"]);
    const Scal rho1(var.Double["rho1"]);
    const Scal rho2(var.Double["rho2"]);
    const Scal mu1(var.Double["mu1"]);
    const Scal mu2(var.Double["mu2"]);

    // Init density and viscosity
    for (auto c : m.AllCells()) {
      const Scal a2 = fcvfsm[c];
      const Scal a1 = 1 - a2;
      fc_rho_[c] = rho1 * a1 + rho2 * a2;
      fc_mu_[c] = mu1 * a1 + mu2 * a2;
    }
    FieldFace<Scal> ff_rho(m);
    for (auto f : m.AllFaces()) {
      const Scal a2 = ffvf[f];
      const Scal a1 = 1 - a2;
      ff_rho[f] = rho1 * a1 + rho2 * a2;
    }

    if (tracer_ && var.Int["tracer_override_mixture"]) {
      const auto& fc_rho_mix = tracer_->GetMixtureDensity();
      const auto& fc_mu_mix = tracer_->GetMixtureViscosity();
      const auto& fc_vf_mix = tracer_->GetMixtureVolumeFraction();
      for (auto c : m.AllCells()) {
        const Scal rho1m = fc_rho_mix[c] + (1 - fc_vf_mix[c]) * rho1;
        const Scal mu1m = fc_mu_mix[c] + (1 - fc_vf_mix[c]) * mu1;
        const Scal a2 = fcvfsm[c];
        const Scal a1 = 1 - a2;
        fc_rho_[c] = rho1m * a1 + rho2 * a2;
        fc_mu_[c] = mu1m * a1 + mu2 * a2;
      }
      if (eb_) {
        ff_rho = UEB::Interpolate(fc_rho_, mebc_vfsm_, *eb_).GetFieldFace();
      } else {
        ff_rho = UEB::Interpolate(fc_rho_, mebc_vfsm_, m);
      }
    }

    // Append gravity to force
    for (auto f : m.AllFaces()) {
      const Vect n = m.GetNormal(f);
      febp_[f] += force.dot(n);
      febp_[f] += grav.dot(n) * ff_rho[f];
    }

    // Surface tension
    if (var.Int["enable_surftens"] && as_) {
      CalcSurfaceTension(
          m, layers, var, fc_force_, febp_.GetFieldFace(), fc_sig_,
          GetBCondZeroGrad<Scal>(mebc_fluid_), fck_, fc_vf0, ffvf, as_.get());
    }

    // zero force in z if 2D
    const size_t edim = var.Int["dim"];
    if (edim < M::dim) {
      for (auto f : m.Faces()) {
        using Dir = typename M::Dir;
        if (m.GetIndexFaces().GetDir(f) == Dir(2)) {
          febp_[f] = 0;
        }
      }
    }

    // Source of uniform bubble growth.
    if (auto* rate = var.Double.Find("uniform_growth_rate")) {
      if (auto as = dynamic_cast<ASV*>(as_.get())) {
        auto& fcn = as->GetNormal();
        auto& fca = as->GetAlpha();
        auto& fcvf = as->GetField();
        const Vect h = m.GetCellSize();
        for (auto c : m.Cells()) {
          if (fcvf[c] > 0 && fcvf[c] < 1 && !IsNan(fca[c])) {
            using R = Reconst<Scal>;
            auto area =
                std::abs(R::GetArea(R::GetCutPoly2(fcn[c], fca[c], h), fcn[c]));
            auto src2 = (*rate) * area / m.GetVolume(c);
            fc_src2_[c] += src2;
            fc_src_[c] += src2;
          }
        }
      }
    }
    // Clear tracer source.
    if (tracer_) {
      const auto& conf = tracer_->GetConf();
      for (auto l : GRange<size_t>(conf.layers)) {
        fc_tracer_source_[l].Reinit(m, 0);
      }
    }
    // Source of bubble growth from tracers.
    if (tracer_) {
      const auto& conf = tracer_->GetConf();
      for (auto l : GRange<size_t>(conf.layers)) {
        const auto sl = std::to_string(l);
        const Scal rate = var.Double("growth_rate_tracer" + sl, 0);
        if (rate) {
          auto apply = [&](const auto& meb, auto as) {
            if (!as) return;
            const Scal dt = st_.dt;
            const auto& fctu = tracer_->GetVolumeFraction()[l];
            const auto& fcvf = as->GetField();
            for (auto c : meb.Cells()) {
              if (meb.IsRegular(c)) {
                Scal src2 = fctu[c] * rate * fcvf[c] / dt;
                src2 = std::min(src2, std::max(0., fctu[c] / dt));
                src2 = std::max(src2, -fcvf[c] / dt);
                fc_src2_[c] += src2;
                fc_tracer_source_[l][c] -= src2;
                fc_src_[c] += src2;
              }
            }
            for (auto c : nucl_cells_) {
              if (meb.IsRegular(c)) {
                Scal src2 = fctu[c] * rate * fcvf[c] / dt;
                src2 = std::min(src2, std::max(0., fctu[c] / dt));
                src2 = std::max(src2, 0.);
                fc_src2_[c] += src2;
                fc_tracer_source_[l][c] -= src2;
                fc_src_[c] += src2;
              }
            }
          };
          if (eb_) {
            apply(*eb_, dynamic_cast<ASVEB*>(as_.get()));
          } else {
            apply(m, dynamic_cast<ASV*>(as_.get()));
          }
        }
      }
    }
    if (tracer_) {
      const auto& conf = tracer_->GetConf();
      for (auto l : GRange<size_t>(conf.layers)) {
        const Scal factor = var.Double("tracer_interface_diffusion_factor", 0);
        if (factor) {
          auto apply = [&](const auto& meb, auto as) {
            if (!as) return;
            const Scal dt = st_.dt;
            const auto& fctu = tracer_->GetVolumeFraction()[l];
            const auto& fcn = as->GetNormal();
            const auto& fca = as->GetAlpha();
            const auto& fcvf = as->GetField();
            const Vect h = m.GetCellSize();
            for (auto c : meb.CellsM()) {
              if (meb.IsRegular(c) && fcvf[c] > 0 && fcvf[c] < 1 &&
                  !IsNan(fca[c])) {
                using R = Reconst<Scal>;
                auto area = std::abs(
                    R::GetArea(R::GetCutPoly2(fcn[c], fca[c], h), fcn[c]));
                const Scal d = conf.diffusion[l] * factor;
                Scal src2 = fctu[c] * d * area / (h[0] * c.volume);
                src2 = std::min(src2, std::max(0., fctu[c] / dt));
                src2 = std::max(src2, -fcvf[c] / dt);
                fc_src2_[c] += src2;
                fc_tracer_source_[l][c] -= src2;
                fc_src_[c] += src2;
              }
            }
          };
          if (eb_) {
            apply(*eb_, dynamic_cast<ASVEB*>(as_.get()));
          } else {
            apply(m, dynamic_cast<ASV*>(as_.get()));
          }
        }
      }
    }
    // Growth of particles from tracer.
    if (tracer_ && particles_) {
      const auto& conf = tracer_->GetConf();
      auto view = particles_->GetView();
      for (auto l : GRange<size_t>(conf.layers)) {
        const auto& fctu = tracer_->GetVolumeFraction()[l];
        const auto h = m.GetCellSize();
        const Scal factor = var.Double["particles_growth_diffusion_factor"];
        const Scal rate = var.Double["particles_growth_rate"];
        const auto& blockall = m.GetAllBlockCells();
        const MIdx blockall_wmin = blockall.GetBegin();
        const MIdx blockall_wmax = blockall.GetEnd() - MIdx(1);
        // Traverse inner and halo particles.
        // Accumulate the particle sources in `view.source[i]`
        // and field sources in `fc_tracer_source_[l][c]`.
        // Particles sources are only valid for inner particles
        // since they depend on the 5x5x5 stencil only available in inner cells.
        // To ensure that `fc_tracer_source_` is correctly updated from halo
        // particles, the particle source must be a sum of contributions from
        // cells, i.e. sum of values for each pair (i,c).
        for (size_t i = 0; i < view.x.size(); ++i) {
          view.source[i] = 0;
          // Range of cells intersecting with the particle.
          MIdx wmin(((view.x[i] - Vect(view.r[i]) - m.flags.global_origin) / h)
                        .floor());
          MIdx wmax(((view.x[i] + Vect(view.r[i]) - m.flags.global_origin) / h)
                        .floor());
          // Cell containing the particle centroid.
          const auto cpart = m(m.GetCellFromPoint(view.x[i]));
          // Limit the range of cells to 5x5x5 stencil around the centroid.
          const MIdx wcpart(cpart);
          wmin = wmin.max(wcpart - MIdx(2));
          wmax = wmax.min(wcpart + MIdx(2));
          // Limit the range of cells to inner and halo cells of current block.
          wmin = wmin.max(blockall_wmin);
          wmax = wmax.min(blockall_wmax);
          // Block of cells from range [wmin, wmax].
          const GBlock<int, M::dim> wcells(wmin, wmax - wmin + MIdx(1));
          const auto& index = m.GetIndexCells();
          // Level-set function positive inside the particle.
          auto ls = [&](const Vect& x) -> Scal {
            return view.r[i] - x.dist(view.x[i]);
          };
          for (MIdx w : wcells) {
            const auto c = m(index.GetIdx(w));
            // Volume fraction of intersection between particle and cell.
            const Scal vf = //
                view.r[i] < h[0]
                    ? view.r[i] / h[0] * 0.5
                    : GetLevelSetVolume<Scal, M::dim>(ls, c.center, h);
            const Scal u = fctu[c]; // Concentration.
            if (u > 0 && vf > 0) {
              const Scal d = conf.diffusion[l] * factor;
              Scal src = d / c.volume * u * vf * h[0] * rate;
              src = std::min(src, std::max(0., u / tracer_dt_));

              if (view.inner[i]) {
                view.source[i] += src * c.volume;
              }
              fc_tracer_source_[l][c] -= src;
              fc_src_[c] += src;
            }
          }
        }
      }
    }
    // Nucleation through concentration fields.
    if (tracer_ && var.Int["enable_nucleation"]) {
      const auto& conf = tracer_->GetConf();
      for (auto l : GRange<size_t>(conf.layers)) {
        const auto sl = std::to_string(l);
        const Scal rate = var.Double("nucleation_rate" + sl, 0);
        if (rate) {
          const Scal cmax = var.Double["nucleation_cmax" + sl];
          const auto& tu = tracer_->GetVolumeFraction()[l];
          std::uniform_real_distribution<Scal> uniform(0, 1);
          const Scal k = rate * st_.dt;
          for (auto c : m.Cells()) {
            if ((tu[c] - cmax) * k > uniform(randgen_)) {
              nucl_cells_.insert(c);
            }
          }
        }
      }
    }
    if (tracer_ && particles_ && var.Int["enable_nucleation_points"]) {
      const auto& conf = tracer_->GetConf();
      std::vector<Vect> p_x;
      std::vector<bool> p_inner;
      std::vector<Vect> p_v;
      std::vector<Scal> p_id;
      std::vector<Scal> p_r;
      std::vector<Scal> p_source;
      std::vector<Scal> p_rho;
      std::vector<Scal> p_termvel;
      std::vector<Scal> p_removed;
      ParticlesView view{p_x,      p_inner, p_v,       p_id,     p_r,
                         p_source, p_rho,   p_termvel, p_removed};
      const Scal density = var.Double["particles_density"];
      for (auto l : GRange<size_t>(conf.layers)) {
        const auto sl = std::to_string(l);
        // Threshold for nucleation.
        const Scal cmax = var.Double["nucleation_cmax" + sl];
        // Post-nucleation concentration.
        const Scal cpost = var.Double["nucleation_cpost" + sl];
        const auto& tu = tracer_->GetVolumeFraction()[l];
        for (auto x : nucl_points_) {
          const auto c = m(m.GetCellFromPoint(x));
          if (tu[c] > cmax) {
            // Volume fraction is reduced to `cpost` and a particle with
            // the equivalent volume is created.
            const Scal du = tu[c] - cpost;
            const Scal vol = du * c.volume();
            fc_tracer_source_[l][c] -= du / st_.dt;
            view.x.push_back(x);
            view.inner.push_back(true);
            view.v.push_back(Vect(0));
            view.id.push_back(0);
            view.r.push_back(std::pow(vol * 3 / (4 * M_PI), 1. / 3));
            view.source.push_back(0);
            view.rho.push_back(density);
            // FIXME: Tracer id written to `termvel`, add a new particle field.
            view.termvel.push_back(l);
            view.removed.push_back(false);
          }
        }
      }
      particles_->Append(view);
    }
    if (tracer_ && electro_) {
      if (const Scal coeff = var.Double["flux_from_current"]) {
        auto& vmebc = tracer_->GetBCondMutable();
        auto& ff_current = electro_->GetFaceCurrent();
        me_group_.LoopPairs([&](auto cf_group) {
          auto cf = cf_group.first;
          size_t group = cf_group.second;
          auto nci = mebc_fluid_[cf].nci;
          const auto& custom = bc_group_custom_[group];
          auto getptr = [&](std::string key) -> const Scal* {
            auto it = custom.find(key);
            if (it != custom.end()) {
              return &it->second;
            }
            return nullptr;
          };

          const auto& conf = tracer_->GetConf();
          for (auto l : GRange<size_t>(conf.layers)) {
            auto sl = std::to_string(l);
            const Scal k = coeff * ff_current[cf];
            if (auto* ptr = getptr("tracer" + sl + "_dirichlet")) {
              vmebc[l][cf] = BCond<Scal>(BCondType::dirichlet, nci, (*ptr) * k);
            }
            if (auto* ptr = getptr("tracer" + sl + "_neumann")) {
              const auto d = conf.diffusion[l];
              vmebc[l][cf] = BCond<Scal>(
                  BCondType::neumann, nci, d > 0 ? (*ptr) * k / d : 0);
            }
          }
        });
      }
    }
    if (electro_ && var.Int["enable_electro_control"]) {
      mebc_electro_.LoopPairs([&](auto& p) {
        const auto cf = p.first;
        auto& bc = p.second;
        if (bc.type == BCondType::dirichlet) {
          if (const Scal* factor = electro_control_factor_.find(cf)) {
            bc.val = (*factor) * st_.electro_control_potential;
          }
        }
      });
      const Scal current = var.Double["electro_control_current"];
      const Scal rate = var.Double["electro_control_rate"];
      st_.electro_control_potential +=
          (current - electro_->GetStat().current) * rate;
    }
    // Compute total source of second component for statistics.
    if (tracer_) {
      st_.tmp_sum_src2 = 0;
      for (auto c : m.Cells()) {
        st_.tmp_sum_src2 += fc_src2_[c] * m.GetVolume(c);
      }
      m.Reduce(&st_.tmp_sum_src2, Reduction::sum);
    }
  }
  if (sem("reduce")) {
    st_.total_src2 += st_.tmp_sum_src2 * st_.dt;
  }
  if (particles_ && sem.Nested()) {
    particles_->SetUniqueIdForAppended(m);
  }
  // FIXME move, but keep inside nested call
  if (!vf_save_state_path_.empty() && sem.Nested("vf_save_state")) {
    if (auto as = dynamic_cast<ASVMEB*>(as_.get())) {
      as->SaveState(vf_save_state_path_);
    }
  }
}

template <class M>
void Hydro<M>::Dump(bool force) {
  auto sem = m.GetSem("dump");
  struct {
    std::vector<std::pair<std::string, std::vector<Scal>>> nucl_data;
  } * ctx(sem);
  auto& t = *ctx;
  if (sem.Nested("fields")) {
    if (dumper_.Try(st_.t, st_.dt) || force) {
      HydroPost<M>::DumpFields(this, m);
    }
  }
  if (dmptraj_.Try(st_.t, st_.dt) || force) {
    if (sem.Nested("trajdump")) {
      if (var.Int["enable_color"]) {
        auto plic = as_->GetPlic();
        if (eb_) {
          DumpTraj<EB>(
              *eb_, true, var, dmptraj_.GetN(), st_.t, layers, plic.vfcu,
              plic.vfccl_stat, plic.vfcim, fs_->GetPressure(),
              fs_->GetVelocity(), fcvm_, st_.dt);
        } else {
          DumpTraj<M>(
              m, true, var, dmptraj_.GetN(), st_.t, layers, plic.vfcu,
              plic.vfccl_stat, plic.vfcim, fs_->GetPressure(),
              fs_->GetVelocity(), fcvm_, st_.dt);
        }
      }
    }
  }
  if (sem("dmptrep")) {
    if (m.IsRoot() && dmptrep_.Try(st_.t, st_.dt)) {
      const std::string path = GetDumpName("trep", ".log", dmptrep_.GetN());
      m.TimerReport(path);
      if (!silent_) {
        std::cerr << std::fixed << std::setprecision(8) << "timer report"
                  << " t=" << st_.t << " to " << path << std::endl;
      }
    }
  }
  if (sem("dumpstat")) {
    if (m.IsRoot() && dumpstat_) {
      if (st_.step % var.Int("stat_step_every", 1) == 0 || force) {
        stat_->WriteValues(fstat_);
      }
    }
  }
  auto dump_part = [&](auto* as) {
    if (!as) {
      return;
    }
    if (auto* curv = dynamic_cast<const curvature::Particles<M>*>(
            curv_estimator_.get())) {
      if (dumper_.Try(st_.t, st_.dt)) {
        if (var.Int["dumppart"] && sem.Nested("part-dump")) {
          curv->GetParticles()->DumpParticles(
              as->GetAlpha(), as->GetNormal(), dumper_.GetN(), st_.t);
        }
        if (var.Int["dumppartinter"] && sem.Nested("partinter-dump")) {
          curv->GetParticles()->DumpPartInter(
              as->GetAlpha(), as->GetNormal(), dumper_.GetN(), st_.t);
        }
      }
    }
  };
  dump_part(dynamic_cast<ASV*>(as_.get()));
  dump_part(dynamic_cast<ASVEB*>(as_.get()));
  dump_part(dynamic_cast<ASVM*>(as_.get()));
  dump_part(dynamic_cast<ASVMEB*>(as_.get()));
  if (tracer_ && var.Int["enable_nucleation"] && var.Int("dump_nucl", 0) &&
      dumper_.Try(st_.t, st_.dt)) {
    const std::string path = GetDumpName("nucl", ".csv", dumper_.GetN(), -1);
    if (sem("dump-nucl-local")) {
      if (m.IsRoot()) {
        std::cerr << util::Format("dump t={:} to {}\n", st_.t, path);
      }
      t.nucl_data.push_back({"x", {}});
      t.nucl_data.push_back({"y", {}});
      t.nucl_data.push_back({"z", {}});
      for (auto c : nucl_cells_) {
        using Vect3 = generic::Vect<Scal, 3>;
        Vect3 x(m.GetCenter(c));
        t.nucl_data[0].second.push_back(x[0]);
        t.nucl_data[1].second.push_back(x[1]);
        t.nucl_data[2].second.push_back(x[2]);
      }
      for (auto xd : nucl_points_) {
        using Vect3 = generic::Vect<Scal, 3>;
        Vect3 x(xd);
        t.nucl_data[0].second.push_back(x[0]);
        t.nucl_data[1].second.push_back(x[1]);
        t.nucl_data[2].second.push_back(x[2]);
      }
    }
    if (sem.Nested()) {
      dump::DumpCsv(t.nucl_data, path, m);
    }
  }
}

template <class M>
void Hydro<M>::Run() {
  auto sem = m.GetSem("run");
  struct {
    Scal nabort;
  } * ctx(sem);

  if (sem.Nested("init") && !initialized_) {
    Init();
  }

  sem.LoopBegin();

  if (sem("events")) {
    if (events_) {
      vf_save_state_path_ = "";
      events_->Execute(st_.t);
    }
  }
  if (sem("loop-check")) {
    if (st_.t + st_.dt * 0.25 > var.Double["tmax"]) {
      if (m.IsRoot() && !silent_) {
        std::cerr << "End of simulation, t > tmax=" << var.Double["tmax"]
                  << std::endl;
      }
      sem.LoopBreak();
      finished_ = true;
    } else if (int(st_.step + 0.5) >= var.Int["max_step"]) {
      if (m.IsRoot() && !silent_) {
        std::cerr << "End of simulation, step > max_step="
                  << var.Int["max_step"] << std::endl;
      }
      sem.LoopBreak();
      finished_ = true;
    } else if (st_.step > 1 && fs_->GetError() < var.Double("stop_diff", 0)) {
      if (m.IsRoot() && !silent_) {
        std::cerr << "End of simulation, diff < stop_diff="
                  << var.Double["stop_diff"] << std::endl;
      }
      sem.LoopBreak();
      finished_ = true;
    } else {
      if (m.IsRoot() && !silent_) {
        if (st_.step % var.Int("report_step_every", 1) == 0 && !silent_) {
          ReportStep();
        }
      }
    }
  }
  CheckAbort(sem, ctx->nabort);

  if (sem("updatepar")) {
    if (auto fs = dynamic_cast<Simple<M>*>(fs_.get())) {
      auto par = ParsePar<Simple<M>>()(var);
      par.meshvel = st_.meshvel;
      fs->SetPar(par);
    }
    if (auto fs = dynamic_cast<Proj<M>*>(fs_.get())) {
      auto par = ParsePar<Proj<M>>()(var);
      par.meshvel = st_.meshvel;
      fs->SetPar(par);
    }
    if (linsolver_symm_) {
      linsolver_symm_->SetConf(linear::ModuleLinear<M>::GetConf(var, "symm"));
    }
    UpdateAdvectionPar();
    fcvm_ = fs_->GetVelocity();
  }
  if (sem.Nested("mixture")) {
    CalcMixture(as_->GetField());
  }
  if (sem.Nested("prestephook")) {
    PreStepHook(this);
  }
  if (sem.Nested("fs-start")) {
    fs_->StartStep();
  }
  if (sem.Nested("fs-iters")) {
    if (var.Int["enable_fluid"]) {
      StepFluid();
    }
  }
  if (sem.Nested("fs-finish")) {
    fs_->FinishStep();
  }
  if (sem.Nested("as-steps")) {
    if (var.Int["enable_advection"]) {
      StepAdvection();
    }
  }
  if (sem.Nested("tracer-step")) {
    if (tracer_) {
      StepTracer();
    }
  }
  if (sem.Nested("particles-step")) {
    if (particles_) {
      StepParticles();
    }
  }
  if (sem.Nested("electro-step")) {
    if (electro_) {
      StepElectro();
    }
  }
  if (sem.Nested("stat")) {
    CalcStat();
  }
  if (sem.Nested("stat")) {
    stat_->Update();
  }
  if (sem.Nested()) {
    CalcDt(); // must be after CalcStat to keep dt for moving mesh velocity
  }
  if (sem.Nested()) {
    Dump(false);
  }
  if (sem.Nested("stephook")) {
    StepHook(this);
  }
  if (par_.step_callback && sem.Nested("step_callback")) {
    par_.step_callback(par_.ptr, this);
  }
  if (sem.Nested("module_post_step") && module_post_step_) {
    (*module_post_step_)(this, m);
  }
  if (sem("inc")) {
    ++st_.step;
    if (var.Int("return_after_each_step", 0)) {
      sem.LoopBreak();
    }
  }
  sem.LoopEnd();

  if (sem.Nested("dumplast") && finished_) {
    if (var.Int["dumplast"]) {
      Dump(true);
    }
  }

  if (sem.Nested("posthook") && finished_) {
    if (eb_) {
      PostHook(var, fs_->GetVelocity(), m, *eb_);
    } else {
      PostHook(var, fs_->GetVelocity(), m);
    }
  }
  if (sem.Nested("finalhook") && finished_) {
    FinalHook(this);
  }
}

template <class M>
void Hydro<M>::ReportStep() {
  std::cerr << std::fixed << std::setprecision(8) << "STEP=" << st_.step
            << " t=" << st_.t << " dt=" << st_.dt
            << " wt=" << timer_.GetSeconds() << std::endl;
}

template <class M>
void Hydro<M>::ReportStepAdv() {
  std::cerr << std::fixed << std::setprecision(8)
            << ".....adv: t=" << as_->GetTime() << " dt=" << as_->GetTimeStep()
            << std::endl;
}

template <class M>
void Hydro<M>::ReportStepTracer() {
  std::cerr << std::fixed << std::setprecision(8)
            << ".....tracer: t=" << tracer_->GetTime() << " dt=" << tracer_dt_
            << std::endl;
}

template <class M>
void Hydro<M>::ReportStepParticles() {
  std::cerr << std::fixed << std::setprecision(8)
            << ".....particles: t=" << particles_->GetTime()
            << " dt=" << particles_dt_ << std::endl;
}

template <class M>
void Hydro<M>::ReportStepElectro() {
  std::cerr << std::fixed << std::setprecision(8)
            << ".....electro: t=" << electro_->GetTime()
            << " dt=" << fs_->GetTimeStep() << std::endl;
}

template <class M>
auto Hydro<M>::CalcPressureDrag(const FieldCell<Scal>& fcp, const Embed<M>& eb)
    -> Vect {
  MapEmbed<BCond<Scal>> me_pressure;
  mebc_fluid_.LoopBCond(eb, [&](auto cf, IdxCell, auto bc) { //
    const auto nci = bc.nci;
    if (bc.type == BCondFluidType::slipwall ||
        bc.type == BCondFluidType::symm) {
      me_pressure[cf] = BCond<Scal>(BCondType::neumann, nci);
    } else {
      me_pressure[cf] = BCond<Scal>(BCondType::extrap, nci);
    }
  });
  auto fep = UEB::Interpolate(fcp, me_pressure, eb);
  Vect sum(0);
  mebc_fluid_.LoopBCond(eb, [&, &m = m](auto cf, IdxCell c, auto bc) { //
    if (m.IsInner(c)) {
      if (bc.type == BCondFluidType::wall) {
        sum += eb.GetSurface(cf) * fep[cf];
      }
    }
  });
  return sum;
}

template <class M>
auto Hydro<M>::CalcViscousDrag(
    const FieldCell<Vect>& fcvel, const FieldCell<Scal>& fcmu,
    const Embed<M>& eb) -> Vect {
  MapEmbed<BCond<Scal>> me_neumann;
  mebc_fluid_.LoopBCond(eb, [&](auto cf, IdxCell, auto bc) { //
    me_neumann[cf] = BCond<Scal>(BCondType::neumann, bc.nci);
  });
  auto feg = UEB::Gradient(fcvel, fs_->GetVelocityCond(), eb);
  auto femu = UEB::Interpolate(fcmu, me_neumann, eb);
  Vect sum(0);
  mebc_fluid_.LoopBCond(eb, [&, &m = m](auto cf, IdxCell c, auto bc) { //
    if (m.IsInner(c)) {
      if (bc.type == BCondFluidType::wall) {
        sum += feg[cf] * (-eb.GetArea(cf) * femu[cf]);
      }
    }
  });
  return sum;
}

template <class M>
void Hydro<M>::ReportIter() {
  std::cerr << std::scientific << std::setprecision(16)
            << ".....iter=" << fs_->GetIter() << ", diff=" << fs_->GetError()
            << std::endl;
}

template <class M>
void Hydro<M>::CheckAbort(Sem& sem, Scal& nabort) {
  if (sem("abort-local")) {
    const Scal abortvel = var.Double["abortvel"];
    CHECKNAN(as_->GetField(), true)
    CHECKNAN(fs_->GetVelocity(), true)
    CHECKNAN(fs_->GetPressure(), true)
    for (auto c : m.Cells()) {
      if (fs_->GetVelocity()[c].norminf() > abortvel) {
        std::cerr << util::Format(
            "abortvel exceeded at x={}\n", m.GetCenter(c));
        nabort += 1;
        break;
      }
    }
    m.Reduce(&nabort, Reduction::sum);
  }

  if (sem("abort-reduce")) {
    if (nabort != 0) {
      if (m.IsRoot()) {
        std::cerr << "nabort = " << nabort << std::endl;
      }
      sem.LoopBreak();
      finished_ = true;
    }
  }
}

template <class M>
void Hydro<M>::StepFluid() {
  auto sem = m.GetSem("iter"); // sem nested
  if (sem("iter")) {
    OverwriteBc();
  }
  sem.LoopBegin();
  if (sem.Nested("iter")) {
    fs_->MakeIteration();
  }
  if (sem("report")) {
    ++st_.iter;
    if (m.IsLead()) {
      this->var_mutable.Int["iter"] = st_.iter;
    }
    if (m.IsRoot()) {
      if (st_.step % var.Int("report_step_every", 1) == 0 && !silent_) {
        ReportIter();
      }
    }
  }
  if (sem("convcheck")) {
    auto it = fs_->GetIter();
    if ((fs_->GetError() < var.Double["tol"] &&
         (int)it >= var.Int["min_iter"]) ||
        (int)it >= var.Int["max_iter"]) {
      sem.LoopBreak();
    }
  }
  // TODO: Suspender loop hangs if (probably) Nested is last
  sem.LoopEnd();
}

template <class M>
void Hydro<M>::StepTracer() {
  auto sem = m.GetSem("tracer-steps"); // sem nested
  sem.LoopBegin();
  if (sem("spawn")) {
    SpawnTracer();
  }
  if (sem.Nested("start")) {
    tracer_->Step(tracer_dt_, fs_->GetVolumeFlux());
  }
  if (sem("report")) {
    if (m.IsRoot()) {
      if (st_.step % var.Int("report_step_every", 1) == 0 && !silent_) {
        ReportStepTracer();
      }
    }
  }
  if (sem("convcheck")) {
    if (tracer_->GetTime() >= fs_->GetTime() - 0.5 * tracer_dt_) {
      sem.LoopBreak();
    }
  }
  sem.LoopEnd();
}

template <class M>
void Hydro<M>::StepParticles() {
  auto sem = m.GetSem("particles-steps"); // sem nested
  sem.LoopBegin();
  if (sem.Nested()) {
    auto velocity_hook = [this](const ParticlesView& view) {
      // Move particles away from the wall.
      for (size_t i = 0; i < view.x.size(); ++i) {
        if (view.inner[i]) {
          const IdxCell c = m.GetCellFromPoint(view.x[i]);
          const Vect xc = m.GetCenter(c);
          const Vect walldist = fc_wall_dist_[c];
          if (IsFinite(walldist) && walldist.sqrnorm() > 0) {
            const Vect walldir = walldist / walldist.norm();
            // Distance from particle center to wall along `walldist`
            const Scal dz = walldir.dot(walldist + xc - view.x[i]);
            if (dz - view.r[i] < 0) {
              // If particle intersects the wall, move it away
              view.v[i] += walldir * (dz - view.r[i]) / particles_dt_;
            }
          }
        }
      }
      if (particles_hook_) {
        particles_hook_(view);
      }
    };
    particles_->Step(
        particles_dt_, fs_->GetVolumeFlux(), fs_->GetVelocityCond(),
        velocity_hook);
  }
  if (var.Int["particles_to_vof"] && sem.Nested("partilces_to_vof")) {
    ConvertParticlesToVof();
  }
  // Coalescence of particles.
  if (var.Int["particles_coal"] && sem("particles_coal")) {
    auto view = particles_->GetView();
    for (size_t i = 0; i < view.x.size(); ++i) {
      if (view.inner[i]) {
        for (size_t j = 0; j < view.x.size(); ++j) {
          const Vect xi = view.x[i];
          const Vect xj = view.x[j];
          if (i == j) {
            continue;
          }
          const Scal ri = view.r[i];
          const Scal rj = view.r[j];
          if (xi.dist(xj) < view.r[i] + view.r[j]) {
            // Remove the lexicographically minimal particle
            // and enlarge the other one. Modify only inner particles.
            const bool less = ri < rj || (ri == rj && xi.lexless(xj));
            if (less) {
              view.removed[i] = true;
            } else  {
              const Scal voli = std::pow(ri, 3);
              const Scal volj = std::pow(rj, 3);
              view.r[i] = std::pow(voli + volj, 1. / 3);
              view.x[i] = (voli * xi + volj * xj) / (voli + volj);
              view.v[i] = (voli * view.v[i] + volj * view.v[j]) / (voli + volj);
            }
          }
        }
      }
    }
  }
  if (sem("spawn")) {
    std::vector<Vect> p_x;
    std::vector<bool> p_inner;
    std::vector<Vect> p_v;
    std::vector<Scal> p_id;
    std::vector<Scal> p_r;
    std::vector<Scal> p_source;
    std::vector<Scal> p_rho;
    std::vector<Scal> p_termvel;
    std::vector<Scal> p_removed;
    ParticlesView view{p_x,      p_inner, p_v,       p_id,     p_r,
                       p_source, p_rho,   p_termvel, p_removed};
    SpawnParticles(view);
    particles_->Append(view);
  }
  if (sem("report")) {
    if (m.IsRoot()) {
      if (st_.step % var.Int("report_step_every", 1) == 0 && !silent_) {
        ReportStepParticles();
      }
    }
  }
  if (sem("convcheck")) {
    if (particles_->GetTime() >= fs_->GetTime() - 0.5 * particles_dt_) {
      sem.LoopBreak();
    }
  }
  sem.LoopEnd();
}

template <class M>
void Hydro<M>::ConvertParticlesToVof() {
  auto sem = m.GetSem();
  struct Particle {
    Vect x;
    Scal r;
  };

  struct {
    Scal nconv = 0;
    std::vector<Scal> intersects_vof; // 1 if particle intersects a cell
                                      // with large volume fraction
  } * ctx(sem);
  auto& t = *ctx;

  if (sem()) {
    const auto& fcu = as_->GetField();
    auto view = particles_->GetView();
    const auto h = m.GetCellSize();
    t.intersects_vof.resize(view.x.size(), 0);
    const auto& blockall = m.GetAllBlockCells();
    const MIdx blockall_wmin = blockall.GetBegin();
    const MIdx blockall_wmax = blockall.GetEnd() - MIdx(1);
    for (size_t i = 0; i < view.x.size(); ++i) {
      if (view.inner[i]) {
        // TODO: Avoid code duplication with growth of particles from tracer.
        // Range of cells containing the particle.
        MIdx wmin(((view.x[i] - Vect(view.r[i]) - m.flags.global_origin) / h)
                      .floor());
        MIdx wmax(((view.x[i] + Vect(view.r[i]) - m.flags.global_origin) / h)
                      .floor());
        // Cell containing the particle centroid.
        const auto cpart = m(m.GetCellFromPoint(view.x[i]));
        // Limit the range of cells to 5x5x5 stencil around the centroid.
        const MIdx wcpart(cpart);
        wmin = wmin.max(wcpart - MIdx(2));
        wmax = wmax.min(wcpart + MIdx(2));
        // Limit the range of cells to inner and halo cells of current block.
        wmin = wmin.max(blockall_wmin);
        wmax = wmax.min(blockall_wmax);
        // Block of cells from range [wmin, wmax].
        const GBlock<int, M::dim> wcells(wmin, wmax - wmin + MIdx(1));
        const auto& index = m.GetIndexCells();
        for (MIdx w : wcells) {
          const auto c = m(index.GetIdx(w));
          const Vect xc = c.center;
          if (fcu[c] > 0.5 &&
              xc.sqrdist(view.x[i]) < sqr(view.r[i] + h[0] * 0.5)) {
            t.intersects_vof[i] = 1;
          }
        }
      }
    }

    auto req = GetCommPartRequest<M>(view);
    req.attr_scal.push_back(&t.intersects_vof);
    m.CommPart(req);
  }

  if (sem()) {
    auto view = particles_->GetView();
    const Scal r0 = var.Double["particles_to_vof_radius"] * m.GetCellSize()[0];
    const Scal dt = as_->GetTimeStep();

    std::vector<Particle> toconvert; // Particles to convert to vof.
    for (size_t i = 0; i < view.x.size(); ++i) {
      // Convert if particle is sufficiently large
      // or intersects a cell with large volume fraction.
      if (view.r[i] > r0 || t.intersects_vof[i]) {
        toconvert.push_back({view.x[i] + view.v[i] * dt, view.r[i]});
        view.removed[i] = true;
      }
    }
    t.nconv += toconvert.size();
    m.Reduce(&t.nconv, Reduction::sum);

    auto mod = [this, toconvert](
                   FieldCell<Scal>& fcu, FieldCell<Scal>& fccl, auto&) {
      (void)fccl;
      auto ls = [&](const Vect& x) -> Scal {
        Scal res = -std::numeric_limits<Scal>::max();
        for (auto& p : toconvert) {
          res = std::max(res, p.r - x.dist(p.x));
        }
        return res;
      };
      for (auto c : m.AllCellsM()) {
        fcu[c] = std::max(
            fcu[c],
            GetLevelSetVolume<Scal, M::dim>(ls, c.center, m.GetCellSize()));
      }
    };
    if (auto* as = dynamic_cast<Vof<M>*>(as_.get())) {
      as->AddModifier(mod);
    }
  }
  if (sem()) {
    if (m.IsRoot() && var.Int("particles_to_vof_verbose", 0) && t.nconv) {
      std::cout << util::Format("Converted {:} particles", t.nconv)
                << std::endl;
    }
  }
}

template <class M>
void Hydro<M>::StepElectro() {
  auto sem = m.GetSem(__func__);
  if (sem.Nested("start")) {
    electro_->Step(fs_->GetTimeStep(), as_->GetField());
  }
  if (sem("report")) {
    if (m.IsRoot()) {
      if (st_.step % var.Int("report_step_every", 1) == 0 && !silent_) {
        ReportStepElectro();
      }
    }
  }
}

template <class M>
void Hydro<M>::StepAdvection() {
  auto sem = m.GetSem("steps"); // sem nested
  sem.LoopBegin();
  if (auto as = dynamic_cast<ASVM*>(as_.get())) {
    const Scal* const voidpenal = var.Double.Find("voidpenal");
    if (voidpenal && sem("void-penal")) {
      auto fccl = as->GetColor();
      auto fcu = as->GetFieldM();
      for (auto f : m.Faces()) {
        const IdxCell cm = m.GetCell(f, 0);
        const IdxCell cp = m.GetCell(f, 1);
        Scal um = 0;
        Scal up = 0;
        for (auto l : layers) {
          if ((*fccl[l])[cm] != kClNone) {
            um += (*fcu[l])[cm];
          }
          if ((*fccl[l])[cp] != kClNone) {
            up += (*fcu[l])[cp];
          }
        }
        um = std::min(1., um);
        up = std::min(1., up);
        FieldFace<Scal>& ffv =
            const_cast<FieldFace<Scal>&>(fs_->GetVolumeFlux().GetFieldFace());
        ffv[f] += -(up - um) * (*voidpenal) * m.GetArea(f);
      }
    }
  }
  if (sem.Nested("start")) {
    as_->StartStep();
  }
  if (sem.Nested("iter")) {
    as_->MakeIteration();
  }
  if (sem.Nested("finish")) {
    as_->FinishStep();
  }
  if (sem("report")) {
    if (m.IsRoot()) {
      if (st_.step % var.Int("report_step_every", 1) == 0 && !silent_) {
        ReportStepAdv();
      }
    }
  }
  if (sem("convcheck")) {
    if (as_->GetTime() >= fs_->GetTime() - 0.5 * as_->GetTimeStep()) {
      sem.LoopBreak();
    }
  }
  sem.LoopEnd();
  if (sem.Nested("as-post")) {
    as_->PostStep();
  }
  if (sem.Nested("curv")) {
    if (eb_) {
      curv_estimator_->CalcCurvature(fck_, as_->GetPlic(), m, *eb_);
    } else {
      curv_estimator_->CalcCurvature(fck_, as_->GetPlic(), m, m);
    }
  }
  if (var.Int["enable_bubgen"]) {
    if (sem.Nested("bubgen")) {
      StepBubgen();
    }
  }
  if (var.Int["enable_erasevf"]) {
    if (sem("erasevf")) {
      StepEraseVolumeFraction("erasevf", erasevf_last_t_);
      StepEraseVolumeFraction("erasevf2", erasevf2_last_t_);
    }
  }
  if (sem.Nested("erasecl")) {
    StepEraseColor("erasecl");
  }
}

template <class M>
void Hydro<M>::InitTracerFields(Multi<FieldCell<Scal>>& vfcu) {
  auto sem = m.GetSem("tracerfields");
  struct {
    Vars vart;
  } * ctx(sem);
  auto& t = *ctx;
  if (sem("init")) {
    vfcu.resize(var.Int["tracer_layers"]);
    vfcu.InitAll(FieldCell<Scal>(m, 0));
  }
  for (int l = 0; l < var.Int["tracer_layers"]; ++l) {
    const std::string prefix = "tracer" + std::to_string(l);
    if (sem("var" + std::to_string(l))) {
      t.vart.String.Set("init_vf", var.String[prefix + "_init"]);
      t.vart.Double.Set("init_vf_value", var.Double(prefix + "_init_value", 0));
      t.vart.String.Set("list_path", var.String[prefix + "_list_path"]);
      t.vart.Int.Set("dim", var.Int["dim"]);
      t.vart.Int.Set("list_ls", 3);
    }
    if (sem.Nested("field" + std::to_string(l))) {
      InitVf(vfcu[l], t.vart, m, !silent_);
    }
    if (sem("factor" + std::to_string(l))) {
      auto k = var.Double[prefix + "_factor"];
      for (auto c : m.AllCells()) {
        vfcu[l][c] *= k;
      }
    }
  }
  if (sem()) {
  }
}

template <class M>
void Hydro<M>::StepBubgen() {
  auto sem = m.GetSem("bubgen");
  struct {
    std::shared_ptr<FieldCell<Scal>> fcvf; // volume fraction
    Vars var;
    bool verbose = false;
  } * ctx(sem);
  auto& t = *ctx;
  if (bubgen_.Try(st_.t, st_.dt)) {
    if (sem("as-bubgen-var")) {
      t.var.String.Set("init_vf", "list");
      t.var.String.Set("list_path", var.String["bubgen_path"]);
      t.var.Int.Set("dim", var.Int["dim"]);
      t.var.Int.Set("list_ls", var.Int["list_ls"]);
      t.fcvf = std::make_shared<FieldCell<Scal>>(m);
      t.verbose = var.Int("verbose_bubgen", 0);
    }
    if (sem.Nested("as-bubgen-initvf")) {
      InitVf(*t.fcvf, t.var, m, t.verbose);
    }
    if (sem("as-bubgen-apply")) {
      const Scal clnew = fs_->GetTime();
      auto modify = [fcvf = t.fcvf, clnew](auto& u, auto& cl, auto& eb) {
        for (auto c : eb.AllCells()) {
          if ((*fcvf)[c] > 0) {
            const Scal v = std::min((*fcvf)[c], eb.GetVolumeFraction(c));
            u[c] = std::max(u[c], v);
            cl[c] = clnew;
          }
        }
      };
      auto modifym = [fcvf = t.fcvf, clnew](auto& u, auto& cl, auto, auto& eb) {
        for (auto c : eb.AllCells()) {
          if ((*fcvf)[c] > 0) {
            const Scal v = std::min((*fcvf)[c], eb.GetVolumeFraction(c));
            (*u[0])[c] = std::max((*u[0])[c], v);
            (*cl[0])[c] = clnew;
          }
        }
      };
      if (auto* as = dynamic_cast<ASVEB*>(as_.get())) {
        as->AddModifier(modify);
      }
      if (auto* as = dynamic_cast<ASVMEB*>(as_.get())) {
        as->AddModifier(modifym);
      }
      if (auto* as = dynamic_cast<ASV*>(as_.get())) {
        as->AddModifier(modify);
      }
      if (auto* as = dynamic_cast<ASVM*>(as_.get())) {
        as->AddModifier(modifym);
      }
    }
    if (sem()) {
      // FIXME: empty stage to finish communication to keep ctx
    }
  }
}

template <class M>
void Hydro<M>::StepEraseVolumeFraction(std::string prefix, Scal& last_t) {
  if (!var.Double.Contains(prefix + "_t0")) {
    return;
  }
  const Vect rect_x0(var.Vect[prefix + "_rect_x0"]);
  const Vect rect_x1(var.Vect[prefix + "_rect_x1"]);
  const Rect<Vect> rect(rect_x0, rect_x1);
  const Scal t0 = var.Double[prefix + "_t0"];
  const Scal tper = var.Double[prefix + "_per"];
  if (st_.t > t0 && st_.t - last_t >= tper) {
    if (m.IsRoot() && tper > fs_->GetTimeStep()) {
      std::cerr << prefix + " t=" << st_.t << std::endl;
    }
    last_t = st_.t;
    auto apply_vof = [this, &rect](auto* as, const auto& eb) {
      if (as) {
        auto& u = const_cast<FieldCell<Scal>&>(as->GetField());
        for (auto c : eb.AllCells()) {
          const auto x = m.GetCenter(c);
          if (rect.IsInside(x) && u[c] > 0) {
            u[c] = 0;
          }
        }
      }
    };
    auto apply_vofm = [this, &rect](auto* as, const auto& eb) {
      if (as) {
        for (auto l : layers) {
          auto& u = const_cast<FieldCell<Scal>&>(*as->GetFieldM()[l]);
          auto& cl = const_cast<FieldCell<Scal>&>(*as->GetColor()[l]);
          for (auto c : eb.AllCells()) {
            const auto x = m.GetCenter(c);
            if (rect.IsInside(x) && u[c] > 0) {
              u[c] = 0;
              cl[c] = kClNone;
            }
          }
        }
      }
    };
    if (eb_) {
      apply_vofm(dynamic_cast<ASVMEB*>(as_.get()), *eb_);
      apply_vof(dynamic_cast<ASVEB*>(as_.get()), *eb_);
    } else {
      apply_vofm(dynamic_cast<ASVM*>(as_.get()), m);
      apply_vof(dynamic_cast<ASV*>(as_.get()), m);
    }
  }
}
template <class M>
void Hydro<M>::ReportSysinfo(std::ostream& out) {
  auto sem = m.GetSem(__func__);
  struct {
    sysinfo::Info info;
    std::vector<std::vector<char>> hostname;
    std::vector<size_t> cuda_mem;
    std::vector<uint16_t> cuda_uuid;
  } * ctx(sem);
  auto& t = *ctx;
  if (var.Int("report_sysinfo", 0)) {
    if (sem()) {
      if (m.IsLead()) {
        t.info = sysinfo::GetInfo(sysinfo::InfoSelect());
        t.hostname.push_back({t.info.hostname.begin(), t.info.hostname.end()});
        if (t.info.cuda_enabled) {
          t.cuda_mem.push_back(t.info.cuda_mem);
          t.cuda_uuid.push_back(t.info.cuda_uuid);
        }
      }
      m.Reduce(&t.hostname, Reduction::concat);
      m.Reduce(&t.cuda_mem, Reduction::concat);
      m.Reduce(&t.cuda_uuid, Reduction::concat);
    }
    if (sem()) {
      if (m.IsRoot()) {
        out << "\nOpenMP num_threads: " << t.info.omp_num_threads;
        out << "\nOpenMP max_threads: " << t.info.omp_max_threads;
        {
          std::map<std::string, int> map;
          for (auto h : t.hostname) {
            ++map[std::string(h.begin(), h.end())];
          }
          out << "\nHostname\n";
          for (auto p : map) {
            out << p.first << ' ' << p.second << '\n';
          }
        }
        if (!t.cuda_uuid.empty()) {
          out << "\nCUDA GPU UUID\n";
          for (auto u : t.cuda_uuid) {
            out << util::Format("{:04X} ", u);
          }
          out << "\nCUDA GPU MEM (GiB)\n";
          for (auto u : t.cuda_mem) {
            out << util::Format("{:.3f} ", double(u) / (1 << 30));
          }
          out << '\n';
        }
        out << '\n';
      }
    }
  }
}

// Finds the minimal color in rectangle [x0,x1]
// and clears the volume fraction in all cells with that color.
template <class M>
void Hydro<M>::StepEraseColor(std::string prefix) {
  if (!var.Vect.Contains(prefix + "_rect_x0")) {
    return;
  }
  const Vect rect_x0(var.Vect[prefix + "_rect_x0"]);
  const Vect rect_x1(var.Vect[prefix + "_rect_x1"]);
  const Rect<Vect> rect(rect_x0, rect_x1);

  auto sem = m.GetSem();
  struct {
    Scal cl = std::numeric_limits<Scal>::max();
  } * ctx(sem);
  auto& t = *ctx;
  auto apply = [this, &rect, &t, &sem](auto* as, const auto& eb) {
    if (as) {
      if (sem()) {
        auto plic = as->GetPlic();
        for (auto l : layers) {
          const auto& u = *plic.vfcu[l];
          const auto& cl = *plic.vfccl_stat[l];
          for (auto c : eb.AllCells()) {
            const auto x = m.GetCenter(c);
            if (rect.IsInside(x) && u[c] > 0 && cl[c] != kClNone) {
              t.cl = cl[c];
            }
          }
        }
        m.Reduce(&t.cl, Reduction::min);
      }
      if (sem()) {
        auto plic = as->GetPlic();
        for (auto l : layers) {
          auto& u = const_cast<FieldCell<Scal>&>(*plic.vfcu[l]);
          auto& cl = const_cast<FieldCell<Scal>&>(*plic.vfccl_stat[l]);
          for (auto c : eb.AllCells()) {
            if (cl[c] == t.cl) {
              u[c] = 0;
              cl[c] = kClNone;
            }
          }
        }
      }
    }
  };
  if (eb_) {
    apply(as_.get(), *eb_);
  } else {
    apply(as_.get(), m);
  }
}
