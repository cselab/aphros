// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <mpi.h>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdio>
#include <ctime>
#include <fstream>
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
#include "dump/dumper.h"
#include "dump/output.h"
#include "func/init.h"
#include "func/init_bc.h"
#include "func/init_contang.h"
#include "geom/mesh.h"
#include "kernelmeshpar.h"
#include "parse/curv.h"
#include "parse/parser.h"
#include "parse/proj.h"
#include "parse/simple.h"
#include "parse/tvd.h"
#include "parse/util.h"
#include "parse/vars.h"
#include "parse/vof.h"
#include "solver/advection.h"
#include "solver/approx.h"
#include "solver/approx_eb.h"
#include "solver/curv.h"
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
#include "solver/tvd.h"
#include "solver/vof.h"
#include "solver/vofm.h"
#include "util/convdiff.h"
#include "util/events.h"
#include "util/hydro.h"
#include "util/metrics.h"
#include "util/posthook.h"
#include "young/young.h"

class GPar {};

template <class M>
FieldCell<typename M::Scal> GetDivergence(
    const FieldFace<typename M::Scal>& ffv, const M& m, const Embed<M>& eb) {
  using Scal = typename M::Scal;
  FieldCell<Scal> fcdiv(m, 0);
  for (auto c : eb.Cells()) {
    Scal div = 0;
    for (auto q : eb.Nci(c)) {
      div += ffv[m.GetFace(c, q)] * m.GetOutwardFactor(c, q);
    }
    fcdiv[c] = div / eb.GetVolume(c);
  }
  return fcdiv;
}

template <class M>
FieldCell<typename M::Vect> GetVort(
    const FieldCell<typename M::Vect>& fcvel,
    const MapEmbed<BCond<typename M::Vect>>& me_vel, Embed<M>& eb) {
  auto& m = eb.GetMesh();
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using UEB = UEmbed<M>;

  std::array<FieldCell<Vect>, 3> grad;
  for (size_t d = 0; d < 3; ++d) {
    grad[d].Reinit(eb, Vect(0));
    const auto mebc = GetScalarCond(me_vel, d, m);
    const FieldCell<Scal> fcu = GetComponent(fcvel, d);
    const FieldEmbed<Scal> ffg = UEB::Gradient(fcu, mebc, eb);
    grad[d] = UEB::AverageGradient(ffg, eb);
  }

  FieldCell<Vect> r(m, Vect(0));
  for (auto c : eb.Cells()) {
    r[c][0] = grad[2][c][1] - grad[1][c][2];
    r[c][1] = grad[0][c][2] - grad[2][c][0];
    r[c][2] = grad[1][c][0] - grad[0][c][1];
  }

  return r;
}

template <class M_>
class Hydro : public KernelMeshPar<M_, GPar> {
 public:
  using P = KernelMeshPar<M_, GPar>; // parent
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;
  using Sem = typename M::Sem;
  using Par = GPar;
  template <class T>
  using Multi = Multi<T>;
  using UEB = UEmbed<M>;
  using ParticlesView = typename ParticlesInterface<M>::ParticlesView;
  static constexpr size_t dim = M::dim;
  friend void StepHook<>(Hydro*);

  // TODO: issue warning if variable in Vars was not used
  // but differs from default (like in CMake)
  Hydro(Vars&, const MyBlockInfo&, Par&);
  void Run() override;
  M& GetMesh() {
    return m;
  }

 protected:
  using P::bi_;
  using P::m;
  using P::var;

 private:
  void Init();
  void InitEmbed();
  void InitTracer(Multi<FieldCell<Scal>>& vfcu);
  void InitTracerFields(Multi<FieldCell<Scal>>& vfcu);
  void SpawnTracer();
  void InitParticles();
  void SpawnParticles(ParticlesView& view);
  void OverwriteBc();
  void InitFluid(const FieldCell<Vect>& fc_vel);
  void InitAdvection(const FieldCell<Scal>& fcvf, const FieldCell<Scal>& fccl);
  void InitStat();
  Vect CalcPressureDrag(const FieldCell<Scal>& fcp, const Embed<M>& eb);
  Vect CalcViscousDrag(
      const FieldCell<Vect>& fcvel, const FieldCell<Scal>& fcmu,
      const Embed<M>& eb);
  void DumpFields();
  void Dump(bool force);
  // Calc rho, mu and force based on volume fraction
  void CalcMixture(const FieldCell<Scal>& vf);
  // Clips v to given range, uses const_cast
  void Clip(const FieldCell<Scal>& v, Scal min, Scal max);
  void CalcStat();
  void CalcDt();
  void ReportStep();
  void ReportStepAdv();
  void ReportStepTracer();
  void ReportStepParticles();
  void ReportIter();
  // Issue sem.LoopBreak if abort conditions met
  void CheckAbort(Sem& sem, Scal& nabort);
  void StepFluid();
  void StepAdvection();
  void StepTracer();
  void StepParticles();
  void StepBubgen();

  using AST = Tvd<M>; // advection TVD
  using ASV = Vof<M>; // advection VOF
  using ASVM = Vofm<M>; // advection multi VOF
  using ASVEB = Vof<Embed<M>>; // advection VOF embed
  using ASVMEB = Vofm<Embed<M>>; // advection multi VOF embed
  using EB = Embed<M>;
  static constexpr Scal kClNone = ASVM::kClNone;

  void UpdateAdvectionPar() {
    if (auto as = dynamic_cast<AST*>(as_.get())) {
      as->SetPar(ParsePar<AST>()(var));
    } else if (auto as = dynamic_cast<ASV*>(as_.get())) {
      as->SetPar(ParsePar<ASV>()(var));
    } else if (auto as = dynamic_cast<ASVM*>(as_.get())) {
      as->SetPar(ParsePar<ASVM>()(var));
    } else if (auto as = dynamic_cast<ASVEB*>(as_.get())) {
      as->SetPar(ParsePar<ASVEB>()(var));
    } else if (auto as = dynamic_cast<ASVMEB*>(as_.get())) {
      as->SetPar(ParsePar<ASVMEB>()(var));
    }
  }
  // Surface tension time step
  Scal GetStDt() {
    const Scal sig = var.Double["sigma"];
    const Scal* cflst = var.Double.Find("cflst");
    if (cflst && sig != 0.) {
      Scal pi = M_PI;
      Scal h3 = m.GetVolume(IdxCell(0));
      Scal r1 = var.Double["rho1"];
      Scal r2 = var.Double["rho2"];
      return (*cflst) * std::sqrt(h3 * (r1 + r2) / (4. * pi * std::abs(sig)));
    }
    return std::numeric_limits<Scal>::max();
  }
  // Viscosity time step
  Scal GetVisDt() {
    const Scal rho1 = var.Double["rho1"];
    const Scal rho2 = var.Double["rho2"];
    const Scal mu1 = var.Double["mu1"];
    const Scal mu2 = var.Double["mu2"];
    const Scal nu1 = mu1 / rho1;
    const Scal nu2 = mu2 / rho2;
    const Scal num = std::max(nu1, nu2);
    const Scal* cflvis = var.Double.Find("cflvis");
    if (cflvis && num != 0.) {
      const Scal h2 = sqr(m.GetCellSize()[0]); // XXX adhoc cubic cell
      return (*cflvis) * h2 / num;
    }
    return std::numeric_limits<Scal>::max();
  }
  YoungParam GetYoungPar() const {
    YoungParam q;
    q.rhov = var.Double["rho1"];
    q.rhou = var.Double["rho2"];
    q.muv = var.Double["mu1"];
    q.muu = var.Double["mu2"];
    q.hv = 1.;
    q.hu = 1.;
    q.gamma0 = var.Double["sigma"];
    q.gamma1 = var.Vect["sig_grad"][0];
    q.T0 = 0.;
    q.T1 = 1.;
    q.R = var.Double["youngbc_r"];
    return q;
  }
  void InitYoung() {
    young_ini(GetYoungPar());
  }
  void CalcVort() {
    auto& fcv = fs_->GetVelocity();
    if (eb_) {
      fcom_ = GetVort(fcv, fs_->GetVelocityCond(), *eb_);
    } else {
      fcom_ = GetVort(fcv, GetCond<Vect>(fs_->GetVelocityCond()), m);
    }
    fcomm_.Reinit(m);
    for (auto c : m.Cells()) {
      fcomm_[c] = fcom_[c].norm();
    }
  }
  void CalcStrain(const FieldCell<Vect> fcvel, FieldCell<Scal>& fc_strain) {
    auto& fcv = fcvel;
    auto& fcs = fc_strain;

    auto ffv = UEmbed<M>::Interpolate(fcv, fs_->GetVelocityCond(), m);

    std::array<FieldCell<Vect>, dim> g; // g[i][c][j] is derivative du_i/dx_j
    for (size_t i = 0; i < dim; ++i) {
      g[i] = Gradient(GetComponent(ffv, i), m);
    }

    fcs.Reinit(m, 0);
    int edim = var.Int["dim"];
    for (auto c : m.Cells()) {
      for (int i = 0; i < edim; ++i) {
        for (int j = 0; j < edim; ++j) {
          fcs[c] += sqr(g[i][c][j]) + g[i][c][j] * g[j][c][i];
        }
      }
      fcs[c] *= 0.5;
    }
  }
  FieldCell<Scal> GetDiv() {
    auto& ffv = fs_->GetVolumeFlux().GetFieldFace();
    if (eb_) {
      return GetDivergence(ffv, m, *eb_);
    }
    FieldCell<Scal> fc(m, 0); // result
    for (auto c : m.Cells()) {
      for (auto q : m.Nci(c)) {
        IdxFace f = m.GetFace(c, q);
        fc[c] += ffv[f] * m.GetOutwardFactor(c, q);
      }
      fc[c] /= m.GetVolume(c);
    }
    return fc;
  }
  Vect GetYoungVel(Vect x) const {
    x -= Vect(0.5);
    // 0: streamwise
    // 1,2: crossstream with symmetry around 0-axis
    Scal r = Vect(0, x[1], x[2]).norm();
    Vect v(0);
    Scal vc; // vel cross
    young_fields(r, x[0], &vc, &v[0], nullptr, nullptr);
    Scal e1 = x[1] / r;
    Scal e2 = x[2] / r;
    v[1] = vc * e1;
    v[2] = vc * e2;
    return v;
  }

  GRange<size_t> layers;
  FieldCell<Scal> fc_mu_; // viscosity
  FieldCell<Scal> fc_rho_; // density
  FieldCell<Scal> fc_src_; // source of mixture volume
  FieldCell<Scal> fc_src2_; // source of second phase volume
  FieldCell<Scal> fc_srcm_; // mass source
  FieldCell<Vect> fc_force_; // force
  FieldCell<Scal> fc_dist_; // distance from eb
  FieldCell<Scal> fc_phi_; // distance from eb
  FieldEmbed<Scal> febp_; // balanced force projections

  MapEmbed<BCondAdvection<Scal>> mf_adv_;
  MapCondFace mf_cond_vfsm_;
  MapEmbed<BCondFluid<Vect>> mebc_fluid_;
  MapEmbed<BCondFluid<Vect>> mebc_fluid_orig_;
  MapEmbed<size_t> me_group_;
  MapEmbed<Scal> me_contang_;
  std::vector<std::string> vdesc_;
  MapCondFaceFluid mf_fluid_;
  MapCell<std::shared_ptr<CondCell>> mc_cond_;
  MapCell<std::shared_ptr<CondCellFluid>> mc_velcond_;

  std::unique_ptr<Embed<M>> eb_;
  std::unique_ptr<AdvectionSolver<M>> as_;
  std::unique_ptr<FluidSolver<M>> fs_;

  FieldCell<Scal> fc_smvf_; // smoothed volume fraction used by CalcMixture()
  FieldFace<Scal> ff_smvf_; // interpolated fc_smvf_

  FieldCell<Scal> fc_sig_; // surface tension sigma
  FieldCell<Scal> fc_contang_; // contact angle

  FieldCell<Vect> fcvm_; // velocity field time_prev // TODO: revise

  FieldCell<Vect> fcyv_; // Young velocity
  FieldCell<Vect> fcom_; // vorticity
  FieldCell<Scal> fcomm_; // vorticity magnitude
  FieldCell<Scal> fcbc_; // boundary condition type by GetBcField()
  FieldCell<Scal> fc_strain_; // double inner product of strain rate tensor
  Multi<FieldCell<Scal>> fck_; // curvature
  typename PartStrMeshM<M>::Par psm_par_;
  std::unique_ptr<PartStrMeshM<M>> psm_;

  Scal bgt_ = -1.; // bubgen last time

  struct Stat {
    struct Vofm {
      Scal cells_vf; // number of cells with vf>0
      Scal cells_cl; // number of cells wih cl != kClNone
      Scal sum_vf; // sum of vf
      Scal hist; // vofm[l].hist is number of cells containing l+1
                 // cells with vf>0
    };

    Scal m1, m2; // volume
    Scal m20; // initial volume
    Scal m2d; // relative volume difference
    Vect c1, c2; // center of mass
    Vect vc1, vc2; // center of mass velocity
    Vect v1, v2; // average velocity
    Scal dtt; // temporary to reduce
    Scal dt; // dt fluid
    Scal dta; // dt advection
    size_t iter; // iter of fluid solver
    Scal dumpt = -1e10; // last dump time (rounded to nearest dtdump)
    Scal t;
    size_t step;
    size_t dumpn;
    Vect meshpos; // mesh position
    Vect meshvel; // mesh velocity
    Scal ekin, ekin1, ekin2; // kinetic energy
    Scal workst; // work by surface tension
    Vect vlm, vl2; // max-norm and l2-norm of velocity minus "vel"
    Scal pmin, pmax, pd; // pressure min,max
    Scal pavg1, pavg2; // pressure average
    Scal boxomm; // integral of vorticity magnitude over box
    Scal boxomm2; // integral of vorticity magnitude over box
    Vect vomm; // velocity weighted by vorticity
    Scal vommw; // integral of vorticity
    Scal enstr; // enstrophy
    Scal area; // interface area
    Scal dissip, dissip1, dissip2; // energy dissipation rate
    Scal edis, edis1, edis2; // dissipated energy
    std::vector<Vofm> vofm;
    Vect eb_pdrag; // pressure drag on embedded boundaries
    Vect eb_vdrag; // viscous drag
    Vect eb_drag; // total drag
    Scal particles_n; // number of particles
    Scal particles_nrecv; // number of particles received on last communciation

    std::map<std::string, Scal> mst; // map stat
    // Add scalar field for stat.
    void Add(const FieldCell<Scal>& fc, std::string name, M& m) {
      Scal min = std::numeric_limits<Scal>::max();
      Scal max = -std::numeric_limits<Scal>::max();
      Scal sum = 0.;
      for (auto c : m.Cells()) {
        Scal v = fc[c];
        min = std::min(min, v);
        max = std::max(max, v);
        sum += v * m.GetVolume(c);
      }
      m.Reduce(&(mst[name + "_min"] = min), "min");
      m.Reduce(&(mst[name + "_max"] = max), "max");
      m.Reduce(&(mst[name + "_sum"] = sum), "sum");
    }
    // Add component comp of vector field for stat.
    void Add(const FieldCell<Vect>& fc, size_t comp, std::string name, M& m) {
      Scal min = std::numeric_limits<Scal>::max();
      Scal max = -std::numeric_limits<Scal>::max();
      Scal sum = 0.;
      for (auto c : m.Cells()) {
        Scal v = fc[c][comp];
        min = std::min(min, v);
        max = std::max(max, v);
        sum += v * m.GetVolume(c);
      }
      m.Reduce(&(mst[name + "_min"] = min), "min");
      m.Reduce(&(mst[name + "_max"] = max), "max");
      m.Reduce(&(mst[name + "_sum"] = sum), "sum");
    }
    void Print(std::ostream& out) {
      std::string dl = "";
      auto fl = out.flags();
      out.precision(16);
      out << std::scientific << std::setprecision(20);
      size_t i = 0;
      for (auto& it : mst) {
        if (i % 3 == 0) {
          out << "> ";
          dl = "";
        }
        out << dl << it.first << "=" << it.second;
        dl = ", ";
        ++i;
        if (i % 3 == 0) {
          out << std::endl;
        }
      }
      out.flags(fl);
    }
    void Clear() {
      mst.clear();
    }
  };
  Stat st_;
  std::shared_ptr<output::Ser> ost_; // output stat
  Dumper dumper_;
  Dumper dmptraj_; // dumper for traj
  Dumper dmptrep_; // dumper for timer report
  std::unique_ptr<Events> events_; // events from var
  SingleTimer timer_;
  std::unique_ptr<TracerInterface<M>> tracer_;
  std::unique_ptr<ParticlesInterface<M>> particles_;
  std::mt19937 randgen_;
  Scal tracer_dt_;
  Scal particles_dt_;
};

template <class M>
void Hydro<M>::InitEmbed() {
  if (var.Int["enable_embed"]) {
    auto sem = m.GetSem("embed");
    struct {
      FieldNode<Scal> fnl;
    } * ctx(sem);
    if (sem("ctor")) {
      eb_.reset(new Embed<M>(m, var.Double["embed_gradlim"]));
      ctx->fnl = UEmbed<M>::InitEmbed(m, var, m.IsRoot());
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
void Hydro<M>::SpawnTracer() {
  auto& conf = tracer_->GetConf();
  const auto trl = GRange<size_t>(conf.layers);
  const Vect sphere_c(var.Vect["tracer_spawn_sphere_c"]);
  const Scal sphere_r = var.Double["tracer_spawn_sphere_r"];
  const size_t dim = m.GetEdim();
  auto vfcu = tracer_->GetVolumeFraction();
  for (auto l : trl) {
    const std::string prefix = "tracer" + std::to_string(l);
    auto k = var.Double[prefix + "_factor"];
    for (auto c : m.AllCells()) {
      const auto xc = m.GetCenter(c);
      Vect dx = xc - sphere_c;
      if (dim == 2) {
        dx[2] = 0;
      }
      if (dx.sqrnorm() < sqr(sphere_r)) {
        vfcu[l][c] = k;
      }
    }
  }
  tracer_->SetVolumeFraction(vfcu);
}

template <class M>
void Hydro<M>::OverwriteBc() {
  const auto factors = var.Vect["overwrite_inlet_factors"];
  const auto times = var.Vect["overwrite_inlet_times"];
  fassert_equal(factors.size(), times.size());
  if (times.size() == 0) {
    return;
  }
  const Scal t = fs_->GetTime();
  size_t i = 0;
  while (i < times.size() && times[i] <= t) { ++i; }

  auto apply = [&](Scal factor) {
    mebc_fluid_.LoopPairs([&](auto cf_bc) {
      auto& curr = mebc_fluid_[cf_bc.first];
      const auto& orig = mebc_fluid_orig_[cf_bc.first];
      if (curr.type == BCondFluidType::inlet) {
        fassert(orig.type == BCondFluidType::inlet);
        curr.velocity = orig.velocity * factor;
      }
    });
  };

  if (i == 0) { // t < times[0]
    return;
  }

  if (i < times.size()) { // times[i - 1] <= t < times[i]
    const auto t0 = times[i - 1];
    const auto t1 = times[i];
    const auto f0 = factors[i - 1];
    const auto f1 = factors[i];
    const auto f = (t0 < t1 ? f0 + (f1 - f0) * (t - t0) / (t1 - t0) : f0);
    apply(f);
  } else { // t > times[i - 1]
    apply(factors[i - 1]);
  }
}


template <class M>
void Hydro<M>::SpawnParticles(ParticlesView& view) {
  const Vect sphere_c(var.Vect["particles_spawn_sphere_c"]);
  const Scal sphere_r = var.Double["particles_spawn_sphere_r"];
  // particles per unit time
  const size_t dim = m.GetEdim();
  const Scal sphere_vol =
      (dim == 3 ? 4. / 3. * M_PI * std::pow(sphere_r, 3)
                : M_PI * sqr(sphere_r));
  const Vect h = m.GetCellSize();
  const Scal cell_vol = (dim == 3 ? h.prod() : h[0] * h[1]);
  const Vect velocity(var.Vect["particles_spawn_velocity"]);
  const Scal density = var.Double["particles_density"];
  const auto spawn_rate = var.Vect["particles_spawn_rate"];
  const auto diameter = var.Vect["particles_diameter"];
  const auto termvel = var.Vect["particles_termvel"];
  const size_t num_rates = spawn_rate.size();
  fassert_equal(diameter.size(), num_rates);
  fassert_equal(termvel.size(), num_rates);
  std::uniform_real_distribution<Scal> u(0, 1);
  std::uniform_real_distribution<Scal> um(-0.5, 0.5);
  auto& g = randgen_;

  for (auto c : m.Cells()) {
    const auto xc = m.GetCenter(c);
    Vect dx = xc - sphere_c;
    if (dim == 2) {
      dx[2] = 0;
    }
    if (dx.sqrnorm() < sqr(sphere_r)) {
      for (size_t i = 0; i < num_rates; ++i) {
        const Scal prob = particles_dt_ * cell_vol * spawn_rate[i] / sphere_vol;
        if (u(g) < prob) {
          view.x.push_back(m.GetCenter(c) + Vect(um(g), um(g), um(g)) * h);
          view.v.push_back(velocity);
          view.r.push_back(diameter[i] * 0.5);
          view.rho.push_back(density);
          view.termvel.push_back(termvel[i]);
        }
      }
    }
  }
}

template <class M>
void Hydro<M>::InitParticles() {
  if (var.Int["enable_particles"]) {
    typename ParticlesInterface<M>::Conf conf;
    conf.mixture_density = var.Double["rho1"];
    conf.mixture_viscosity = var.Double["mu1"];
    conf.gravity = Vect(var.Vect["gravity"]);

    std::vector<Vect> p_x;
    std::vector<Vect> p_v;
    std::vector<Scal> p_r;
    std::vector<Scal> p_rho;
    std::vector<Scal> p_termvel;
    ParticlesView init{p_x, p_v, p_r, p_rho, p_termvel};
    SpawnParticles(init);
    if (eb_) {
      particles_.reset(new Particles<EB>(m, *eb_, init, fs_->GetTime(), conf));
    } else {
      particles_.reset(new Particles<M>(m, m, init, fs_->GetTime(), conf));
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
    const auto trl = GRange<size_t>(conf.layers);

    conf.density = multi(var.Vect["tracer_density"]);
    fassert(conf.density.size() >= conf.layers);

    conf.viscosity = multi(var.Vect["tracer_viscosity"]);
    fassert(conf.viscosity.size() >= conf.layers);
    conf.gravity = Vect(var.Vect["gravity"]);

    auto termvel = multi(var.Vect["tracer_termvel"]);
    conf.diameter.resize(conf.layers);
    if (var.Int["tracer_use_termvel"]) {
      const Scal mixture_density = var.Double["rho1"];
      for (auto l : trl) {
        conf.diameter[l] = std::sqrt(std::abs(
            18 * conf.viscosity[l] * termvel[l] /
            (conf.gravity.norm() * (conf.density[l] - mixture_density))));
      }
    } else {
      conf.diameter = multi(var.Vect["tracer_diameter"]);
    }
    fassert(conf.diameter.size() >= conf.layers);

    conf.scheme = GetConvSc(var.String["tracer_scheme"]);


    using SlipType = typename TracerInterface<M>::SlipType;
    conf.slip.resize(conf.layers);
    for (auto l : trl) {
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
        slip.velocity =  conf.gravity * (termvel[l] / conf.gravity.norm());
      } else if (type == "constant") {
        slip.type = SlipType::constant;
        arg >> slip.velocity;
      } else {
        throw std::runtime_error(FILELINE + "Unknown slip='" + type + "'");
      }
    }
    Multi<MapEmbed<BCond<Scal>>> vmebc(conf.layers); // boundary conditions
    mebc_fluid_.LoopPairs([&](auto cf_bc) {
      for (auto l : trl) {
        if (l == 0) {
          vmebc[l][cf_bc.first] =
              BCond<Scal>(BCondType::dirichlet, cf_bc.second.nci, 1.);
        } else {
          vmebc[l][cf_bc.first] =
              BCond<Scal>(BCondType::dirichlet, cf_bc.second.nci, 0.);
        }
      }
    });

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

  // XXX ahoc: young velocity
  if (var.Int["youngbc"]) {
    fcyv_.Reinit(m);
    InitYoung();
    for (auto c : m.Cells()) {
      Vect x = m.GetCenter(c);
      fcyv_[c] = GetYoungVel(x);
    }
  }

  std::string fs = var.String["fluid_solver"];
  if (fs == "simple") {
    auto p = ParsePar<Simple<M>>()(var);
    if (eb_) {
      fs_.reset(new Simple<Embed<M>>(
          m, *eb_, fc_vel, mebc_fluid_, mc_velcond_, &fc_rho_, &fc_mu_,
          &fc_force_, &febp_, &fc_src_, &fc_srcm_, 0., st_.dt, p));
    } else {
      fs_.reset(new Simple<M>(
          m, m, fc_vel, mebc_fluid_, mc_velcond_, &fc_rho_, &fc_mu_, &fc_force_,
          &febp_, &fc_src_, &fc_srcm_, 0., st_.dt, p));
    }
  } else if (fs == "proj") {
    auto p = ParsePar<Proj<M>>()(var);
    if (eb_) {
      fs_.reset(new Proj<Embed<M>>(
          m, *eb_, fc_vel, mebc_fluid_, mc_velcond_, &fc_rho_, &fc_mu_,
          &fc_force_, &febp_, &fc_src_, &fc_srcm_, 0., st_.dt, p));
    } else {
      fs_.reset(new Proj<M>(
          m, m, fc_vel, mebc_fluid_, mc_velcond_, &fc_rho_, &fc_mu_, &fc_force_,
          &febp_, &fc_src_, &fc_srcm_, 0., st_.dt, p));
    }
  } else if (fs == "dummy") {
    if (eb_) {
      fs_.reset(new FluidDummy<Embed<M>>(
          m, *eb_, fc_vel, &fc_rho_, &fc_mu_, &fc_force_, &febp_, &fc_src_,
          &fc_srcm_, 0., st_.dt, var));
    } else {
      fs_.reset(new FluidDummy<M>(
          m, m, fc_vel, &fc_rho_, &fc_mu_, &fc_force_, &febp_, &fc_src_,
          &fc_srcm_, 0., st_.dt, var));
    }
  } else {
    throw std::runtime_error("Unknown fluid_solver=" + fs);
  }

  fcbc_ = GetBcField(mf_fluid_, m);
}

template <class M>
void Hydro<M>::InitAdvection(
    const FieldCell<Scal>& fcvf, const FieldCell<Scal>& fccl) {
  std::string as = var.String["advection_solver"];
  if (as == "tvd") {
    auto p = ParsePar<AST>()(var);
    as_.reset(new AST(
        m, fcvf, mf_adv_, &fs_->GetVolumeFlux(Step::time_curr), &fc_src2_, 0.,
        st_.dta, p));
  } else if (as == "vof") {
    if (eb_) {
      auto p = ParsePar<ASVEB>()(var);
      as_.reset(new ASVEB(
          m, *eb_, fcvf, fccl, mf_adv_, &fs_->GetVolumeFlux(Step::time_curr),
          &fc_src2_, 0., st_.dta, p));
    } else {
      auto p = ParsePar<ASV>()(var);
      as_.reset(new ASV(
          m, m, fcvf, fccl, mf_adv_, &fs_->GetVolumeFlux(Step::time_curr),
          &fc_src2_, 0., st_.dta, p));
    }
    layers = GRange<size_t>(1);
  } else if (as == "vofm") {
    if (eb_) {
      auto p = ParsePar<ASVMEB>()(var);
      auto as = new ASVMEB(
          m, *eb_, fcvf, fccl, mf_adv_, &fs_->GetVolumeFlux(Step::time_curr),
          &fc_src2_, 0., st_.dta, p);
      as_.reset(as);
      layers = GRange<size_t>(as->GetNumLayers());
    } else {
      auto p = ParsePar<ASVM>()(var);
      auto as = new ASVM(
          m, m, fcvf, fccl, mf_adv_, &fs_->GetVolumeFlux(Step::time_curr),
          &fc_src2_, 0., st_.dta, p);
      as_.reset(as);
      layers = GRange<size_t>(as->GetNumLayers());
    }
  } else {
    throw std::runtime_error("Unknown advection_solver=" + as);
  }

  st_.vofm.resize(layers.size());
  fck_.resize(layers);
  fck_.InitAll(FieldCell<Scal>(m, GetNan<Scal>()));
  auto ps = ParsePar<PartStr<Scal>>()(m.GetCellSize().norminf(), var);
  psm_par_ = ParsePar<PartStrMeshM<M>>()(ps, var);
}

template <class M>
void Hydro<M>::InitStat() {
  if (m.IsRoot()) {
    // Stat: var.Double[p] with name n
    /*
    auto on = [this](std::string n, std::string p) {
      return std::make_shared<output::OutScalFunc<Scal>>(
          n, [&,p](){ return var.Double[p]; });
    };
    */

    // Stat: *p with name n
    auto op = [](std::string n, Scal* p) {
      return std::make_shared<output::OutScalFunc<Scal>>(
          n, [p]() { return *p; });
    };

    auto& s = st_;
    output::VOut con = {
        op("t", &s.t),
        std::make_shared<output::OutScalFunc<int>>(
            "iter", [this]() { return st_.iter; }),
        op("dt", &s.dt),
        op("dta", &s.dta),
        std::make_shared<output::OutScalFunc<Scal>>(
            "diff", [this]() { return fs_->GetError(); }),
        std::make_shared<output::OutScalFunc<Scal>>(
            "wt", [this]() { return timer_.GetSeconds(); }),
        op("m1", &s.m1),
        op("m2", &s.m2),
        op("m20", &s.m20),
        op("m2d", &s.m2d),
        op("ekin", &s.ekin),
        op("ekin1", &s.ekin1),
        op("ekin2", &s.ekin2),
        op("c1x", &s.c1[0]),
        op("c1y", &s.c1[1]),
        op("c1z", &s.c1[2]),
        op("c2x", &s.c2[0]),
        op("c2y", &s.c2[1]),
        op("c2z", &s.c2[2]),
        op("vc1x", &s.vc1[0]),
        op("vc1y", &s.vc1[1]),
        op("vc1z", &s.vc1[2]),
        op("vc2x", &s.vc2[0]),
        op("vc2y", &s.vc2[1]),
        op("vc2z", &s.vc2[2]),
        op("v1x", &s.v1[0]),
        op("v1y", &s.v1[1]),
        op("v1z", &s.v1[2]),
        op("v2x", &s.v2[0]),
        op("v2y", &s.v2[1]),
        op("v2z", &s.v2[2]),
        op("meshposx", &s.meshpos[0]),
        op("meshvelx", &s.meshvel[0]),
        op("meshvelz", &s.meshvel[2]),
    };
    if (var.Int["statbox"]) {
      con.push_back(op("boxomm", &s.boxomm));
      con.push_back(op("boxomm2", &s.boxomm2));
      con.push_back(op("vommz", &s.vomm[2]));
    }
    if (var.Int["statvel"]) {
      con.push_back(op("vlmx", &s.vlm[0]));
      con.push_back(op("vlmy", &s.vlm[1]));
      con.push_back(op("vlmz", &s.vlm[2]));
      con.push_back(op("vl2x", &s.vl2[0]));
      con.push_back(op("vl2y", &s.vl2[1]));
      con.push_back(op("vl2z", &s.vl2[2]));
    }
    con.push_back(op("pmin", &s.pmin));
    con.push_back(op("pmax", &s.pmax));
    con.push_back(op("pd", &s.pd));
    con.push_back(op("pavg1", &s.pavg1));
    con.push_back(op("pavg2", &s.pavg2));
    if (var.Int["enstrophy"]) {
      con.push_back(op("enstr", &s.enstr));
    }
    if (var.Int["stat_area"]) {
      con.push_back(op("area", &s.area));
    }
    if (var.Int["stat_dissip"]) {
      con.push_back(op("dissip", &s.dissip));
      con.push_back(op("dissip1", &s.dissip1));
      con.push_back(op("dissip2", &s.dissip2));
      con.push_back(op("edis", &s.edis));
      con.push_back(op("edis1", &s.edis1));
      con.push_back(op("edis2", &s.edis2));
    }
    if (var.Int["stat_vofm"]) {
      s.vofm.resize(layers.size());
      for (auto l : layers) {
        auto sl = std::to_string(l);
        con.push_back(op("vofm_cells_vf" + sl, &s.vofm[l].cells_vf));
        con.push_back(op("vofm_cells_cl" + sl, &s.vofm[l].cells_cl));
        con.push_back(op("vofm_sum_vf" + sl, &s.vofm[l].sum_vf));
        con.push_back(op("vofm_hist" + std::to_string(l + 1), &s.vofm[l].hist));
      }
    }
    if (eb_) {
      con.push_back(op("pdragx", &s.eb_pdrag[0]));
      con.push_back(op("pdragy", &s.eb_pdrag[1]));
      con.push_back(op("pdragz", &s.eb_pdrag[2]));
      con.push_back(op("vdragx", &s.eb_vdrag[0]));
      con.push_back(op("vdragy", &s.eb_vdrag[1]));
      con.push_back(op("vdragz", &s.eb_vdrag[2]));
      con.push_back(op("dragx", &s.eb_drag[0]));
      con.push_back(op("dragy", &s.eb_drag[1]));
      con.push_back(op("dragz", &s.eb_drag[2]));
    }
    if (particles_) {
      con.push_back(op("particles_n", &s.particles_n));
      con.push_back(op("particles_nrecv", &s.particles_nrecv));
    }
    ost_ = std::make_shared<output::SerScalPlain<Scal>>(con, "stat.dat");
  }
}

template <class M>
void Hydro<M>::Init() {
  using namespace fluid_condition;
  auto sem = m.GetSem("init");
  struct {
    FieldCell<Vect> fcvel; // initial velocity
    FieldCell<Scal> fcvf; // initial volume fraction
    FieldCell<Scal> fccl; // initial color
    FieldCell<Scal> fcbody;
    FieldCell<bool> fcbodymask;
    Multi<FieldCell<Scal>> tracer_vfcu;
    Vars varbody;
  } * ctx(sem);
  auto& fcvel = ctx->fcvel;
  auto& fcvf = ctx->fcvf;
  auto& fccl = ctx->fccl;
  if (sem("flags")) {
    m.flags.linreport = var.Int["linreport"];
    m.flags.check_symmetry = var.Int["check_symmetry"];
    m.flags.check_symmetry_dump_threshold =
        var.Double["check_symmetry_dump_threshold"];
    m.flags.is_periodic[0] = var.Int["hypre_periodic_x"];
    m.flags.is_periodic[1] = var.Int["hypre_periodic_y"];
    m.flags.is_periodic[2] = var.Int["hypre_periodic_z"];
    randgen_.seed(m.GetId() + 1);
  }
  if (sem.Nested("embed")) {
    InitEmbed();
  }
  if (sem.Nested()) {
    InitVf(fcvf, var, m);
  }
  if (sem.Nested()) {
    if (var.Int["enable_tracer"]) {
      InitTracerFields(ctx->tracer_vfcu);
    }
  }
  if (sem("fields")) {
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
    auto isig = CreateInitSig<M>(var);
    isig(fc_sig_, m);
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
        throw std::runtime_error(FILELINE + ": Unknown init_contang=" + name);
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

    // global mesh size
    MIdx gs = m.GetGlobalSize();

    if (m.IsRoot()) {
      std::cout << "global mesh=" << gs << std::endl;
      std::cout << "surface tension dt=" << GetStDt() << std::endl;
      std::cout << "viscosity dt=" << GetVisDt() << std::endl;
    }

    // boundary conditions
    if (var.Int["enable_bc"]) {
      if (eb_) {
        auto p = InitBc(var, *eb_);
        mebc_fluid_ = std::get<0>(p);
        mf_adv_ = std::get<1>(p);
        me_group_ = std::get<2>(p);
        vdesc_ = std::get<3>(p);
        mf_adv_.LoopBCond(*eb_, [&](auto cf, auto c, auto) { //
          me_contang_[cf] = fc_contang_[c];
          mf_adv_[cf].contang = fc_contang_[c];
        });
      } else {
        auto p = InitBc(var, m);
        mebc_fluid_ = std::get<0>(p);
        mf_adv_ = std::get<1>(p);
        me_group_ = std::get<2>(p);
        vdesc_ = std::get<3>(p);
        for (auto& p : mf_adv_.GetMapFace()) {
          const auto& f = p.first;
          auto& bc = p.second;
          const auto c = m.GetCell(f, bc.nci);
          me_contang_[f] = fc_contang_[c];
          bc.contang = fc_contang_[c];
        }
      }
      mf_fluid_ = GetCondFluid<M>(mebc_fluid_);
    } else {
      GetFluidFaceCond(var, m, mf_fluid_, mf_adv_);
      mebc_fluid_ = GetBCondFluid<M>(mf_fluid_);
    }
    mebc_fluid_orig_ = mebc_fluid_;

    // boundary conditions for smoothing of volume fraction
    for (auto& p : mebc_fluid_.GetMapFace()) {
      const IdxFace f = p.first;
      auto& bc = p.second;
      const auto nci = bc.nci;
      if (bc.type == BCondFluidType::symm) {
        mf_cond_vfsm_[f].template Set<CondFaceReflect>(nci);
      } else {
        mf_cond_vfsm_[f].template Set<CondFaceGradFixed<Scal>>(Scal(0), nci);
      }
    }
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

  if (var.Int["initvort"] && sem.Nested("initvort")) {
    InitVort(fcvel, fcvel, mf_fluid_, var.Int["vort_report"], m);
  }

  if (var.Int["wavelamb_vort"] && sem("wavelamb")) {
    FieldCell<Vect> fc(m);
    Vars vr;
    vr.String.Set("vel_init", "wavelamb");
    vr.Double.Set("wavelamb_a0", var.Double["wavelamb_a0"]);
    vr.Double.Set("wavelamb_xc", var.Double["wavelamb_xc"]);
    vr.Double.Set("wavelamb_h", var.Double["wavelamb_h"]);
    vr.Double.Set("wavelamb_k", var.Double["wavelamb_k"]);
    vr.Vect.Set("gravity", var.Vect["gravity"]);
    InitVel(fc, vr, m);
    for (auto c : m.AllCells()) {
      if (fc[c].sqrnorm() != 0) {
        fcvel[c] = fc[c];
      }
    }
  }

  if (var.Int["vel_init_noise"] && sem("noise")) {
    Vect vel0(var.Vect["noise_vel0"]);
    Vect vel1(var.Vect["noise_vel1"]);
    Vect vel2(var.Vect["noise_vel2"]);
    Vect per0(var.Vect["noise_per0"]);
    Vect per1(var.Vect["noise_per1"]);
    Vect per2(var.Vect["noise_per2"]);
    Vect k0 = Vect(2 * M_PI) / (per0 * m.GetCellSize());
    Vect k1 = Vect(2 * M_PI) / (per1 * m.GetCellSize());
    Vect k2 = Vect(2 * M_PI) / (per2 * m.GetCellSize());
    for (auto c : m.AllCells()) {
      auto x = m.GetCenter(c);
      fcvel[c] += vel0 * std::sin(k0.dot(x));
      fcvel[c] += vel1 * std::sin(k1.dot(x));
      fcvel[c] += vel2 * std::sin(k2.dot(x));
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
    Smoothen(fcvf, mf_cond_vfsm_, m, var.Int["vf_init_sm"]);
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

  if (sem("body-mask")) {
    auto& varbody = ctx->varbody;
    varbody.String.Set("init_vf", var.String["body_init"]);
    varbody.String.Set("list_path", var.String["body_list_path"]);
    varbody.Int.Set("dim", var.Int["dim"]);
    varbody.Int.Set("list_ls", 3);
  }
  if (sem.Nested("body-mask")) {
    InitVf(ctx->fcbody, ctx->varbody, m);
  }
  if (sem("body-bc")) {
    // Step-wise approximation of bodies
    const Scal clear0 = var.Double["bcc_clear0"];
    const Scal clear1 = var.Double["bcc_clear1"];
    const Scal inletcl = var.Double["inletcl"];
    const Scal fill_vf = var.Double["bcc_fill"];
    auto& fc = ctx->fcbodymask;
    fc.Reinit(m, false);
    for (auto c : m.AllCells()) {
      fc[c] = (ctx->fcbody[c] > 0.5);
    }
    if (var.Int["body_init_inverse"]) {
      for (auto c : m.AllCells()) {
        fc[c] = !fc[c];
      }
    }
    AppendBodyCond<M>(
        fc, var.String["body_bc"], m, clear0, clear1, inletcl, fill_vf, nullptr,
        mf_fluid_, mf_adv_);
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

    InitTracer(ctx->tracer_vfcu);

    InitParticles();

    st_.iter = 0;
    st_.t = fs_->GetTime();

    if (m.IsLead()) {
      this->var_mutable.Int.Set("iter", st_.iter);
    }

    InitStat();

    if (var.Int["fill_halo_nan"]) {
      std::vector<std::pair<IdxFace, size_t>> vf;
      for (auto& p : mf_adv_.GetMapFace()) {
        vf.emplace_back(p.first, p.second.GetNci());
      }
      m.SetNanFaces(vf);
      m.flags.nan_faces_value = var.Double["fill_halo_nan_value"];
    }

    if (m.IsLead()) {
      events_ = std::unique_ptr<Events>(
          new Events(this->var_mutable, m.IsRoot(), m.IsLead()));
      events_->Parse();
    }
  }
  if (var.Int["dumpbc"]) {
    if (vdesc_.size()) {
      if (sem("dump-bcgroups")) {
        if (m.IsRoot()) {
          std::ofstream fdesc("bc_groups.dat");
          for (size_t i = 0; i < vdesc_.size(); ++i) {
            fdesc << i << " " << vdesc_[i] << std::endl;
          }
        }
      }
      if (sem.Nested("bcdump")) {
        if (eb_) {
          UInitEmbedBc<M>::DumpPoly("bc.vtk", me_group_, me_contang_, *eb_, m);
        } else {
          UInitEmbedBc<M>::DumpPoly("bc.vtk", me_group_, me_contang_, m, m);
        }
      }
    } else {
      if (sem.Nested()) {
        DumpBcFaces(mf_adv_, mf_fluid_, "bc.vtk", m);
      }
    }
  }
  if (eb_ && sem.Nested()) {
    eb_->DumpPoly(var.Int["vtkbin"], var.Int["vtkmerge"]);
  }
  if (var.Int["dumpinit"]) {
    if (sem.Nested()) {
      Dump(true);
    }
  }
}

template <class M>
Hydro<M>::Hydro(Vars& var0, const MyBlockInfo& bi, Par& par)
    : KernelMeshPar<M, Par>(var0, bi, par)
    , st_{}
    , dumper_(var, "dump_field_")
    , dmptraj_(var, "dump_traj_")
    , dmptrep_(var, "dump_trep_") {}

template <class M>
void Hydro<M>::CalcStat() {
  auto sem = m.GetSem("stat");

  auto& s = st_;
  auto& fa = as_->GetField();
  auto& fv = fs_->GetVelocity();
  auto& fp = fs_->GetPressure();

  if (sem("stat-add")) {
    s.Add(fa, "vf", m);
    s.Add(fv, 0, "vx", m);
    s.Add(fv, 1, "vy", m);
    s.Add(fv, 2, "vz", m);
    s.Add(fp, "p", m);
  }
  if (sem("stat-print")) {
    if (var.Int["report_stat"] && m.IsRoot()) {
      s.Print(std::cout);
    }
    s.Clear();
  }

  if (sem("local")) {
    // Store vc1 and vc2
    s.vc1 = s.c1;
    s.vc2 = s.c2;

    // mass, center, velocity, kinetic energy
    s.m1 = 0;
    s.m2 = 0;
    s.c1 = Vect(0);
    s.c2 = Vect(0);
    s.v1 = Vect(0);
    s.v2 = Vect(0);
    s.ekin = 0;
    s.ekin1 = 0;
    s.ekin2 = 0;
    for (auto i : m.Cells()) {
      Scal o = m.GetVolume(i);
      Scal a2 = fa[i];
      Scal a1 = 1. - a2;
      Vect v = fv[i];
      Vect x = m.GetCenter(i);

      s.m1 += a1 * o;
      s.m2 += a2 * o;
      s.c1 += x * (a1 * o);
      s.c2 += x * (a2 * o);
      s.v1 += v * (a1 * o);
      s.v2 += v * (a2 * o);
      s.ekin += 0.5 * v.dot(v) * fc_rho_[i] * o;
      s.ekin1 += 0.5 * v.dot(v) * fc_rho_[i] * o * a1;
      s.ekin2 += 0.5 * v.dot(v) * fc_rho_[i] * o * a2;
    }

    m.Reduce(&s.m1, "sum");
    m.Reduce(&s.m2, "sum");
    m.Reduce(&s.ekin, "sum");
    m.Reduce(&s.ekin1, "sum");
    m.Reduce(&s.ekin2, "sum");
    for (size_t d = 0; d < dim; ++d) {
      m.Reduce(&s.c1[d], "sum");
      m.Reduce(&s.c2[d], "sum");
      m.Reduce(&s.v1[d], "sum");
      m.Reduce(&s.v2[d], "sum");
    }

    // pressure
    s.pmin = 1e10;
    s.pmax = -1e10;
    s.pavg1 = 0.;
    s.pavg2 = 0.;
    for (auto c : m.Cells()) {
      Scal o = m.GetVolume(c);
      Scal a2 = fa[c];
      Scal p = fp[c];
      s.pmin = std::min(s.pmin, p);
      s.pmax = std::max(s.pmax, p);
      s.pavg1 += p * (1. - a2) * o;
      s.pavg2 += p * a2 * o;
    }
    m.Reduce(&s.pmin, "min");
    m.Reduce(&s.pmax, "max");
    m.Reduce(&s.pavg1, "sum");
    m.Reduce(&s.pavg2, "sum");

    if (var.Int["statvel"]) {
      s.vlm = Vect(0);
      s.vl2 = Vect(0);
      Vect v0(var.Vect["vel"]);
      for (auto c : m.Cells()) {
        Scal o = m.GetVolume(c);
        Scal a2 = fa[c];
        Vect v = fv[c];
        auto dv = (v - v0).abs();
        for (size_t d = 0; d < dim; ++d) {
          s.vlm[d] = std::max(s.vlm[d], dv[d] * a2);
          s.vl2[d] += sqr(dv[d]) * o * a2;
        }
      }
      for (size_t d = 0; d < dim; ++d) {
        m.Reduce(&s.vlm[d], "max");
        m.Reduce(&s.vl2[d], "sum");
      }
    }

    if (var.Int["statbox"] || var.Int["enstrophy"]) {
      CalcVort();
    }

    // XXX: adhoc: also controls s.vomm
    if (var.Int["statbox"]) {
      // XXX requires CalcVort above
      Vect xa(var.Vect["statboxa"]);
      Vect xb(var.Vect["statboxb"]);
      Vect x2a(var.Vect["statbox2a"]);
      Vect x2b(var.Vect["statbox2b"]);
      Vect h = m.GetCellSize();
      size_t dm = (xb - xa).argmin();
      size_t dm2 = (x2b - x2a).argmin();
      // box size at least h
      for (size_t d = 0; d < dim; ++d) {
        xb[d] = std::max<Scal>(xb[d], xa[d] + h[d]);
        x2b[d] = std::max<Scal>(x2b[d], x2a[d] + h[d]);
      }
      // integrate
      s.boxomm = 0.;
      s.boxomm2 = 0.;
      for (auto c : m.Cells()) {
        auto xc = m.GetCenter(c);
        if (xa <= xc && xc <= xb) {
          s.boxomm += m.GetVolume(c) * fcom_[c][dm];
        }
        if (x2a <= xc && xc <= x2b) {
          s.boxomm2 += m.GetVolume(c) * fcom_[c][dm2];
        }
      }
      // divide by smallest dimension to get integral over slice
      s.boxomm /= h[dm];
      s.boxomm2 /= h[dm2];
      m.Reduce(&s.boxomm, "sum");
      m.Reduce(&s.boxomm2, "sum");

      // vomm
      s.vomm = Vect(0);
      s.vommw = 0;
      for (auto c : m.Cells()) {
        s.vomm += fv[c] * fcomm_[c];
        s.vommw += fcomm_[c];
      }
      for (size_t d = 0; d < dim; ++d) {
        m.Reduce(&s.vomm[d], "sum");
      }
      m.Reduce(&s.vommw, "sum");
    }

    if (var.Int["enstrophy"]) {
      // XXX requires CalcVort above
      s.enstr = 0.;
      for (auto c : m.Cells()) {
        s.enstr += 0.5 * sqr(fcomm_[c]) * fc_rho_[c] * m.GetVolume(c);
      }
      m.Reduce(&s.enstr, "sum");
    }
    // surface area
    if (var.Int["stat_area"]) {
      s.area = 0;
      if (auto as = dynamic_cast<ASVM*>(as_.get())) {
        using R = Reconst<Scal>;
        auto& fcn = *as->GetNormal()[0];
        auto& fca = *as->GetAlpha()[0];
        auto& fcvf = as->GetField(Step::time_curr, 0);
        Vect h = m.GetCellSize();
        for (auto c : m.Cells()) {
          if (fcvf[c] > 0. && fcvf[c] < 1. && !IsNan(fca[c])) {
            auto xx = R::GetCutPoly2(fcn[c], fca[c], h);
            Scal ar = std::abs(R::GetArea(xx, fcn[c]));
            s.area += ar;
          }
        }
        if (IsNan(s.area)) {
          s.area = 0;
        }
        m.Reduce(&s.area, "sum");
      }
    }
    if (var.Int["stat_dissip"]) {
      s.dissip = 0.;
      s.dissip1 = 0.;
      s.dissip2 = 0.;
      CalcStrain(fs_->GetVelocity(), fc_strain_);
      auto& fcmu = fc_mu_;
      auto& fcd = fc_strain_;
      auto& fcvf = as_->GetField();
      for (auto c : m.Cells()) {
        Scal o = m.GetVolume(c);
        Scal mu = fcmu[c];
        Scal d = fcd[c];
        Scal vf = fcvf[c];
        s.dissip += 2. * mu * d * o;
        s.dissip1 += 2. * mu * d * (1. - vf) * o;
        s.dissip2 += 2. * mu * d * vf * o;
      }
      m.Reduce(&s.dissip, "sum");
      m.Reduce(&s.dissip1, "sum");
      m.Reduce(&s.dissip2, "sum");
    }
    if (var.Int["stat_vofm"]) {
      if (auto as = dynamic_cast<ASVM*>(as_.get())) {
        // clear
        for (auto l : layers) {
          auto& v = s.vofm[l];
          v.cells_vf = 0;
          v.cells_cl = 0;
          v.sum_vf = 0;
          v.hist = 0;
        }
        for (auto l : layers) {
          auto& v = s.vofm[l];
          auto& vf = *as->GetFieldM()[l];
          auto& cl = *as->GetColor()[l];
          for (auto c : m.Cells()) {
            v.cells_vf += (vf[c] > 0 ? 1 : 0);
            v.cells_cl += (cl[c] != kClNone ? 1 : 0);
            v.sum_vf += vf[c] * m.GetVolume(c);
          }
          m.Reduce(&v.cells_vf, "sum");
          m.Reduce(&v.cells_cl, "sum");
          m.Reduce(&v.sum_vf, "sum");
        }
        for (auto c : m.Cells()) {
          size_t cnt = 0;
          for (auto l : layers) {
            if ((*as->GetFieldM()[l])[c] > 0) {
              ++cnt;
            }
          }
          if (cnt > 0) {
            assert(cnt - 1 >= 0 && cnt - 1 < layers.size());
            s.vofm[cnt - 1].hist += 1;
          }
        }
        for (auto l : layers) {
          m.Reduce(&s.vofm[l].hist, "sum");
        }
      }
    }
    if (eb_) {
      st_.eb_pdrag = CalcPressureDrag(fs_->GetPressure(), *eb_);
      st_.eb_vdrag = CalcViscousDrag(fs_->GetVelocity(), fc_mu_, *eb_);
      st_.eb_drag = st_.eb_pdrag + st_.eb_vdrag;
      if (m.GetEdim() == 2) {
        const Scal h = m.GetCellSize()[2];
        st_.eb_pdrag /= h;
        st_.eb_vdrag /= h;
        st_.eb_drag /= h;
      }
      for (size_t d = 0; d < dim; ++d) {
        m.Reduce(&st_.eb_pdrag[d], "sum");
        m.Reduce(&st_.eb_vdrag[d], "sum");
        m.Reduce(&st_.eb_drag[d], "sum");
      }
    }
    if (particles_) {
      st_.particles_n = particles_->GetView().x.size();
      st_.particles_nrecv = particles_->GetNumRecv();
      m.Reduce(&st_.particles_n, "sum");
      m.Reduce(&st_.particles_nrecv, "sum");
    }
  }

  if (sem("reduce")) {
    Scal im1 = (s.m1 == 0 ? 0. : 1. / s.m1);
    Scal im2 = (s.m2 == 0 ? 0. : 1. / s.m2);
    s.c1 *= im1;
    s.c2 *= im2;
    s.v1 *= im1;
    s.v2 *= im2;
    s.pavg1 *= im1;
    s.pavg2 *= im2;

    if (s.m20 == 0.) {
      s.m20 = s.m2;
    } else {
      s.m2d = (s.m2 - s.m20) / s.m20;
    }

    if (s.vommw != 0) {
      s.vomm /= s.vommw;
    }

    Scal dt = fs_->GetTimeStep();

    // Moving mesh
    s.c1 += st_.meshpos;
    s.c2 += st_.meshpos;
    s.vc1 = (s.c1 - s.vc1) / dt;
    s.vc2 = (s.c2 - s.vc2) / dt;

    // dissipated energy
    s.edis += s.dissip * dt;
    s.edis1 += s.dissip1 * dt;
    s.edis2 += s.dissip2 * dt;

    if (const std::string* s = var.String.Find("meshvel_auto")) {
      auto upd = [this, s](Vect& meshvel) {
        Vect vel(0);
        if (*s == "v") {
          vel = st_.v2;
        } else if (*s == "vc") {
          vel = st_.vc2;
        } else if (*s == "vomm") {
          vel = st_.vomm;
        } else {
          throw std::runtime_error("Unknown meshvel_auto=" + *s);
        }
        Vect mask(var.Vect["meshvel_mask"]);
        vel *= mask;
        double w = var.Double["meshvel_weight"];
        meshvel = vel * w + meshvel * (1. - w);
        st_.meshvel = meshvel;
        st_.meshpos += st_.meshvel * st_.dt;
      };
      if (auto fs = dynamic_cast<Simple<M>*>(fs_.get())) {
        auto par = fs->GetPar();
        upd(par.meshvel);
        fs->SetPar(par);
      } else if (auto fs = dynamic_cast<Proj<M>*>(fs_.get())) {
        auto par = fs->GetPar();
        upd(par.meshvel);
        fs->SetPar(par);
      }
    }

    if (var.Int["statvel"]) {
      for (size_t d = 0; d < dim; ++d) {
        s.vl2[d] = std::sqrt(s.vl2[d] * im2);
      }
      s.pd = s.pmax - s.pmin;
    }
  }

  if (sem("vfslip")) {
    auto kslip = var.Double["kslip"];
    if (kslip != 0) {
      Vect slipvel(var.Vect["slipvel"]);
      // XXX: adhoc, overwrite wall conditions
      auto& fa = as_->GetField();
      for (auto& p : mebc_fluid_.GetMapFace()) {
        const IdxFace f = p.first;
        auto& bc = p.second;
        if (bc.type == BCondFluidType::wall) {
          const IdxCell c = m.GetCell(f, bc.nci);
          bc.velocity = slipvel * kslip * fa[c];
        }
      }
    }

    // Slip velocity penalization.
    const Scal penalslip = var.Double["penalslip"];
    if (penalslip != 0) {
      const Scal dt = fs_->GetTimeStep();
      const Vect slipvel(var.Vect["slipvel"]);
      // const auto& fa = as_->GetField();
      const auto& fa = fc_smvf_;
      const auto& fv = fs_->GetVelocity();
      for (auto& p : mebc_fluid_.GetMapFace()) {
        const IdxFace f = p.first;
        const auto& bc = p.second;
        if (bc.type == BCondFluidType::wall) {
          const IdxCell c = m.GetCell(f, bc.nci);
          Scal sgn = (slipvel - fv[c]).dot(slipvel);
          if (sgn > 0) {
            fc_force_[c] +=
                (slipvel - fv[c]) * (fc_rho_[c] * penalslip * fa[c] / dt);
          }
        }
      }
    }
    // Repulsive force from walls.
    const Scal slipnormal = var.Double["slipnormal"];
    if (slipnormal != 0) {
      const Scal slipnormal_dist = var.Double["slipnormal_dist"];
      const Scal dt = fs_->GetTimeStep();
      const Scal h = m.GetCellSize()[0];
      const auto& fcvf = fc_smvf_;
      if (eb_) {
        auto& eb = *eb_;
        fc_dist_.Reinit(eb, 0);
        for (auto c : eb.SuCells()) {
          fc_dist_[c] = eb.GetSignedDistance(c);
        }
        fc_phi_.Reinit(eb, 0); // potential [length]
        for (auto c : eb.SuCells()) {
          if (!IsNan(fc_dist_[c])) {
            const Scal d0 = slipnormal_dist * h;
            const Scal d = std::max(0., d0 - fc_dist_[c]);
            fc_phi_[c] += slipnormal * d;
          }
        }
        const auto ffg = UEB::Gradient(fc_phi_, {}, eb); // potential grad [-]
        const auto ff_rho = UEB::InterpolateHarmonic(fc_rho_, {}, eb);
        if (auto as = dynamic_cast<ASVMEB*>(as_.get())) {
          for (auto f : eb.Faces()) {
            const IdxCell cm = m.GetCell(f, 0);
            const IdxCell cp = m.GetCell(f, 1);
            if (eb.GetType(cm) == EB::Type::excluded ||
                eb.GetType(cp) == EB::Type::excluded) {
              continue;
            }
            const auto& fccl = as->GetColor();
            const auto& fcu = as->GetFieldM();
            std::set<Scal> colors;
            for (auto l : layers) {
              const Scal clm = (*fccl[l])[cm];
              const Scal clp = (*fccl[l])[cp];
              if (clm != kClNone) colors.insert(clm);
              if (clp != kClNone) colors.insert(clp);
            }
            for (auto cl : colors) {
              Scal um = 0;
              Scal up = 0;
              for (auto l : layers) {
                if ((*fccl[l])[cm] == cl) um = (*fcu[l])[cm];
                if ((*fccl[l])[cp] == cl) up = (*fcu[l])[cp];
              }
              const Scal uf = (um + up) * 0.5;
              if (uf != 0) {
                febp_[f] -= ff_rho[f] * ffg[f] * uf * h / sqr(dt);
              }
            }
          }
        }
      } else {
        for (auto& p : mebc_fluid_.GetMapFace()) {
          const IdxFace f = p.first;
          const auto& bc = p.second;
          const auto nci = bc.nci;
          if (bc.type == BCondFluidType::wall) {
            const IdxCell c = m.GetCell(f, bc.nci);
            fc_force_[c] += m.GetNormal(f) * ((nci == 1 ? 1 : -1) * fc_rho_[c] *
                                              slipnormal * fcvf[c]);
          }
        }
      }
    }
  }
  if (sem("young")) {
    if (var.Int["youngbc"]) {
      InitYoung();
      for (auto& p : mebc_fluid_.GetMapFace()) {
        const IdxFace f = p.first;
        auto& bc = p.second;
        if (bc.type == BCondFluidType::wall) {
          const Vect x = m.GetCenter(f);
          bc.velocity = GetYoungVel(x);
        }
      }
    }
  }
}

template <class M>
void Hydro<M>::CalcDt() {
  auto sem = m.GetSem("calcdt");
  struct {
    Scal dtmin;
  } * ctx(sem);

  if (sem("local")) {
    st_.t = fs_->GetTime();
    ctx->dtmin = fs_->GetAutoTimeStep();
    m.Reduce(&ctx->dtmin, "min");
  }
  if (sem("reduce")) {
    // set from cfl if defined
    if (auto* cfl = var.Double.Find("cfl")) {
      st_.dt = ctx->dtmin * (*cfl);
      st_.dt = std::min<Scal>(st_.dt, var.Double["dtmax"]);
    }

    // constraint from surface tension
    st_.dt = std::min<Scal>(st_.dt, GetStDt());

    // constraint from viscosity
    st_.dt = std::min<Scal>(st_.dt, GetVisDt());

    fs_->SetTimeStep(st_.dt);

    // set from cfla if defined
    if (auto* cfla = var.Double.Find("cfla")) {
      st_.dta = ctx->dtmin * (*cfla);
      st_.dta = std::min<Scal>(st_.dta, var.Double["dtmax"]);
    }
    // round up dta to such that dt / dta is integer
    const Scal dt = fs_->GetTime() + fs_->GetTimeStep() - as_->GetTime();
    st_.dta = dt / std::max(1, int(dt / st_.dta + 0.5));
    as_->SetTimeStep(st_.dta);

    if (tracer_) {
      if (auto* cflt = var.Double.Find("cflt")) {
        tracer_dt_ = ctx->dtmin * (*cflt);
        // round up dta to such that dt / dta is integer
        const Scal dtwhole =
            fs_->GetTime() + fs_->GetTimeStep() - tracer_->GetTime();
        tracer_dt_ = dtwhole / std::max(1, int(dtwhole / tracer_dt_ + 0.5));
      } else {
        tracer_dt_ = fs_->GetTimeStep();
      }
    }

    if (particles_) {
      if (auto* cflp = var.Double.Find("cflp")) {
        particles_dt_ = ctx->dtmin * (*cflp);
        // round up dta to such that dt / dta is integer
        const Scal dtwhole =
            fs_->GetTime() + fs_->GetTimeStep() - particles_->GetTime();
        particles_dt_ =
            dtwhole / std::max(1, int(dtwhole / particles_dt_ + 0.5));
      } else {
        particles_dt_ = fs_->GetTimeStep();
      }
    }
  }
  if (sem()) {
    // FIXME: empty stage
  }
}

template <class M>
void Hydro<M>::CalcMixture(const FieldCell<Scal>& fc_vf0) {
  auto sem = m.GetSem("mixture");

  if (sem("init")) {
    fc_mu_.Reinit(m);
    fc_rho_.Reinit(m);
    fc_force_.Reinit(m, Vect(0));
    febp_.Reinit(m, 0);
    fc_smvf_ = fc_vf0;

    // XXX: oscillating source
    Scal source_mag = var.Double["source_mag"];
    if (source_mag != 0) {
      Scal mag = source_mag;
      Scal mag0 = var.Double["source_mag0"]; // constant component
      Scal freq = var.Double["source_freq"];
      Scal wly = var.Double["source_wly"];
      Scal mexp = var.Double["source_mexp"];
      Scal mexp0 = var.Double["source_mexp0"];
      Scal pi = M_PI;
      Scal s = std::sin(st_.t * freq * 2. * pi);
      fc_src_.Reinit(m, 0);
      for (auto c : m.Cells()) {
        const Vect x = m.GetCenter(c);
        const Scal sy = std::cos(2. * pi * x[1] / wly);
        const Scal vf = fc_vf0[c];
        Scal q = s * sy * mag * vf;
        Scal q0 = mag0 * vf;
        if (st_.m20 > 0 && st_.m2 > 0) {
          q *= std::pow(st_.m2 / st_.m20, mexp);
          q0 *= std::pow(st_.m2 / st_.m20, mexp0);
        }
        fc_src_[c] = q + q0;
      }
    }

    // XXX: oscillating force
    Scal force_mag = var.Double["force_mag"];
    if (force_mag != 0) {
      Scal force_freq = var.Double["force_freq"];
      Scal force_wly = var.Double["force_wly"];
      Scal pi = M_PI;
      Scal s = std::sin(st_.t * force_freq * 2. * pi);
      for (auto c : m.Cells()) {
        Vect x = m.GetCenter(c);
        Scal sy = std::cos(2. * pi * x[1] / force_wly);
        fc_force_[c][0] += s * sy * force_mag * fc_vf0[c];
      }
    }
  }

  if (sem.Nested("smooth")) {
    Smoothen(fc_smvf_, mf_cond_vfsm_, m, var.Int["vfsmooth"]);
  }

  if (sem.Nested("smoothnode")) {
    SmoothenNode(fc_smvf_, m, var.Int["vfsmoothnode"]);
  }

  if (sem("calc")) {
    FieldCell<Scal>& a = fc_smvf_;
    FieldFace<Scal>& af = ff_smvf_;
    if (eb_) {
      auto& eb = *eb_;
      // Compute volume fraction relative to cut cell volume
      // as needed for density and viscosity.
      for (auto c : eb.AllCells()) {
        const Scal ebv = eb.GetVolumeFraction(c);
        if (ebv != 0) {
          a[c] = std::min(a[c], ebv) / ebv;
        } else {
          a[c] = 0;
        }
      }
      af = UEmbed<M>::Interpolate(a, {}, eb).GetFieldFace();
    } else {
      af = Interpolate(a, mf_cond_vfsm_, m);
    }

    const Vect force(var.Vect["force"]);
    const Vect grav(var.Vect["gravity"]);
    const Scal r1(var.Double["rho1"]);
    const Scal r2(var.Double["rho2"]);
    const Scal m1(var.Double["mu1"]);
    const Scal m2(var.Double["mu2"]);

    // Init density and viscosity
    for (auto c : m.AllCells()) {
      const Scal v2 = a[c];
      const Scal v1 = 1. - v2;
      fc_rho_[c] = r1 * v1 + r2 * v2;
      fc_mu_[c] = m1 * v1 + m2 * v2;
    }
    FieldFace<Scal> ff_rho(m);
    for (auto f : m.AllFaces()) {
      const Scal v2 = af[f];
      const Scal v1 = 1. - v2;
      ff_rho[f] = r1 * v1 + r2 * v2;
    }

    if (tracer_ && var.Int["tracer_override_mixture"]) {
      const auto& fc_rho_mix = tracer_->GetMixtureDensity();
      const auto& fc_mu_mix = tracer_->GetMixtureViscosity();
      for (auto c : m.AllCells()) {
        const Scal v2 = a[c];
        const Scal v1 = 1. - v2;
        fc_rho_[c] = fc_rho_mix[c] * v1 + r2 * v2;
        fc_mu_[c] = fc_mu_mix[c] * v1 + m2 * v2;
      }
      FieldFace<Scal> ff_rho_mix(m);
      if (eb_) {
        ff_rho_mix =
            UEmbed<M>::Interpolate(fc_rho_mix, {}, *eb_).GetFieldFace();
      } else {
        ff_rho_mix = Interpolate(fc_rho_mix, mf_cond_vfsm_, m);
      }
      for (auto f : m.AllFaces()) {
        const Scal v2 = af[f];
        const Scal v1 = 1. - v2;
        ff_rho[f] = ff_rho_mix[f] * v1 + r2 * v2;
      }
    }

    // Append gravity to force
    for (auto f : m.AllFaces()) {
      const Vect n = m.GetNormal(f);
      febp_[f] += force.dot(n);
      febp_[f] += grav.dot(n) * ff_rho[f];
    }

    if (eb_) {
      auto& eb = *eb_;
      // Compute gravity as gradient of potential
      // to ensure exact hydrostatic solution.
      MapEmbed<BCond<Scal>> me_neumann;
      mebc_fluid_.LoopBCond(eb, [&](auto cf, IdxCell, auto bc) { //
        me_neumann[cf] = BCond<Scal>(BCondType::neumann, bc.nci);
      });

      auto fe_dens = UEB::InterpolateHarmonic(fc_rho_, me_neumann, eb);
      FieldCell<Scal> fcp(eb);
      for (auto c : eb.AllCells()) {
        const auto x = m.GetCenter(c);
        fcp[c] = grav.dot(x);
      }
      febp_ = UEB::Gradient(fcp, me_neumann, eb);
      eb.LoopFaces([&](auto cf) { //
        febp_[cf] *= fe_dens[cf];
      });
      febp_.SetName(FILELINE + "febp");
      febp_.SetHalo(0);
    }

    // Surface tension
    if (var.Int["enable_surftens"] && as_) {
      CalcSurfaceTension(
          m, layers, var, fc_force_, febp_.GetFieldFace(), fc_sig_,
          GetCondZeroGrad<Scal>(mf_fluid_), fck_, fc_vf0, af, as_.get());
    }

    // vortex force
    Scal force_vort = var.Double["force_vort_g"];
    if (force_vort != 0) {
      Scal r = var.Double["force_vort_r"];
      Vect xc(var.Vect["force_vort_c"]);
      for (auto c : m.Cells()) {
        Vect x = m.GetCenter(c);
        Vect dx = x - xc;
        Scal q = std::exp(-dx.sqrnorm() / sqr(r)) * force_vort;
        Scal fx = -dx[1];
        Scal fy = dx[0];
        fc_force_[c] = Vect(fx, fy, 0.) * q * fc_rho_[c];
      }
    }

    // moving force on the interface
    auto fmov = [this, &a](std::string pre) {
      Vect force_mov(var.Vect[pre]);
      if (force_mov.sqrnorm()) {
        Vect x0(var.Vect[pre + "_x0"]);
        Vect x1(var.Vect[pre + "_x1"]);
        Scal t0 = var.Double[pre + "_t0"];
        Scal t1 = var.Double[pre + "_t1"];
        Vect sig(var.Vect[pre + "_sig"]);
        Scal pi = M_PI;
        int edim = var.Int["dim"];

        Scal t = (st_.t - t0) / (t1 - t0);
        if (t >= 0 && t <= 1) {
          if (var.Int[pre + "_parab"]) {
            t = t * t;
          }
          for (auto c : m.Cells()) {
            Vect xt = x0 + (x1 - x0) * t;
            Vect r = (xt - m.GetCenter(c)) / sig;
            Scal k = std::exp(-r.sqrnorm() * 0.5) /
                     (std::pow(2 * pi, 1. / edim) * sig.prod());
            Scal vf = a[c];
            fc_force_[c] += force_mov * (k * fc_rho_[c] * vf * (1 - a[c]) * 2);
          }
        }
      }
    };

    fmov("force_mov");
    fmov("force_mov2");

    // Kolmogorov forcing
    if (var.Int["force_kolm"]) {
      for (auto f : m.AllFaces()) {
        Vect n = m.GetNormal(f);
        Vect x = m.GetCenter(f);
        febp_[f] += Vect(std::sin(x[1]), 0., 0.).dot(n);
      }
    }

    // Kolmogorov forcing as acceleration
    if (var.Int["force_kolm_accel"]) {
      for (auto f : m.AllFaces()) {
        Vect n = m.GetNormal(f);
        Vect x = m.GetCenter(f);
        febp_[f] += Vect(std::sin(x[1]), 0., 0.).dot(n) * ff_rho[f];
      }
    }

    // zero force in z if 2D
    if (var.Int["dim"] <= 2) {
      for (auto f : m.Faces()) {
        using Dir = typename M::Dir;
        if (m.GetIndexFaces().GetDir(f) == Dir::k) {
          febp_[f] = 0.; // XXX: zero in z
        }
      }
    }
  }
}

template <class M>
void Hydro<M>::DumpFields() {
  auto sem = m.GetSem("dumpfields");
  struct {
    std::array<Multi<FieldCell<Scal>>, dim> im;
    FieldCell<Scal> fc_cellcond;
    FieldCell<Scal> fcdiv; // divergence of velocity
    FieldCell<Scal> fcdis; // energy dissipation
    FieldCell<Scal> fc_ebvf; // embedded boundaries volume fraction
    FieldCell<Scal> fc_tracer_sum; // sum of tracer_ fields starting from 1
  } * ctx(sem);
  if (sem("dump")) {
    if (m.IsRoot()) {
      dumper_.Report();
    }

    auto dl = GetWords(var.String["dumplist"]);

    auto dump = [&dl, this](const FieldCell<Scal>& fc, std::string name) {
      if (dl.count(name)) {
        m.Dump(&fc, name);
      }
    };
    auto dumpv = [&dl, this](
                     const FieldCell<Vect>& fc, size_t i, std::string name) {
      if (dl.count(name)) {
        m.Dump(&fc, i, name);
      }
    };

    auto& fcv = fs_->GetVelocity();
    dumpv(fcv, 0, "vx");
    dumpv(fcv, 1, "vy");
    dumpv(fcv, 2, "vz");
    dump(fs_->GetPressure(), "p");
    dump(as_->GetField(), "vf");
    dump(fc_rho_, "rho");
    dump(fc_mu_, "mu");
    dump(fc_sig_, "sig");
    dump(fc_contang_, "contang");
    dump(fcbc_, "bc");
    if (dl.count("cellcond")) {
      auto& fc = ctx->fc_cellcond;
      fc.Reinit(m, 0);
      for (auto& it : mc_velcond_) {
        fc[it.first] = 1;
      }
      m.Dump(&fc, "cellcond");
    }
    if (var.Int["youngbc"]) {
      dumpv(fcyv_, 0, "yvx");
      dumpv(fcyv_, 1, "yvy");
      dumpv(fcyv_, 2, "yvz");
    }
    if (dl.count("omx") || dl.count("omy") || dl.count("omz") ||
        dl.count("omm") || dl.count("omcalc")) {
      CalcVort();
      dumpv(fcom_, 0, "omx");
      dumpv(fcom_, 1, "omy");
      dumpv(fcom_, 2, "omz");
      dump(fcomm_, "omm");
    }
    if (dl.count("dis") || dl.count("strain")) {
      CalcStrain(fs_->GetVelocity(), fc_strain_);
      if (dl.count("strain")) m.Dump(&fc_strain_, "strain");
      if (dl.count("dis")) {
        ctx->fcdis = fc_strain_;
        for (auto c : m.Cells()) {
          ctx->fcdis[c] *= 2. * fc_mu_[c];
        }
        m.Dump(&ctx->fcdis, "dis");
      }
    }
    if (dl.count("div")) {
      ctx->fcdiv = GetDiv();
      m.Dump(&ctx->fcdiv, "div");
    }
    if (auto as = dynamic_cast<ASV*>(as_.get())) {
      dumpv(as->GetNormal(), 0, "nx");
      dumpv(as->GetNormal(), 1, "ny");
      dumpv(as->GetNormal(), 2, "nz");
      dump(as->GetColor(), "cls");
      dump(fck_[0], "k");
    }
    // TODO reuse ASV code
    if (auto as = dynamic_cast<ASVEB*>(as_.get())) {
      dumpv(as->GetNormal(), 0, "nx");
      dumpv(as->GetNormal(), 1, "ny");
      dumpv(as->GetNormal(), 2, "nz");
      dump(as->GetColor(), "cls");
      dump(fck_[0], "k");
    }
    if (auto as = dynamic_cast<ASVM*>(as_.get())) {
      for (auto l : layers) {
        auto sl = std::to_string(l);
        dump(*as->GetFieldM()[l], "vf" + sl);
        dump(*as->GetColor()[l], "cl" + sl);
        dumpv(*as->GetNormal()[l], 0, "nx" + sl);
        dumpv(*as->GetNormal()[l], 1, "ny" + sl);
        dumpv(*as->GetNormal()[l], 2, "nz" + sl);
        dump(fck_[l], "k" + sl);
      }

      // combined colors
      dump(as->GetColorSum(), "cls");

      // image
      auto conv = [&](size_t d, size_t l,
                      Multi<FieldCell<Scal>>& fc) -> const FieldCell<Scal>& {
        fc.resize(layers);
        fc[l].Reinit(m);
        for (auto c : m.Cells()) {
          fc[l][c] = as->GetImage(l, c)[d];
        }
        return fc[l];
      };
      for (auto d : {0, 1, 2}) {
        for (auto l : layers) {
          std::stringstream st;
          st << "im"
             << "xyz"[d] << l;
          std::string s = st.str();
          dump(conv(d, l, ctx->im[d]), s);
        }
      }
    }
    // TODO add ASVMEB

    if (eb_) {
      auto& eb = *eb_;
      if (dl.count("ebvf")) {
        auto& fc = ctx->fc_ebvf;
        fc.Reinit(m, 0);
        for (auto c : eb.Cells()) {
          fc[c] = eb.GetVolumeFraction(c);
        }
        m.Dump(&fc, "ebvf");
      }
      if (fc_dist_.size()) {
        dump(fc_dist_, "ebdist");
      }
      if (fc_phi_.size()) {
        dump(fc_phi_, "ebphi");
      }
    }

    if (tracer_) {
      for (auto l : tracer_->GetView().layers) {
        dump(tracer_->GetVolumeFraction()[l], "tu" + std::to_string(l));
      }
      if (dl.count("tusum")) {
        ctx->fc_tracer_sum.Reinit(m, 0);
        for (auto l : tracer_->GetView().layers) {
          if (l > 0) {
            const auto& fc = tracer_->GetVolumeFraction()[l];
            for (auto c : m.Cells()) {
              ctx->fc_tracer_sum[c] += fc[c];
            }
          }
        }
        dump(ctx->fc_tracer_sum, "tusum");
      }
    }
  }
  if (var.Int["enable_advection"]) {
    if (var.Int["dumppoly"] && sem.Nested()) {
      as_->DumpInterface(GetDumpName("s", ".vtk", dumper_.GetN()));
    }
    if (var.Int["dumppolymarch"] && sem.Nested()) {
      as_->DumpInterfaceMarch(GetDumpName("sm", ".vtk", dumper_.GetN()));
    }
  }
  if (particles_ && var.Int["dump_particles"]) {
    const std::string path = GetDumpName("part", ".csv", dumper_.GetN());
    if (sem()) {
      if (m.IsRoot()) {
        std::cout << std::fixed << std::setprecision(8) << "dump"
                  << " t=" << particles_->GetTime() << " to " << path
                  << std::endl;
      }
    }
    if (sem.Nested()) {
      particles_->DumpCsv(path);
    }
  }
  if (sem()) {
    // XXX: empty stage, otherwise ctx is destroyed before dump
  }
}

template <class M>
void Hydro<M>::Dump(bool force) {
  auto sem = m.GetSem("dump");
  struct {
    Multi<FieldCell<MIdx>> fcim;
  } * ctx(sem);
  if (sem.Nested("fields")) {
    if (dumper_.Try(st_.t, st_.dt) || force) {
      DumpFields();
    }
  }
  if (dmptraj_.Try(st_.t, st_.dt) || force) {
    if (sem("copyimage")) {
      ctx->fcim.resize(layers);
      ctx->fcim.InitAll(FieldCell<MIdx>(m));
      if (auto as = dynamic_cast<ASVM*>(as_.get())) {
        for (auto c : m.AllCells()) {
          for (auto l : layers) {
            ctx->fcim[l][c] = as->GetImage(l, c);
          }
        }
      } else if (auto as = dynamic_cast<ASVMEB*>(as_.get())) {
        for (auto c : m.AllCells()) {
          for (auto l : layers) {
            ctx->fcim[l][c] = as->GetImage(l, c);
          }
        }
      } else if (auto as = dynamic_cast<ASV*>(as_.get())) {
        for (auto c : m.AllCells()) {
          ctx->fcim[0][c] = as->GetImage(c);
        }
      } else if (auto as = dynamic_cast<ASVEB*>(as_.get())) {
        for (auto c : m.AllCells()) {
          ctx->fcim[0][c] = as->GetImage(c);
        }
      }
    }
    if (sem.Nested("trajdump")) {
      if (var.Int["enable_color"]) {
        Multi<const FieldCell<Scal>*> fcu(layers);
        Multi<const FieldCell<Scal>*> fccl(layers);
        if (auto as = dynamic_cast<ASVM*>(as_.get())) {
          fcu = as->GetFieldM();
          fccl = as->GetColor();
        } else if (auto as = dynamic_cast<ASVMEB*>(as_.get())) {
          // TODO reuse ASVM code
          fcu = as->GetFieldM();
          fccl = as->GetColor();
        } else if (auto as = dynamic_cast<ASVEB*>(as_.get())) {
          // TODO reuse ASV code
          fcu[0] = &as->GetField();
          fccl[0] = &as->GetColor();
        } else if (auto as = dynamic_cast<ASV*>(as_.get())) {
          fcu[0] = &as->GetField();
          fccl[0] = &as->GetColor();
        }
        if (eb_) {
          DumpTraj<EB>(
              *eb_, true, var, dmptraj_.GetN(), st_.t, layers, fcu, fccl,
              ctx->fcim, fs_->GetPressure(), fs_->GetVelocity(), fcvm_, st_.dt);
        } else {
          DumpTraj<M>(
              m, true, var, dmptraj_.GetN(), st_.t, layers, fcu, fccl,
              ctx->fcim, fs_->GetPressure(), fs_->GetVelocity(), fcvm_, st_.dt);
        }
      }
    }
  }
  if (sem("dmptrep")) {
    if (m.IsRoot() && dmptrep_.Try(st_.t, st_.dt)) {
      std::string s = GetDumpName("trep", ".log", dmptrep_.GetN());
      m.TimerReport(s);
      std::cout << std::fixed << std::setprecision(8) << "timer report"
                << " t=" << st_.t << " to " << s << std::endl;
    }
  }
  if (sem("dumpstat")) {
    if (m.IsRoot()) {
      if (st_.step % var.Int("stat_step_every", 1) == 0 || force) {
        ost_->Write(0., "");
      }
    }
  }
  if (auto as = dynamic_cast<ASV*>(as_.get())) {
    if (psm_ && dumper_.Try(st_.t, st_.dt)) {
      if (var.Int["dumppart"] && sem.Nested("part-dump")) {
        psm_->DumpParticles(
            &as->GetAlpha(), &as->GetNormal(), dumper_.GetN(), st_.t);
      }
      if (var.Int["dumppartinter"] && sem.Nested("partinter-dump")) {
        psm_->DumpPartInter(
            &as->GetAlpha(), &as->GetNormal(), dumper_.GetN(), st_.t);
      }
    }
  }
  // TODO reuse ASV code
  if (auto as = dynamic_cast<ASVEB*>(as_.get())) {
    if (psm_ && dumper_.Try(st_.t, st_.dt)) {
      if (var.Int["dumppart"] && sem.Nested("part-dump")) {
        psm_->DumpParticles(
            &as->GetAlpha(), &as->GetNormal(), dumper_.GetN(), st_.t);
      }
      if (var.Int["dumppartinter"] && sem.Nested("partinter-dump")) {
        psm_->DumpPartInter(
            &as->GetAlpha(), &as->GetNormal(), dumper_.GetN(), st_.t);
      }
    }
  }
  if (auto as = dynamic_cast<ASVM*>(as_.get())) {
    if (psm_ && dumper_.Try(st_.t, st_.dt)) {
      if (var.Int["dumppart"] && sem.Nested("part-dump")) {
        psm_->DumpParticles(
            as->GetAlpha(), as->GetNormal(), dumper_.GetN(), st_.t);
      }
      if (var.Int["dumppartinter"] && sem.Nested("partinter-dump")) {
        psm_->DumpPartInter(
            as->GetAlpha(), as->GetNormal(), dumper_.GetN(), st_.t);
      }
    }
  }
  // TODO reuse ASVM code
  if (auto as = dynamic_cast<ASVMEB*>(as_.get())) {
    if (psm_ && dumper_.Try(st_.t, st_.dt)) {
      if (var.Int["dumppart"] && sem.Nested("part-dump")) {
        psm_->DumpParticles(
            as->GetAlpha(), as->GetNormal(), dumper_.GetN(), st_.t);
      }
      if (var.Int["dumppartinter"] && sem.Nested("partinter-dump")) {
        psm_->DumpPartInter(
            as->GetAlpha(), as->GetNormal(), dumper_.GetN(), st_.t);
      }
    }
  }
}

template <class M>
void Hydro<M>::Run() {
  auto sem = m.GetSem("run");
  struct {
    Scal nabort;
  } * ctx(sem);

  if (sem.Nested("init")) {
    Init();
  }

  sem.LoopBegin();

  if (sem("events")) {
    if (events_) {
      events_->Exec(st_.t);
    }
  }
  if (sem("loop-check")) {
    if (st_.t + st_.dt * 0.25 > var.Double["tmax"]) {
      if (m.IsRoot()) {
        std::cout << "End of simulation, t > tmax=" << var.Double["tmax"]
                  << std::endl;
      }
      sem.LoopBreak();
    } else if (int(st_.step + 0.5) >= var.Int["max_step"]) {
      if (m.IsRoot()) {
        std::cout << "End of simulation, step > max_step="
                  << var.Int["max_step"] << std::endl;
      }
      sem.LoopBreak();
    } else if (st_.step > 1 && fs_->GetError() < var.Double("stop_diff", 0)) {
      if (m.IsRoot()) {
        std::cout << "End of simulation, diff < stop_diff="
                  << var.Double["stop_diff"] << std::endl;
      }
      sem.LoopBreak();
    } else {
      if (m.IsRoot()) {
        if (st_.step % var.Int("report_step_every", 1) == 0) {
          ReportStep();
        }
      }
      m.SeedSample();
    }
  }

  CheckAbort(sem, ctx->nabort);

  if (sem("updatepar")) {
    if (auto fs = dynamic_cast<Simple<M>*>(fs_.get())) {
      fs->SetPar(ParsePar<Simple<M>>()(var));
    } else if (auto fs = dynamic_cast<Proj<M>*>(fs_.get())) {
      fs->SetPar(ParsePar<Proj<M>>()(var));
    }
    UpdateAdvectionPar();
    fcvm_ = fs_->GetVelocity();
  }
  if (sem.Nested("mixture")) {
    CalcMixture(as_->GetField());
  }
  if (sem.Nested("stat")) {
    CalcStat();
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

  if (sem.Nested()) {
    CalcDt(); // must be after CalcStat to keep dt for moving mesh velocity
  }

  if (sem.Nested()) {
    Dump(false);
  }

  if (sem.Nested("stephook")) {
    StepHook(this);
  }

  if (sem("inc")) {
    ++st_.step;
    m.CollectSample("Hydro::Step");
  }

  sem.LoopEnd();

  if (sem.Nested("dumplast")) {
    if (var.Int["dumplast"]) {
      Dump(true);
    }
  }

  if (sem.Nested("posthook")) {
    if (eb_) {
      PostHook(var, fs_->GetVelocity(), m, *eb_);
    } else {
      PostHook(var, fs_->GetVelocity(), m);
    }
  }
}

template <class M>
void Hydro<M>::ReportStep() {
  std::cout << std::fixed << std::setprecision(8) << "STEP=" << st_.step
            << " t=" << st_.t << " dt=" << st_.dt << " ta=" << as_->GetTime()
            << " dta=" << as_->GetTimeStep() << " wt=" << timer_.GetSeconds()
            << std::endl;
}

template <class M>
void Hydro<M>::ReportStepAdv() {
  std::cout << std::fixed << std::setprecision(8)
            << ".....adv: t=" << as_->GetTime() << " dt=" << as_->GetTimeStep()
            << std::endl;
}

template <class M>
void Hydro<M>::ReportStepTracer() {
  std::cout << std::fixed << std::setprecision(8)
            << ".....tracer: t=" << tracer_->GetTime() << " dt=" << tracer_dt_
            << std::endl;
}

template <class M>
void Hydro<M>::ReportStepParticles() {
  std::cout << std::fixed << std::setprecision(8)
            << ".....particles: t=" << particles_->GetTime()
            << " dt=" << particles_dt_ << std::endl;
}

template <class M>
auto Hydro<M>::CalcPressureDrag(const FieldCell<Scal>& fcp, const Embed<M>& eb)
    -> Vect {
  MapEmbed<BCond<Scal>> me_pressure;
  auto& m = eb.GetMesh();
  mebc_fluid_.LoopBCond(eb, [&](auto cf, IdxCell, auto bc) { //
    const auto nci = bc.nci;
    if (bc.type == BCondFluidType::slipwall ||
        bc.type == BCondFluidType::symm) {
      me_pressure[cf] = BCond<Scal>(BCondType::neumann, nci);
    } else {
      me_pressure[cf] = BCond<Scal>(BCondType::extrap, nci);
    }
  });
  auto fep = UEmbed<M>::Interpolate(fcp, me_pressure, eb);
  Vect sum(0);
  mebc_fluid_.LoopBCond(eb, [&](auto cf, IdxCell c, auto bc) { //
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
  auto& m = eb.GetMesh();
  MapEmbed<BCond<Scal>> me_neumann;
  mebc_fluid_.LoopBCond(eb, [&](auto cf, IdxCell, auto bc) { //
    me_neumann[cf] = BCond<Scal>(BCondType::neumann, bc.nci);
  });
  auto feg = UEmbed<M>::Gradient(fcvel, fs_->GetVelocityCond(), eb);
  auto femu = UEmbed<M>::Interpolate(fcmu, me_neumann, eb);
  Vect sum(0);
  mebc_fluid_.LoopBCond(eb, [&](auto cf, IdxCell c, auto bc) { //
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
  std::cout << std::scientific << std::setprecision(16)
            << ".....iter=" << fs_->GetIter() << ", diff=" << fs_->GetError()
            << std::endl;
}

template <class M>
void Hydro<M>::CheckAbort(Sem& sem, Scal& nabort) {
  if (sem("abort-local")) {
    nabort = 0.;
    try {
      CHECKNAN(as_->GetField(), true)
      CHECKNAN(fs_->GetVelocity(), true)
      CHECKNAN(fs_->GetPressure(), true)
      // check abort TODO: revise,move
      for (auto c : m.Cells()) {
        if (fs_->GetVelocity()[c].sqrnorm() > sqr(var.Double["abortvel"])) {
          std::stringstream g;
          g << "abortvel exceeded at x=" << m.GetCenter(c);
          throw std::runtime_error(g.str());
        }
      }
    } catch (const std::runtime_error& e) {
      std::cout << e.what() << std::endl;
      nabort += 1.;
    }
    m.Reduce(&nabort, "sum");
  }

  if (sem("abort-reduce")) {
    if (nabort != 0.) {
      if (m.IsRoot()) {
        std::cout << "nabort = " << nabort << std::endl;
      }
      sem.LoopBreak();
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
      if (st_.step % var.Int("report_step_every", 1) == 0) {
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
      if (st_.step % var.Int("report_step_every", 1) == 0) {
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
  if (sem.Nested("start")) {
    particles_->Step(particles_dt_, fs_->GetVolumeFlux());
  }
  if (sem("spawn")) {
    std::vector<Vect> p_x;
    std::vector<Vect> p_v;
    std::vector<Scal> p_r;
    std::vector<Scal> p_rho;
    std::vector<Scal> p_termvel;
    ParticlesView view{p_x, p_v, p_r, p_rho, p_termvel};
    SpawnParticles(view);
    particles_->Append(view);
  }
  if (sem("report")) {
    if (m.IsRoot()) {
      if (st_.step % var.Int("report_step_every", 1) == 0) {
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
      if (st_.step % var.Int("report_step_every", 1) == 0) {
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
    psm_ = UCurv<M>::CalcCurvPart(as_.get(), psm_par_, fck_, m);
  }
  if (var.Int["enable_bubgen"]) {
    if (sem.Nested("bubgen")) {
      StepBubgen();
    }
  }
}

template <class M>
void Hydro<M>::InitTracerFields(Multi<FieldCell<Scal>>& vfcu) {
  auto sem = m.GetSem("tracerfields");
  struct {
    Vars vart;
  } * ctx(sem);
  if (sem("init")) {
    vfcu.resize(var.Int["tracer_layers"]);
    vfcu.InitAll(FieldCell<Scal>(m, 0));
  }
  for (int l = 0; l < var.Int["tracer_layers"]; ++l) {
    const std::string prefix = "tracer" + std::to_string(l);
    if (sem("var" + std::to_string(l))) {
      ctx->vart.String.Set("init_vf", var.String[prefix + "_init"]);
      ctx->vart.String.Set("list_path", var.String[prefix + "_list_path"]);
      ctx->vart.Int.Set("dim", var.Int["dim"]);
      ctx->vart.Int.Set("list_ls", 3);
    }
    if (sem.Nested("field" + std::to_string(l))) {
      InitVf(vfcu[l], ctx->vart, m);
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
    FieldCell<Scal> fcvf; // volume fraction
    Vars var;
  } * ctx(sem);
  auto& fcvf = ctx->fcvf;
  const Scal t0 = var.Double["bubgen_t0"];
  const Scal tper = var.Double["bubgen_per"];
  bool bg = (st_.t > t0 && st_.t - bgt_ >= tper);
  if (bg) {
    if (sem("as-bubgen-var")) {
      ctx->var.String.Set("init_vf", "list");
      ctx->var.String.Set("list_path", var.String["bubgen_path"]);
      ctx->var.Int.Set("dim", var.Int["dim"]);
      ctx->var.Int.Set("list_ls", var.Int["list_ls"]);
    }
    if (sem.Nested("as-bubgen-initvf")) {
      InitVf(fcvf, ctx->var, m);
    }
    if (sem("as-bubgen-apply")) {
      bgt_ = st_.t;
      auto apply_vof = [&](auto* as, const auto& eb) {
        if (as) {
          auto& u = const_cast<FieldCell<Scal>&>(as->GetField());
          for (auto c : eb.AllCells()) {
            if (fcvf[c] > 0) {
              u[c] = std::max(u[c], fcvf[c]);
            }
          }
        }
      };
      auto apply_vofm = [&](auto* as, const auto& eb) {
        if (as) {
          auto& u = const_cast<FieldCell<Scal>&>(*as->GetFieldM()[0]);
          auto& cl = const_cast<FieldCell<Scal>&>(*as->GetColor()[0]);
          for (auto c : eb.AllCells()) {
            if (fcvf[c] > 0) {
              u[c] = std::max(u[c], fcvf[c]);
              cl[c] = 1.;
            }
          }
        }
      };
      apply_vofm(dynamic_cast<ASVM*>(as_.get()), m);
      apply_vof(dynamic_cast<ASV*>(as_.get()), m);
      if (eb_) {
        auto& eb = *eb_;
        for (auto c : m.AllCells()) {
          fcvf[c] = std::min(fcvf[c], eb.GetVolumeFraction(c));
        }
        apply_vofm(dynamic_cast<ASVMEB*>(as_.get()), *eb_);
        apply_vof(dynamic_cast<ASVEB*>(as_.get()), *eb_);
      }
    }
    if (sem()) {
      // FIXME: empty stage to finish communication to keep ctx
    }
  }
}
