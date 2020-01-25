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
#include "geom/mesh.h"
#include "kernelmesh.h"
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
#include "solver/curv.h"
#include "solver/embed.h"
#include "solver/multi.h"
#include "solver/normal.h"
#include "solver/proj.h"
#include "solver/proj_eb.h"
#include "solver/reconst.h"
#include "solver/simple.h"
#include "solver/solver.h"
#include "solver/tvd.h"
#include "solver/vof.h"
#include "solver/vof_eb.h"
#include "solver/vofm.h"
#include "util/convdiff.h"
#include "util/events.h"
#include "util/hydro.h"
#include "util/metrics.h"
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
  static constexpr size_t dim = M::dim;

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
  void InitFluid(const FieldCell<Vect>& fc_vel);
  void InitAdvection(const FieldCell<Scal>& fcvf, const FieldCell<Scal>& fccl);
  void InitStat();
  static Vect CalcPressureDrag(const FieldCell<Scal>& fcp, const Embed<M>& eb);
  static Vect CalcViscousDrag(
      const FieldCell<Vect>& fcvel, const FieldCell<Scal>& fcmu,
      const Embed<M>& eb);
  void DumpFields();
  void Dump();
  // Calc rho, mu and force based on volume fraction
  void CalcMixture(const FieldCell<Scal>& vf);
  // Clips v to given range, uses const_cast
  void Clip(const FieldCell<Scal>& v, Scal min, Scal max);
  void CalcStat();
  void CalcDt();
  void ReportStep();
  void ReportStepAdv();
  void ReportIter();
  // Issue sem.LoopBreak if abort conditions met
  void CheckAbort(Sem& sem, Scal& nabort);
  void StepFluid();
  void StepAdvection();
  void StepBubgen();

  using AST = Tvd<M>; // advection TVD
  using ASV = Vof<M>; // advection VOF
  using ASVM = Vofm<M>; // advection VOF
  using ASVEB = VofEmbed<M>; // advection VOF embed
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
    fcom_ = GetVort(fcv, fs_->GetVelocityCond(), m);
    fcomm_.Reinit(m);
    for (auto c : m.Cells()) {
      fcomm_[c] = fcom_[c].norm();
    }
  }
  void CalcStrain(const FieldCell<Vect> fcvel, FieldCell<Scal>& fc_strain) {
    auto& fcv = fcvel;
    auto& fcs = fc_strain;

    auto ffv = Interpolate(fcv, fs_->GetVelocityCond(), m);

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
    auto& ffv = fs_->GetVolumeFlux();
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
  FieldFace<Scal> ffbp_; // balanced force projections

  MapCondFaceAdvection<Scal> mf_adv_;
  MapCondFace mf_cond_vfsm_;
  MapCondFaceFluid mf_fluid_; // fluid cond
  MapCell<std::shared_ptr<CondCell>> mc_cond_;
  MapCell<std::shared_ptr<CondCellFluid>> mc_velcond_;

  std::unique_ptr<Embed<M>> eb_;
  std::unique_ptr<AdvectionSolver<M>> as_;
  std::unique_ptr<FluidSolver<M>> fs_;

  FieldCell<Scal> fc_smvf_; // smoothed volume fraction used by CalcMixture()
  FieldFace<Scal> ff_smvf_; // interpolated fc_smvf_
  Scal diff_; // convergence indicator

  FieldCell<Scal> fc_sig_; // surface tension sigma

  FieldCell<Vect> fcvm_; // velocity field time_prev // TODO: revise

  FieldCell<Vect> fcyv_; // Young velocity
  FieldCell<Vect> fcom_; // vorticity
  FieldCell<Scal> fcomm_; // vorticity magnitude
  FieldCell<Scal> fcbc_; // boundary condition type by GetBcField()
  FieldCell<Scal> fc_strain_; // double inner product of strain rate tensor
  Multi<FieldCell<Scal>> fck_; // curvature
  typename PartStrMeshM<M>::Par psm_par_;
  std::unique_ptr<PartStrMeshM<M>> psm_;

  std::function<void(FieldCell<typename M::Scal>&, const M&)> bgf_;
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
      m.Reduce(&(mst[name + "_mn"] = min), "min");
      m.Reduce(&(mst[name + "_mx"] = max), "max");
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
      m.Reduce(&(mst[name + "_mn"] = min), "min");
      m.Reduce(&(mst[name + "_mx"] = max), "max");
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
};

template <class M>
void Hydro<M>::InitEmbed() {
  if (var.Int["enable_embed"]) {
    auto sem = m.GetSem("embed");
    struct {
      FieldNode<Scal> fnl;
    } * ctx(sem);
    if (sem("ctor")) {
      eb_.reset(new Embed<M>(m));
      ctx->fnl = ::InitEmbed(m, var, m.IsRoot());
    }
    if (sem.Nested("init")) {
      eb_->Init(ctx->fnl);
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
  if (eb_) {
    auto p = ParsePar<Proj<M>>()(var);
    fs_.reset(new ProjEmbed<M>(
        m, fc_vel, mf_fluid_, *eb_, 0, Vect(0), mc_velcond_, &fc_rho_, &fc_mu_,
        &fc_force_, &ffbp_, &fc_src_, &fc_srcm_, 0., st_.dt, p));
  } else if (fs == "simple") {
    auto p = ParsePar<Simple<M>>()(var);
    fs_.reset(new Simple<M>(
        m, fc_vel, mf_fluid_, mc_velcond_, &fc_rho_, &fc_mu_, &fc_force_,
        &ffbp_, &fc_src_, &fc_srcm_, 0., st_.dt, p));
  } else if (fs == "proj") {
    auto p = ParsePar<Proj<M>>()(var);
    fs_.reset(new Proj<M>(
        m, fc_vel, mf_fluid_, mc_velcond_, &fc_rho_, &fc_mu_, &fc_force_,
        &ffbp_, &fc_src_, &fc_srcm_, 0., st_.dt, p));
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
          m, fcvf, fccl, mf_adv_, &fs_->GetVolumeFlux(Step::time_curr),
          &fc_src2_, 0., st_.dta, p));
    }
    layers = GRange<size_t>(1);
  } else if (as == "vofm") {
    auto p = ParsePar<ASVM>()(var);
    auto as = new ASVM(
        m, fcvf, fccl, mf_adv_, &fs_->GetVolumeFlux(Step::time_curr), &fc_src2_,
        0., st_.dta, p);
    as_.reset(as);
    layers = GRange<size_t>(as->GetNumLayers());
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
        op("diff", &diff_),
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
    Vars varbody;
  } * ctx(sem);
  auto& fcvel = ctx->fcvel;
  auto& fcvf = ctx->fcvf;
  auto& fccl = ctx->fccl;
  if (sem.Nested()) {
    InitVf(fcvf, var, m);
  }
  if (sem("fields")) {
    fc_src_.Reinit(m, 0.);
    fc_src2_.Reinit(m, 0.);
    fc_srcm_.Reinit(m, 0.);

    // initial surface tension sigma
    fc_sig_.Reinit(m, 0);
    auto isig = CreateInitSig<M>(var);
    isig(fc_sig_, m);
    m.Comm(&fc_sig_);

    // initial velocity
    fcvel.Reinit(m, Vect(0));
    InitVel(fcvel, var, m);
    m.Comm(&fcvel);

    // global mesh size
    MIdx gs = m.GetGlobalSize();

    if (m.IsRoot()) {
      std::cout << "global mesh=" << gs << std::endl;
      std::cout << "surface tension dt=" << GetStDt() << std::endl;
      std::cout << "viscosity dt=" << GetVisDt() << std::endl;
    }

    // boundary conditions
    GetFluidFaceCond(var, m, mf_fluid_, mf_adv_);

    // boundary conditions for smoothing of volume fraction
    for (auto& it : mf_fluid_) {
      const IdxFace f = it.first;
      auto& cb = it.second;
      size_t nci = cb->GetNci();
      if (cb.template Get<Symm<M>>()) {
        mf_cond_vfsm_[f].template Set<CondFaceReflect>(nci);
      } else {
        mf_cond_vfsm_[f].template Set<CondFaceGradFixed<Scal>>(Scal(0), nci);
      }
    }
  }

  if (var.Int["bc_wall_init_vel"] && sem("bc_wall_init_vel")) {
    // velocity on walls from neighbour cells
    for (auto& it : mf_fluid_) {
      const IdxFace f = it.first;
      auto& cb = it.second;
      if (auto cd = cb.template Get<NoSlipWallFixed<M>>()) {
        IdxCell c = m.GetCell(f, cd->GetNci());
        cd->SetVelocity(fcvel[c]);
      }
    }
  }

  if (sem.Nested("initvort")) {
    if (var.Int["initvort"]) {
      InitVort(fcvel, fcvel, mf_fluid_, var.Int["vort_report"], m);
    }
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
    varbody.Int.Set("list_ls", 0);
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
    AppendBodyCond<M>(
        fc, var.String["body_bc"], m, clear0, clear1, inletcl, fill_vf, nullptr,
        mf_fluid_, mf_adv_);
  }

  if (sem("dt")) {
    const Scal dt = var.Double["dt0"];
    st_.dt = dt;
    st_.dta = dt;
    if (m.IsLead()) {
      this->var_mutable.Double.Set("dt", st_.dt);
      this->var_mutable.Double.Set("dta", st_.dta);
    }
  }
  if (sem.Nested("embed")) {
    InitEmbed();
  }

  if (sem("solv")) {
    InitFluid(fcvel);

    InitAdvection(fcvf, fccl);

    st_.iter = 0;
    st_.t = fs_->GetTime();

    if (m.IsLead()) {
      this->var_mutable.Int.Set("iter", st_.iter);
      this->var_mutable.Double.Set("t", st_.t);
    }

    InitStat();

    if (var.Int["fill_halo_nan"]) {
      std::vector<std::pair<IdxFace, size_t>> vf;
      for (auto& p : mf_adv_) {
        vf.emplace_back(p.first, p.second.GetNci());
      }
      m.SetNanFaces(vf);
    }

    if (m.IsLead()) {
      events_ = std::unique_ptr<Events>(
          new Events(this->var_mutable, m.IsRoot(), m.IsLead()));
      events_->Parse();
    }
  }
  if (var.Int["dumpbc"] && sem.Nested()) {
    DumpBcFaces(mf_adv_, mf_fluid_, "bc.vtk", m);
  }
  if (eb_ && sem.Nested()) {
    eb_->DumpPoly(var.Int["vtkbin"], var.Int["vtkmerge"]);
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
      for (size_t d = 0; d < dim; ++d) {
        m.Reduce(&st_.eb_pdrag[d], "sum");
        m.Reduce(&st_.eb_vdrag[d], "sum");
        m.Reduce(&st_.eb_drag[d], "sum");
      }
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
      for (auto& it : mf_fluid_) {
        const IdxFace f = it.first;
        CondFaceFluid* cb = it.second.Get();
        if (auto cd = dynamic_cast<fluid_condition::NoSlipWallFixed<M>*>(cb)) {
          size_t nci = cd->GetNci();
          IdxCell c = m.GetCell(f, nci);
          cd->SetVelocity(slipvel * kslip * fa[c]);
        }
      }
    }

    // XXX: slip velocity penalization
    const Scal penalslip = var.Double["penalslip"];
    if (penalslip != 0) {
      Scal dt = fs_->GetTimeStep();
      Vect slipvel(var.Vect["slipvel"]);
      // const auto& fa = as_->GetField();
      const auto& fa = fc_smvf_;
      const auto& fv = fs_->GetVelocity();
      for (auto& it : mf_fluid_) {
        const IdxFace f = it.first;
        CondFaceFluid* cb = it.second.Get();
        if (auto cd = dynamic_cast<fluid_condition::NoSlipWallFixed<M>*>(cb)) {
          size_t nci = cd->GetNci();
          IdxCell c = m.GetCell(f, nci);
          Scal sgn = (slipvel - fv[c]).dot(slipvel);
          if (sgn > 0) {
            fc_force_[c] +=
                (slipvel - fv[c]) * (fc_rho_[c] * penalslip * fa[c] / dt);
          }
        }
      }
    }
    const Scal slipnormal = var.Double["slipnormal"];
    if (slipnormal != 0) {
      const auto& fa = fc_smvf_;
      for (auto& it : mf_fluid_) {
        const IdxFace f = it.first;
        CondFaceFluid* cb = it.second.Get();
        if (auto cd = dynamic_cast<fluid_condition::NoSlipWallFixed<M>*>(cb)) {
          size_t nci = cd->GetNci();
          IdxCell c = m.GetCell(f, nci);
          fc_force_[c] += m.GetNormal(f) * ((nci == 1 ? 1 : -1) * fc_rho_[c] *
                                            slipnormal * fa[c]);
        }
      }
    }
  }

  if (sem("young")) {
    if (var.Int["youngbc"]) {
      InitYoung();
      for (auto& it : mf_fluid_) {
        const IdxFace f = it.first;
        Vect x = m.GetCenter(f);
        CondFaceFluid* cb = it.second.Get();
        if (auto cd = dynamic_cast<fluid_condition::NoSlipWallFixed<M>*>(cb)) {
          cd->SetVelocity(GetYoungVel(x));
        }
      }
    }
  }
}

template <class M>
void Hydro<M>::CalcDt() {
  auto sem = m.GetSem("dt");

  if (sem("local")) {
    st_.t = fs_->GetTime();
    if (m.IsLead()) {
      this->var_mutable.Double["t"] = st_.t;
    }

    st_.dtt = fs_->GetAutoTimeStep();
    m.Reduce(&st_.dtt, "min");
  }
  if (sem("reduce")) {
    // set from cfl if defined
    if (auto* cfl = var.Double.Find("cfl")) {
      st_.dt = st_.dtt * (*cfl);
      st_.dt = std::min<Scal>(st_.dt, var.Double["dtmax"]);
    }

    // constraint from surface tension
    st_.dt = std::min<Scal>(st_.dt, GetStDt());

    // constraint from viscosity
    st_.dt = std::min<Scal>(st_.dt, GetVisDt());

    fs_->SetTimeStep(st_.dt);
    if (m.IsLead()) {
      this->var_mutable.Double["dt"] = st_.dt;
    }

    // set from cfla if defined
    if (auto* cfla = var.Double.Find("cfla")) {
      st_.dta = st_.dtt * (*cfla);
      st_.dta = std::min<Scal>(st_.dta, var.Double["dtmax"]);
    }

    // round up dta to such that dt / dta is integer
    Scal dt = fs_->GetTime() + fs_->GetTimeStep() - as_->GetTime();
    st_.dta = dt / std::max(1, int(dt / st_.dta + 0.5));

    as_->SetTimeStep(st_.dta);
    if (m.IsLead()) {
      this->var_mutable.Double["dta"] = st_.dta;
    }
  }
}

template <class M>
void Hydro<M>::CalcMixture(const FieldCell<Scal>& fc_vf0) {
  auto sem = m.GetSem("mixture");

  if (sem("init")) {
    fc_mu_.Reinit(m);
    fc_rho_.Reinit(m);
    fc_force_.Reinit(m, Vect(0));
    ffbp_.Reinit(m, 0);
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
    af = Interpolate(a, mf_cond_vfsm_, m);

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
    // Init density and viscosity
    for (auto f : m.AllFaces()) {
      const Scal v2 = af[f];
      const Scal v1 = 1. - v2;
      ff_rho[f] = r1 * v1 + r2 * v2;
    }

    // Append gravity to force
    for (auto f : m.AllFaces()) {
      Vect n = m.GetNormal(f);
      ffbp_[f] += force.dot(n);
      ffbp_[f] += grav.dot(n) * ff_rho[f];
    }

    // Surface tension
    if (var.Int["enable_surftens"] && as_) {
      CalcSurfaceTension(
          m, layers, var, fc_force_, ffbp_, fc_sig_,
          GetCondZeroGrad<Scal>(mf_fluid_), fck_, fc_vf0, af, as_.get());
    }

    if (0&&eb_) {
      auto& eb = *eb_;
      for (auto f : m.AllFaces()) {
        if (eb.GetType(f) == Embed<M>::Type::excluded) {
          ffbp_[f] = 0;
        }
      }
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
        ffbp_[f] += Vect(std::sin(x[1]), 0., 0.).dot(n);
      }
    }

    // Kolmogorov forcing as acceleration
    if (var.Int["force_kolm_accel"]) {
      for (auto f : m.AllFaces()) {
        Vect n = m.GetNormal(f);
        Vect x = m.GetCenter(f);
        ffbp_[f] += Vect(std::sin(x[1]), 0., 0.).dot(n) * ff_rho[f];
      }
    }

    // zero force in z if 2D
    if (var.Int["dim"] <= 2) {
      for (auto f : m.Faces()) {
        using Dir = typename M::Dir;
        if (m.GetIndexFaces().GetDir(f) == Dir::k) {
          ffbp_[f] = 0.; // XXX: zero in z
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
  } * ctx(sem);
  if (sem("dump")) {
    if (m.IsRoot()) {
      dumper_.Report();
    }

    auto dl = GetWords(var.String["dumplist"]);
    auto& fcv = fs_->GetVelocity();
    if (dl.count("vx")) m.Dump(&fcv, 0, "vx");
    if (dl.count("vy")) m.Dump(&fcv, 1, "vy");
    if (dl.count("vz")) m.Dump(&fcv, 2, "vz");
    if (dl.count("p")) m.Dump(&fs_->GetPressure(), "p");
    if (dl.count("vf")) m.Dump(&as_->GetField(), "vf");
    if (dl.count("rho")) m.Dump(&fc_rho_, "rho");
    if (dl.count("mu")) m.Dump(&fc_mu_, "mu");
    if (dl.count("sig")) m.Dump(&fc_sig_, "sig");
    if (dl.count("bc")) m.Dump(&fcbc_, "bc");
    if (dl.count("cellcond")) {
      auto& fc = ctx->fc_cellcond;
      fc.Reinit(m, 0);
      for (auto& it : mc_velcond_) {
        fc[it.first] = 1;
      }
      m.Dump(&fc, "cellcond");
    }
    if (var.Int["youngbc"]) {
      if (dl.count("yvx")) m.Dump(&fcyv_, 0, "yvx");
      if (dl.count("yvy")) m.Dump(&fcyv_, 1, "yvy");
      if (dl.count("yvz")) m.Dump(&fcyv_, 2, "yvz");
    }
    if (dl.count("omx") || dl.count("omy") || dl.count("omz") ||
        dl.count("omm") || dl.count("omcalc")) {
      CalcVort();
      if (dl.count("omx")) m.Dump(&fcom_, 0, "omx");
      if (dl.count("omy")) m.Dump(&fcom_, 1, "omy");
      if (dl.count("omz")) m.Dump(&fcom_, 2, "omz");
      if (dl.count("omm")) m.Dump(&fcomm_, "omm");
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
      if (dl.count("nx")) m.Dump(&as->GetNormal(), 0, "nx");
      if (dl.count("ny")) m.Dump(&as->GetNormal(), 1, "ny");
      if (dl.count("nz")) m.Dump(&as->GetNormal(), 2, "nz");
      if (dl.count("cls")) m.Dump(&as->GetColor(), "cls");
      if (dl.count("k")) m.Dump(&fck_[0], "k");
    }
    // TODO reuse ASV code
    if (auto as = dynamic_cast<ASVEB*>(as_.get())) {
      if (dl.count("nx")) m.Dump(&as->GetNormal(), 0, "nx");
      if (dl.count("ny")) m.Dump(&as->GetNormal(), 1, "ny");
      if (dl.count("nz")) m.Dump(&as->GetNormal(), 2, "nz");
      if (dl.count("cls")) m.Dump(&as->GetColor(), "cls");
      if (dl.count("k")) m.Dump(&fck_[0], "k");
    }
    if (auto as = dynamic_cast<ASVM*>(as_.get())) {
      for (auto l : layers) {
        auto sl = std::to_string(l);
        if (dl.count("vf" + sl)) m.Dump(as->GetFieldM()[l], "vf" + sl);
        if (dl.count("cl" + sl)) m.Dump(as->GetColor()[l], "cl" + sl);
        if (dl.count("nx" + sl)) m.Dump(as->GetNormal()[l], 0, "nx" + sl);
        if (dl.count("ny" + sl)) m.Dump(as->GetNormal()[l], 1, "ny" + sl);
        if (dl.count("nz" + sl)) m.Dump(as->GetNormal()[l], 2, "nz" + sl);
        if (dl.count("k" + sl)) m.Dump(&fck_[l], "k" + sl);
      }

      // combined colors
      if (dl.count("cls")) m.Dump(&as->GetColorSum(), "cls");

      // image
      auto conv = [&](size_t d, size_t l, Multi<FieldCell<Scal>>& fc) {
        fc.resize(layers);
        fc[l].Reinit(m);
        for (auto c : m.Cells()) {
          fc[l][c] = as->GetImage(l, c)[d];
        }
        return &fc[l];
      };
      for (auto d : {0, 1, 2}) {
        for (auto l : layers) {
          std::stringstream st;
          st << "im"
             << "xyz"[d] << l;
          std::string s = st.str();
          if (dl.count(s)) {
            m.Dump(conv(d, l, ctx->im[d]), s);
          }
        }
      }
    }
  }
  if (sem()) {
  } // XXX: empty stage, otherwise ctx is destroyed before dump
  if (var.Int["enable_advection"]) {
    if (var.Int["dumppoly"] && sem.Nested()) {
      as_->DumpInterface(GetDumpName("s", ".vtk", dumper_.GetN()));
    }
    if (var.Int["dumppolymarch"] && sem.Nested()) {
      as_->DumpInterfaceMarch(GetDumpName("sm", ".vtk", dumper_.GetN()));
    }
  }
}

template <class M>
void Hydro<M>::Dump() {
  auto sem = m.GetSem("dump");
  struct {
    Multi<FieldCell<MIdx>> fcim;
  } * ctx(sem);
  if (sem.Nested("fields")) {
    if (dumper_.Try(st_.t, st_.dt)) {
      DumpFields();
    }
  }
  if (dmptraj_.Try(st_.t, st_.dt)) {
    if (sem("copyimage")) {
      ctx->fcim.resize(layers);
      ctx->fcim.InitAll(FieldCell<MIdx>(m));
      if (auto as = dynamic_cast<ASVM*>(as_.get())) {
        for (auto c : m.AllCells()) {
          for (auto l : layers) {
            ctx->fcim[l][c] = as->GetImage(l, c);
          }
        }
      } else if (auto as = dynamic_cast<ASV*>(as_.get())) {
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
        } else if (auto as = dynamic_cast<ASVEB*>(as_.get())) {
          // TODO reuse ASV code
          fcu[0] = &as->GetField();
          fccl[0] = &as->GetColor();
        } else if (auto as = dynamic_cast<ASV*>(as_.get())) {
          fcu[0] = &as->GetField();
          fccl[0] = &as->GetColor();
        }
        DumpTraj<M>(
            m, true, var, dmptraj_.GetN(), st_.t, layers, fcu, fccl, ctx->fcim,
            fs_->GetPressure(), fs_->GetVelocity(), fcvm_, st_.dt);
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
      if (st_.step % var.Int("stat_step_every", 1) == 0) {
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

  if (var.Int["dumpinit"]) {
    if (sem.Nested()) {
      Dump();
    }
  }

  sem.LoopBegin();

  if (sem("events")) {
    if (events_) {
      events_->Exec(st_.t);
    }
  }
  if (sem("loop-check")) {
    if (st_.t + st_.dt * 0.25 > var.Double["tmax"] ||
        int(st_.step) >= var.Int["max_step"]) {
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
  if (var.Int["enable_fluid"]) {
    if (sem.Nested("fs-iters")) {
      StepFluid();
    }
  }
  if (sem.Nested("fs-finish")) {
    fs_->FinishStep();
  }

  if (var.Int["enable_advection"]) {
    if (sem.Nested("as-steps")) {
      StepAdvection();
    }
    if (sem.Nested("as-post")) {
      as_->PostStep();
    }
    if (sem.Nested("curv")) {
      psm_ = UCurv<M>::CalcCurvPart(layers, as_.get(), psm_par_, fck_, m);
    }
    if (var.Int["enable_bubgen"]) {
      if (sem.Nested("bubgen")) {
        StepBubgen();
      }
    }
  }

  if (sem.Nested("dt")) {
    CalcDt(); // must be after CalcStat to keep dt for moving mesh velocity
  }

  if (sem.Nested()) {
    Dump();
  }

  if (sem("inc")) {
    ++st_.step;
    m.CollectSample("Hydro::Step");
  }

  sem.LoopEnd();
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
auto Hydro<M>::CalcPressureDrag(const FieldCell<Scal>& fcp, const Embed<M>& eb)
    -> Vect {
  auto fep = eb.Interpolate(fcp, MapCondFace(), 1, 0.);
  Vect sum(0);
  for (auto c : eb.CFaces()) {
    sum += eb.GetSurface(c) * fep[c];
  }
  return sum;
}

template <class M>
auto Hydro<M>::CalcViscousDrag(
    const FieldCell<Vect>& fcvel, const FieldCell<Scal>& fcmu,
    const Embed<M>& eb) -> Vect {
  auto feg = eb.Gradient(fcvel, MapCondFace(), 0, Vect(0));
  auto femu = eb.Interpolate(fcmu, MapCondFace(), 1, 0.);
  Vect sum(0);
  for (auto c : eb.CFaces()) {
    sum += feg[c] * (-eb.GetArea(c) * femu[c]);
  }
  return sum;
}

template <class M>
void Hydro<M>::ReportIter() {
  std::cout << std::scientific << std::setprecision(16)
            << ".....iter=" << fs_->GetIter() << ", diff=" << diff_
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
  sem.LoopBegin();
  if (sem.Nested("iter")) {
    fs_->MakeIteration();
  }
  if (sem("reduce")) {
    diff_ = fs_->GetError();
    m.Reduce(&diff_, "max");
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
    assert(fs_->GetError() <= diff_);
    auto it = fs_->GetIter();
    if ((diff_ < var.Double["tol"] && (int)it >= var.Int["min_iter"]) ||
        (int)it >= var.Int["max_iter"]) {
      sem.LoopBreak();
    }
  }
  // TODO: Suspender loop hangs if (probably) Nested is last
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
            const_cast<FieldFace<Scal>&>(fs_->GetVolumeFlux());
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
}

template <class M>
void Hydro<M>::StepBubgen() {
  auto sem = m.GetSem("bubgen");
  struct {
    FieldCell<Scal> fcvf; // volume fraction
  } * ctx(sem);
  auto& fcvf = ctx->fcvf;
  const Scal t0 = var.Double["bubgen_t0"];
  const Scal tper = var.Double["bubgen_per"];
  bool bg = (st_.t > t0 && st_.t - bgt_ >= tper);
  if (bg) {
    if (sem("as-bubgen-init")) {
      if (!bgf_) {
        Vars vr;
        vr.String.Set("init_vf", "list");
        vr.String.Set("list_path", var.String["bubgen_path"]);
        vr.Int.Set("dim", var.Int["dim"]);
        vr.Int.Set("list_ls", var.Int["list_ls"]);
        bgf_ = CreateInitU<M>(vr, m.IsRoot());
      }
      fcvf.Reinit(m, 0);
      bgf_(fcvf, m);
      m.Comm(&fcvf);
    }
    if (sem("as-bubgen-apply")) {
      bgt_ = st_.t;
      if (auto as = dynamic_cast<ASVM*>(as_.get())) {
        auto& u = const_cast<FieldCell<Scal>&>(*as->GetFieldM()[0]);
        auto& cl = const_cast<FieldCell<Scal>&>(*as->GetColor()[0]);
        for (auto c : m.AllCells()) {
          if (fcvf[c] > 0) {
            u[c] = std::max(u[c], fcvf[c]);
            cl[c] = 1.;
          }
        }
      }
      if (auto as = dynamic_cast<ASV*>(as_.get())) {
        auto& u = const_cast<FieldCell<Scal>&>(as->GetField());
        for (auto c : m.AllCells()) {
          if (fcvf[c] > 0) {
            u[c] = std::max(u[c], fcvf[c]);
          }
        }
      }
    }
  }
}
