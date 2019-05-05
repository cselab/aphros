#pragma once

#include <cassert>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>
#include <array>
#include <list>
#include <chrono>
#include <thread>
#include <mpi.h>
#include <stdexcept>
#include <memory>
#include <limits>
#include <set>

#include "geom/mesh.h"
#include "solver/solver.h"
#include "solver/advection.h"
#include "solver/vof.h"
#include "parse/vof.h"
#include "solver/tvd.h"
#include "parse/tvd.h"
#include "solver/tracker.h"
#include "solver/sphavg.h"
#include "solver/simple.h"
#include "parse/simple.h"
#include "kernel.h"
#include "kernelmesh.h"
#include "parse/vars.h"
#include "parse/parser.h"
#include "parse/util.h"
#include "dump/output.h"
#include "dump/dumper.h"
#include "func/init_u.h"
#include "func/init_sig.h"
#include "func/init_cl.h"
#include "debug/isnan.h"
#include "solver/reconst.h"
#include "young/young.h"
#include "solver/pois.h"
#include "util/fluid.h"
#include "util/events.h"
#include "util/convdiff.h"

class GPar {};

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
  static constexpr size_t dim = M::dim;

  Hydro(Vars&, const MyBlockInfo&, Par&) ;
  void Run() override;
  M& GetMesh() { return m; }

 protected:
  using P::var;
  using P::bi_;
  using P::m;

 private:
  void Init();
  void InitFluid();
  void InitAdvection();
  void InitStat();
  void InitVort();
  // zero-gradient bc vect
  MapFace<std::shared_ptr<solver::CondFace>> GetBcVz() const;
  // zero-gradient bc scal
  MapFace<std::shared_ptr<solver::CondFace>> GetBcSz() const;
  void DumpFields();
  void Dump(Sem& sem);
  void DumpTraj(bool dm);
  // Calc rho, mu and force based on volume fraction
  void CalcMixture(const FieldCell<Scal>& vf);
  // fcvf: volume fraction on cells [a]
  // ffvfsm: smoothed volume fraction on faces [a]
  void CalcSurfaceTension(const FieldCell<Scal>& fcvf,
                          const FieldFace<Scal>& ffvfsm);
  // Clips v to given range, uses const_cast
  void Clip(const FieldCell<Scal>& v, Scal min, Scal max);
  void CalcStat();
  void CalcDt();
  void ReportStep();
  // Issue sem.LoopBreak if abort conditions met
  void CheckAbort(Sem& sem);
  void StepFluid();
  void StepAdvection();
  void StepBubgen();

  using FS = solver::Simple<M>;
  using AST = solver::Tvd<M>; // advection TVD
  using ASV = solver::Vof<M>; // advection VOF
  using TR = solver::Tracker<M>; // color tracker
  using SA = solver::Sphavg<M>; // spherical averages

  void UpdateAdvectionPar() {
    if (auto as = dynamic_cast<AST*>(as_.get())) {
      Parse<M>(as->GetPar(), var);
    } else if (auto as = dynamic_cast<ASV*>(as_.get())) {
      Parse<M>(as->GetPar(), var);
    }
  }
  // Surface tension time step
  Scal GetStDt() {
    Scal sig = var.Double["sigma"];
    Scal* cflst = var.Double("cflst");
    if (cflst && sig != 0.) {
      Scal pi = M_PI;
      Scal h3 = m.GetVolume(IdxCell(0));
      Scal r1 = var.Double["rho1"];
      Scal r2 = var.Double["rho2"];
      return (*cflst) * std::sqrt(h3 * (r1 + r2) / (4. * pi * std::abs(sig)));
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
  void CalcStrain() {
    auto& fcv = fs_->GetVelocity();
    auto& fcs = fc_strain_;

    auto ffv = solver::Interpolate(fcv, fs_->GetVelocityCond(), m);

    std::array<FieldCell<Vect>, dim> g; // g[i][c][j] is derivative du_i/dx_j
    for (size_t i = 0; i < dim; ++i) {
      g[i] = solver::Gradient(GetComponent(ffv, i), m);
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
    FieldCell<Scal> fc(m, 0); // result
    auto& ffv = fs_->GetVolumeFlux();

    for (auto c : m.Cells()) {
      for (auto q : m.Nci(c)) {
        IdxFace f = m.GetNeighbourFace(c, q);
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

  FieldCell<Scal> fc_mu_; // viscosity
  FieldCell<Scal> fc_rho_; // density
  FieldCell<Scal> fc_src_; // source of mixture volume
  FieldCell<Scal> fc_src2_; // source of second phase volume
  FieldCell<Scal> fc_srcm_; // mass source
  FieldCell<Vect> fc_force_;  // force 
  FieldFace<Scal> ffbp_;  // balanced force projections
  FieldFace<Scal> ffk_;  // curvature on faces
  MapFace<std::shared_ptr<solver::CondFace>> mf_cond_;
  MapFace<std::shared_ptr<solver::CondFaceFluid>> mf_velcond_;
  MapCell<std::shared_ptr<solver::CondCell>> mc_cond_;
  MapCell<std::shared_ptr<solver::CondCellFluid>> mc_velcond_;
  std::unique_ptr<solver::AdvectionSolver<M>> as_; // advection solver
  std::unique_ptr<FS> fs_; // fluid solver
  std::unique_ptr<TR> tr_; // color tracker
  std::unique_ptr<SA> sa_; // spherical averages
  FieldCell<Scal> fc_vf_; // volume fraction used by constructor 
  FieldCell<Scal> fccl_; // color used by constructor  
  FieldCell<Vect> fc_vel_; // velocity used by constructor
  FieldCell<Scal> fc_smvf_; // smoothed volume fraction used by CalcMixture()
  FieldFace<Scal> ff_smvf_; // interpolated fc_smvf_
  Scal diff_;  // convergence indicator
  std::pair<Scal, int> pdist_; // distance to pfixed cell

  FieldCell<Scal> fc_sig_; // surface tension sigma
  FieldFace<Scal> ff_sig_; // surface tension sigma

  std::vector<Scal> clr_cl_; // color reduce: cl
  std::vector<std::vector<Scal>> clr_v_; // color reduce: vector
  std::vector<std::string> clr_nm_; // color reduce: variable name
  std::map<Scal, std::vector<Scal>> clmp_; // tmp map color to vector (only root)
  FieldCell<Vect> fcvm_; // velocity field time_prev // TODO: revise

  FieldCell<Vect> fcyv_; // Young velocity
  FieldCell<Vect> fcom_; // vorticity
  FieldCell<Scal> fcomm_; // vorticity magnitude
  std::shared_ptr<solver::PoisSolver<M>> ps_; // Poisson solver for InitVort
  FieldCell<Scal> fcbc_; // boundary condition type by GetBcField()
  FieldCell<Scal> fc_strain_; // double inner product of strain rate tensor
  FieldCell<Scal> fctmp_;    // temporary scalar field
  FieldCell<Vect> fctmpv_;   // temporary vector field
  FieldCell<Vect> fcvcurl_;  // curl-component of velocity
  FieldCell<Scal> fcdiv_;  // divergence of velocity
  bool vcurl_;  // compute curl

  using Sph = typename SA::Sph;
  std::vector<Sph> sa_ss_;

  Scal nabort_; // number of abort requests, used by Reduce in checknan Run()

  std::function<void(FieldCell<typename M::Scal>&,const M&)> bgf_; // bubgen
  Scal bgt_ = -1.; // bubgen last time 

  struct Stat {
    Scal m1, m2;            // volume
    Scal m20;               // initial volume
    Scal m2d;               // relative volume difference
    Vect c1, c2;            // center of mass 
    Vect vc1, vc2;          // center of mass velocity
    Vect v1, v2;            // average velocity
    Scal dtt;               // temporary to reduce
    Scal dt;                // dt fluid 
    Scal dta;               // dt advection
    size_t iter;            // iter of fluid solver
    Scal dumpt = -1e10;     // last dump time (rounded to nearest dtdump)
    Scal t;
    size_t step;
    size_t dumpn;
    Vect meshpos;           // mesh position
    Vect meshvel;           // mesh velocity
    Scal ekin, ekin1, ekin2;// kinetic energy
    Scal workst;            // work by surface tension
    Vect vlm, vl2;          // max-norm and l2-norm of velocity minus "vel"
    Scal pmin, pmax, pd;    // pressure min,max
    Scal pavg1, pavg2;      // pressure average
    Scal boxomm;            // integral of vorticity magnitude over box
    Scal boxomm2;           // integral of vorticity magnitude over box
    Vect vomm;              // velocity weighted by vorticity
    Scal vommw;             // integral of vorticity
    Scal enstr;             // enstrophy
    Scal area;              // interface area
    Scal dissip, dissip1, dissip2;   // energy dissipation rate
    Scal edis, edis1, edis2;   // dissipated energy
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
      for (auto it : mst) {
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
  Events events_;  // events from var
};

template <class M>
auto Hydro<M>::GetBcVz() const -> MapFace<std::shared_ptr<solver::CondFace>> {
  // zero-derivative bc for Vect
  MapFace<std::shared_ptr<solver::CondFace>> r;
  for (auto it : mf_velcond_) {
    IdxFace f = it.GetIdx();
    r[f] = std::make_shared<solver::
      CondFaceGradFixed<Vect>>(Vect(0), it.GetValue()->GetNci());
  }
  return r;
}

template <class M>
auto Hydro<M>::GetBcSz() const -> MapFace<std::shared_ptr<solver::CondFace>> {
  // zero-derivative bc for Vect
  MapFace<std::shared_ptr<solver::CondFace>> r;
  for (auto it : mf_velcond_) {
    IdxFace f = it.GetIdx();
    r[f] = std::make_shared<solver::
      CondFaceGradFixed<Scal>>(0., it.GetValue()->GetNci());
  }
  return r;
}

// Compute velocity fc_vel_ from vorticity stored in fc_vel_
template <class M>
void Hydro<M>::InitVort() {
  auto sem = m.GetSem("initvort");
  auto& fct = fctmp_;  // temporary fields
  auto& fctv = fctmpv_;
  if (sem("initpois")) {
    m.Comm(&fc_vel_);
    fctv.Reinit(m);
  }
  for (size_t d = 0; d < M::dim; ++d) {
    std::string dn = std::to_string(d);
    if (sem("init-" + dn)) {
      ps_ = std::make_shared<solver::PoisSolver<M>>(
          GetScalarCond(fs_->GetVelocityCond(), d, m), m);
      fct = GetComponent(fc_vel_, d);
      for (auto c : m.Cells()) {
        fct[c] *= -1;
      }
    }
    if (sem.Nested("solve-" + dn)) {
      ps_->Solve(fct);
    }
    if (sem("post-" + dn)) {
      SetComponent(fctv, d, ps_->GetField());
      if (m.IsRoot() && var.Int["vort_report"]) {
        std::cout 
            << "om" << ("xyz"[d]) << ":" 
            << " res=" << m.GetResidual()
            << " iter=" << m.GetIter()
            << std::endl;
      }
    }
  }
  if (sem("vel")) {
    fc_vel_ = GetVort(fctv, fs_->GetVelocityCond(), m);
    m.Comm(&fc_vel_);
    fctv.Free();
    fct.Free();
  }
}

template <class M>
void Hydro<M>::InitFluid() {
  auto p = std::make_shared<typename FS::Par>();
  Parse<M>(p.get(), var);

  fcvm_ = fc_vel_;


  // XXX ahoc: young velocity
  if (var.Int["youngbc"]) {
    fcyv_.Reinit(m);
    InitYoung();
    for (auto c : m.Cells()) {
      Vect x = m.GetCenter(c);
      fcyv_[c] = GetYoungVel(x);
    }
  }

  fs_.reset(new FS(
        m, fc_vel_, mf_velcond_, mc_velcond_, 
        &fc_rho_, &fc_mu_, &fc_force_, &ffbp_,
        &fc_src_, &fc_srcm_, 0., st_.dt, p));

  fcbc_ = GetBcField(mf_velcond_, m);
}

template <class M>
void Hydro<M>::InitAdvection() {
  std::string as = var.String["advection_solver"];
  if (as == "tvd") {
    auto p = std::make_shared<typename AST::Par>();
    Parse<M>(p.get(), var);
    as_.reset(new AST(
          m, fc_vf_, mf_cond_, 
          &fs_->GetVolumeFlux(solver::Layers::time_curr),
          &fc_src2_, 0., st_.dta, p));
  } else if (as == "vof") {
    auto p = std::make_shared<typename ASV::Par>();
    Parse<M>(p.get(), var);
    p->dmp = std::unique_ptr<Dumper>(new Dumper(var, "dump_part_"));
    as_.reset(new ASV(
          m, fc_vf_, mf_cond_, 
          &fs_->GetVolumeFlux(solver::Layers::time_curr),
          &fc_src2_, 0., st_.dta, p));
  } else {
    throw std::runtime_error("Unknown advection_solver=" + as);
  }
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
    auto op = [](std::string n,  Scal* p) {
      return std::make_shared<output::OutScalFunc<Scal>>(
          n, [p](){ return *p; });
    };

    auto& s = st_;
    output::VOut con = {
        op("t", &s.t),
        std::make_shared<output::OutScalFunc<int>>(
            "iter", [this](){ return st_.iter; }),
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
        op("c1x", &s.c1[0]), op("c1y", &s.c1[1]), op("c1z", &s.c1[2]),
        op("c2x", &s.c2[0]), op("c2y", &s.c2[1]), op("c2z", &s.c2[2]),
        op("vc1x", &s.vc1[0]), op("vc1y", &s.vc1[1]), op("vc1z", &s.vc1[2]),
        op("vc2x", &s.vc2[0]), op("vc2y", &s.vc2[1]), op("vc2z", &s.vc2[2]),
        op("v1x", &s.v1[0]), op("v1y", &s.v1[1]), op("v1z", &s.v1[2]),
        op("v2x", &s.v2[0]), op("v2y", &s.v2[1]), op("v2z", &s.v2[2]),
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
    ost_ = std::make_shared<output::SerScalPlain<Scal>>(con, "stat.dat");
  }
}

template <class M>
void Hydro<M>::Init() {
  auto sem = m.GetSem("init");

  if (sem("fields")) {
    fc_src_.Reinit(m, 0.);
    fc_src2_.Reinit(m, 0.);
    fc_srcm_.Reinit(m, 0.);

    // initial volume fraction
    fc_vf_.Reinit(m, 0);
    auto ivf = CreateInitU<M>(var, m.IsRoot());
    ivf(fc_vf_, m);
    m.Comm(&fc_vf_);

    // initial surface tension sigma
    fc_sig_.Reinit(m, 0);
    auto isig = CreateInitSig<M>(var);
    isig(fc_sig_, m);
    m.Comm(&fc_sig_);

    // initial velocity
    fc_vel_.Reinit(m, Vect(0));
    InitVel(fc_vel_, var, m);
    m.Comm(&fc_vel_);

    // global mesh size
    MIdx gs = m.GetGlobalSize();

    if (m.IsRoot()) {
      std::cout << "global mesh=" << gs << std::endl;
      std::cout << "surface tension dt=" << GetStDt() << std::endl;
    }

    // boundary conditions
    GetFluidFaceCond(var, m, mf_velcond_, mf_cond_);
  }

  if (sem.Nested("initvort")) {
    if (var.Int["initvort"]) {
      InitVort();
    }
  }

  if (sem.Nested("smooth")) {
    solver::Smoothen(fc_vf_, mf_cond_, m, var.Int["vf_init_sm"]);
  }

  if (sem.Nested("mixture")) {
    CalcMixture(fc_vf_);
  }

  if (sem("color-ini")) {
    if (var.Int["enable_color"]) {
      // initial color
      auto icl = CreateInitCl<M>(var, m.IsRoot());
      icl(fccl_, fc_vf_, m);
      m.Comm(&fccl_);
    } else {
      fccl_.Reinit(m, 0.);
    }
  }

  if (sem.Nested("cellcond")) {
    GetFluidCellCond(var, m, mc_velcond_, pdist_);
  }

  if (sem("solv")) {
    // time step
    const Scal dt = var.Double["dt0"];
    st_.dt = dt;
    st_.dta = dt;
    var.Double.Set("dt", st_.dt);
    var.Double.Set("dta", st_.dta);


    InitFluid();

    InitAdvection();

    // Init color tracker
    if (var.Int["enable_color"] && as_) {
      tr_.reset(new TR(m, fccl_, var.Double["color_th"], var.Int["dim"]));
    }

    // Init sphavg
    if (var.Int["enable_shell"] && tr_) {
      sa_.reset(new SA(m, var.Int["dim"]));
    }

    st_.iter = 0;
    var.Int.Set("iter", st_.iter);

    st_.t = fs_->GetTime();
    var.Double.Set("t", st_.t);

    InitStat();

    // enable curl component of velocity 
    {
      auto dl = GetWords(var.String["dumplist"]);
      if (dl.count("vcurlx") || dl.count("vcurly") || dl.count("vcurlz")) { 
        vcurl_ = true;
      } else {
        vcurl_ = false;
      }
    }

    events_.Parse();
  }
}


template <class M>
Hydro<M>::Hydro(Vars& var, const MyBlockInfo& bi, Par& par) 
    : KernelMeshPar<M,Par>(var, bi, par)
    , st_{}
    , dumper_(var, "dump_field_")
    , dmptraj_(var, "dump_traj_")
    , dmptrep_(var, "dump_trep_")
    , events_(var, m.IsRoot(), m.IsLead())
{}

template <class M>
void Hydro<M>::CalcStat() {
  auto sem = m.GetSem("stat");

  auto& s = st_;
  auto& fa = as_->GetField();
  auto& fv = fs_->GetVelocity();
  auto& fp = fs_->GetPressure();

  if (sem("stat-add")) {
    s.Add(fa, "vf", m);
    s.Add(as_->GetCurv(), "k", m);
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

    // XXX: adhoc: also controls s.vomm
    if (var.Int["statbox"]) {
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
      CalcVort();
      s.enstr = 0.;
      for (auto c : m.Cells()) {
        s.enstr += 0.5 * sqr(fcomm_[c]) * fc_rho_[c] * m.GetVolume(c);
      }
      m.Reduce(&s.enstr, "sum");
    }
    // surface area
    if (var.Int["stat_area"]) {
      s.area = 0;
      if (auto as = dynamic_cast<solver::Vof<M>*>(as_.get())) {
        using R = Reconst<Scal>;
        auto &fcn = as->GetNormal();
        auto &fca = as->GetAlpha();
        auto& fcvf = as->GetField();
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
      CalcStrain();
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

    if (std::string* s = var.String("meshvel_auto")) {
      Vect v(0);
      if (*s == "v") {
        v = st_.v2;
      } else if (*s == "vc") {
        v = st_.vc2;
      } else if (*s == "vomm") {
        v = st_.vomm;
      } else {
        throw std::runtime_error("Unknown meshvel_auto=" + *s);
      }
      Vect mask(var.Vect["meshvel_mask"]); // components 0 or 1
      v *= mask;
      double w = var.Double["meshvel_weight"];
      Vect vp = fs_->GetPar()->meshvel;
      fs_->GetPar()->meshvel = v * w + vp * (1. - w);

      st_.meshvel = fs_->GetPar()->meshvel;
      st_.meshpos += st_.meshvel * st_.dt;
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
      // XXX: adhoc, overwrite wall conditions
      auto& fa = as_->GetField();
      auto& fv = fs_->GetVelocity();
      for (auto it : mf_velcond_) {
        IdxFace f = it.GetIdx();
        solver::CondFaceFluid* cb = it.GetValue().get();
        if (auto cd = dynamic_cast<solver::fluid_condition::
            NoSlipWallFixed<M>*>(cb)) {
          size_t nci = cd->GetNci();
          Vect n = m.GetNormal(f);
          IdxCell c = m.GetNeighbourCell(f, nci);
          auto v = fv[c];
          cd->SetVelocity((v - n * n.dot(v)) * std::min(1., fa[c] * kslip));
        } 
      }
    }
  }

  if (sem("young")) {
    if (var.Int["youngbc"]) {
      InitYoung();
      for (auto it : mf_velcond_) {
        IdxFace f = it.GetIdx();
        Vect x = m.GetCenter(f);
        solver::CondFaceFluid* cb = it.GetValue().get();
        if (auto cd = dynamic_cast<solver::fluid_condition::
            NoSlipWallFixed<M>*>(cb)) {
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
    var.Double["t"] = st_.t;

    st_.dtt = fs_->GetAutoTimeStep();
    m.Reduce(&st_.dtt, "min");
  }
  if (sem("reduce")) {
    // set from cfl if defined
    if (auto* cfl = var.Double("cfl")) {
      st_.dt = st_.dtt * (*cfl);
      st_.dt = std::min<Scal>(st_.dt, var.Double["dtmax"]);
    }

    // constraint from surface tension
    st_.dt = std::min<Scal>(st_.dt, GetStDt());

    fs_->SetTimeStep(st_.dt);
    var.Double["dt"] = st_.dt;

    // set from cfla if defined
    if (auto* cfla = var.Double("cfla")) {
      st_.dta = st_.dtt * (*cfla); 
      st_.dta = std::min<Scal>(st_.dta, var.Double["dtmax"]);
    }

    // round up dta to such that dt / dta is integer
    Scal dt = fs_->GetTime() + fs_->GetTimeStep() - as_->GetTime();
    st_.dta = dt / std::max(1, int(dt / st_.dta + 0.5));

    as_->SetTimeStep(st_.dta);
    var.Double["dta"] = st_.dta;
  }
}

template <class M>
void Hydro<M>::Clip(const FieldCell<Scal>& f, Scal a, Scal b) {
  auto& g = const_cast<FieldCell<Scal>&>(f);
  for (auto i : m.Cells()) {
    g[i] = std::max(a, std::min(b, g[i]));
  }
}

template <class M>
void Hydro<M>::CalcSurfaceTension(const FieldCell<Scal>& fcvf,
                                  const FieldFace<Scal>& ffvfsm) {
  // volume fration gradient on cells
  FieldCell<Vect> gc = solver::Gradient(ffvfsm, m); // [s]
  // volume fration gradient on faces
  FieldFace<Vect> gf = solver::Interpolate(gc, GetBcVz(), m); // [i]

  auto st = var.String["surftens"];
  if (st == "div") { // divergence of tensor (Hu,Adam 2001)
    auto stdiag = var.Double["stdiag"];
    for (auto c : m.Cells()) {
      Vect r(0);
      for (auto q : m.Nci(c)) {
        IdxFace f = m.GetNeighbourFace(c, q);
        const auto& g = gf[f];
        // TODO: revise 1e-6
        auto n = g / (g.norm() + 1e-6);  // inner normal
        auto s = m.GetOutwardSurface(c, q);
        r += s * (g.norm() * stdiag) - g * s.dot(n);
      }
      r /= m.GetVolume(c);
      // here: r = stdiag*div(|g|I) - div(g g/|g|)
      fc_force_[c] += r * fc_sig_[c];
    }
  } else if (st == "kn") {  // curvature * normal
    auto& fck = as_->GetCurv(); // [a]
    FieldFace<Scal> ff_st(m, 0.);  // surface tension projections

    ffk_.Reinit(m, 0);
    ff_sig_.Reinit(m, 0);
    ff_sig_ = solver::Interpolate(fc_sig_, GetBcSz(), m);
    // interpolate curvature
    for (auto f : m.Faces()) {
      IdxCell cm = m.GetNeighbourCell(f, 0);
      IdxCell cp = m.GetNeighbourCell(f, 1);
      if (std::abs(fcvf[cm] - 0.5) < std::abs(fcvf[cp] - 0.5)) {
        ffk_[f] = fck[cm];
      } else {
        ffk_[f] = fck[cp];
      }
    }
    // neighbour cell for boundaries
    for (auto it : mf_velcond_) {
      IdxFace f = it.GetIdx();
      IdxCell c = m.GetNeighbourCell(f, it.GetValue()->GetNci());
      ffk_[f] = fck[c];
    }
    // compute force projections
    size_t nan = 0;
    IdxFace fnan(0);
    bool report_knan = var.Int["report_knan"];
    for (auto f : m.Faces()) {
      IdxCell cm = m.GetNeighbourCell(f, 0);
      IdxCell cp = m.GetNeighbourCell(f, 1);
      // XXX: adhoc uniform (should be half volume cp,cm)
      Scal hr = m.GetArea(f) / m.GetVolume(cp);
      Scal ga = (fcvf[cp] - fcvf[cm]) * hr; 
      if (ga != 0.) {
        if (IsNan(ffk_[f])) {
          ++nan;
          fnan = f;
          ffk_[f] = 0.;
        }
        ff_st[f] += ga * ffk_[f] * ff_sig_[f];
      }
    }
    if (report_knan && nan) {
      auto f = fnan;
      IdxCell cm = m.GetNeighbourCell(f, 0);
      IdxCell cp = m.GetNeighbourCell(f, 1);
      std::stringstream s;
      s.precision(16);
      s << "nan curvature in " << nan 
          << " faces, one at x=" << m.GetCenter(f)
          << " vf[cm]=" << fcvf[cm]
          << " vf[cp]=" << fcvf[cp]
          << " fck[cm]=" << fck[cm]
          << " fck[cp]=" << fck[cp];
      std::cout << s.str() << std::endl;
    }
    // zero on boundaries
    // TODO: test with bubble jump
    for (auto it : mf_velcond_) {
      IdxFace f = it.GetIdx();
      ff_st[f] = 0.;
    }

    // Surface tension decay between x0 and x1 
    // XXX: adhoc TODO: revise
    const Scal x0 = var.Double["zerostx0"];
    const Scal x1 = var.Double["zerostx1"];
    // apply
    for (auto f : m.Faces()) {
      Scal x = m.GetCenter(f)[0];
      if (x > x0) {
        ff_st[f] *= std::max(0., (x1 - x) / (x1 - x0));
      }
    }

    // Append to force
    for (auto f : m.Faces()) {
      ffbp_[f] += ff_st[f];
    }

    // Append Marangoni stress
    if (var.Int["marangoni"]) {
      if (auto as = dynamic_cast<solver::Vof<M>*>(as_.get())) {
        using R = Reconst<Scal>;
        auto fc_gsig = solver::Gradient(ff_sig_, m);
        auto &fcn = as->GetNormal();
        auto &fca = as->GetAlpha();
        Vect h = m.GetCellSize();
        for (auto c : m.Cells()) {
          if (fcvf[c] > 0. && fcvf[c] < 1. && !IsNan(fca[c])) {
            Vect g = fc_gsig[c]; // sigma gradient
            Vect n = fcn[c] / fcn[c].norm(); // unit normal to interface
            Vect gt = g - n * g.dot(n); 
            auto xx = R::GetCutPoly2(fcn[c], fca[c], h);
            Scal ar = std::abs(R::GetArea(xx, fcn[c]));
            Scal vol = h.prod();
            fc_force_[c] += gt * (ar / vol);
          }
        }
      }
    }
  } else {
    throw std::runtime_error("Unknown surftens=" + st);
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
      fc_src_.Reinit(m, 0);
      for (auto c : m.Cells()) {
        Vect x = m.GetCenter(c);
        Scal sy = std::cos(2. * pi * x[1] / force_wly);
        fc_force_[c][0] += s * sy * force_mag * fc_vf0[c];
      }
    }
  }

  if (sem.Nested("smooth")) {
    solver::Smoothen(fc_smvf_, mf_cond_, m, var.Int["vfsmooth"]);
  }

  if (sem("calc")) {
    FieldCell<Scal>& a = fc_smvf_;
    FieldFace<Scal>& af = ff_smvf_;
    af = solver::Interpolate(a, mf_cond_, m);

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
      CalcSurfaceTension(fc_vf0, af);
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
  if (vcurl_) {
    if (sem("vcurl-pre")) {
      CalcVort();
      fc_vel_ = fcom_;
    }
    if (sem.Nested("vcurl-solve")) {
      InitVort();
    }
    if (sem("vcurl-copy")) {
      fcvcurl_ = fc_vel_;
    }
  }
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
    if (dl.count("k")) m.Dump(&as_->GetCurv(), "k");
    if (dl.count("rho")) m.Dump(&fc_rho_, "rho");
    if (dl.count("mu")) m.Dump(&fc_mu_, "mu");
    if (dl.count("sig")) m.Dump(&fc_sig_, "sig");
    if (dl.count("bc")) m.Dump(&fcbc_, "bc");
    if (tr_) {
      if (dl.count("cl")) m.Dump(&tr_->GetColor(), "cl");
      auto& im = tr_->GetImage();
      if (dl.count("imx")) m.Dump(&im, 0, "imx");
      if (dl.count("imy")) m.Dump(&im, 1, "imy");
      if (dl.count("imz")) m.Dump(&im, 2, "imz");
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
    if (dl.count("vcurlx")) m.Dump(&fcvcurl_, 0, "vcurlx");
    if (dl.count("vcurly")) m.Dump(&fcvcurl_, 1, "vcurly");
    if (dl.count("vcurlz")) m.Dump(&fcvcurl_, 2, "vcurlz");
    if (dl.count("dis") || dl.count("strain")) {
      CalcStrain();
      if (dl.count("strain")) m.Dump(&fcomm_, "strain");
      if (dl.count("dis")) {
        fctmp_ = fc_strain_;
        for (auto c : m.Cells()) {
          fctmp_[c] *= 2. * fc_mu_[c];
        }
        m.Dump(&fctmp_, "dis");
      }
    }
    if (dl.count("div")) {
     fcdiv_ = GetDiv();
     m.Dump(&fcdiv_, "div");
    }
  }
}

template <class M>
void Hydro<M>::Dump(Sem& sem) {
  if (sem.Nested("fields")) {
    if (dumper_.Try(st_.t, st_.dt)) {
      DumpFields();
    }
  }
  if (tr_) {
    if (sem.Nested("trajdump")) {
      if (dmptraj_.Try(st_.t, st_.dt)) {
        DumpTraj(1);
      }
    }
  }
  if (sem("dmptrep")) {
    if (m.IsRoot() && dmptrep_.Try(st_.t, st_.dt)) {
      std::string s = GetDumpName("trep", ".log", dmptrep_.GetN());
      m.TimerReport(s);
      std::cout << std::fixed << std::setprecision(8)
          << "timer report" 
          << " t=" << st_.t
          << " to " << s << std::endl;
    }
  }
  if (sem("dumpstat")) {
    if (m.IsRoot()) {
      ost_->Write(0., "");
    }
  }
}

template <class M>
void Hydro<M>::DumpTraj(bool dm) {
  auto sem = m.GetSem("dumptraj");
  if (sem("color-calc")) {
    std::map<Scal, std::vector<Scal>> mp; // map color to vector
    auto kNone = TR::kNone;
    auto& cl = tr_->GetColor();
    auto& im = tr_->GetImage();
    auto& vf = as_->GetField();
    auto& vel = fs_->GetVelocity();
    auto& p = fs_->GetPressure();
    Vect gh = m.GetGlobalLength(); // global domain length

    // add scalar name
    auto nma = [this](const std::string nm) {
      clr_nm_.push_back(nm);
    };
    // add vector name
    auto nmav = [this](const std::string nm) {
      clr_nm_.push_back(nm + "x");
      clr_nm_.push_back(nm + "y");
      clr_nm_.push_back(nm + "z");
    };

    // XXX: adhoc, the following order assumed in post: vs,r,x,y,z,...
    
    // list of vars // TODO: revise
    clr_nm_.clear();
    nma("vf");
    nma("r");
    nmav("");
    nmav("v");
    nma("p");

    // traverse cells, append to mp
    for (auto c : m.Cells()) {
      if (cl[c] != kNone) {
        auto& v = mp[cl[c]]; // vector for data
        auto x = m.GetCenter(c); // cell center
        x += im[c] * gh;  // translation by image vector

        auto w = vf[c] * m.GetVolume(c); // volume

        size_t i = 0;
        // append scalar value
        auto add = [&v,&i](Scal a) {
          if (i >= v.size()) {
            v.resize(i + 1);
          }
          v[i] += a;
          ++i;
        };
        // append vector value 
        auto addv = [&](Vect a) {
          add(a[0]);
          add(a[1]);
          add(a[2]);
        };

        // list of vars, XXX: keep consistent with clr_nm_ 
        add(w); // vf,  XXX: adhoc, vf must be first, divided on dump
        add(0.); // r,  XXX: adhoc, r must be second, computed on dump
        addv(x * w); // x
        addv(vel[c] * w); // v
        add(p[c] * w); // p
      }
    }
    // copy to vector
    clr_cl_.clear();
    clr_v_.clear();
    for (auto it : mp) {
      clr_cl_.push_back(it.first); // color
      clr_v_.push_back(it.second); // vector
    }
    using TS = typename M::template OpCatT<Scal>; 
    using TVS = typename M::template OpCatVT<Scal>; 
    m.Reduce(std::make_shared<TS>(&clr_cl_));
    m.Reduce(std::make_shared<TVS>(&clr_v_));
  }
  if (sem("color-post")) {
    if (m.IsRoot()) {
      // root has concatenation of all clr_cl_ and clr_v_
      if (clr_cl_.size() != clr_v_.size()) {
        throw std::runtime_error(
            "color-reduce: clr_cl_.size() != clr_v_.size()");
      }

      auto& mp = clmp_; // map color to vector
      mp.clear();
      // reduce to map
      for (size_t k = 0; k < clr_cl_.size(); ++k) {
        auto cl = clr_cl_[k];
        auto& v = clr_v_[k];
        auto& vm = mp[cl];
        vm.resize(v.size(), 0.);
        for (size_t i = 0; i < v.size(); ++i) {
          vm[i] += v[i];
        }
      }

      // divide by vf
      for (auto& it : mp) {
        auto& v = it.second;
        Scal vf = v[0]; // XXX: assume vf is first
        Scal pi = M_PI;
        // XXX: assume r is second
        v[1] = std::pow(3. / (4. * pi) * vf, 1. / 3.) ; 
        // divide remaining by vf
        for (size_t i = 2; i < v.size(); ++i) {
          v[i] /= vf;
        }
      }

      clr_cl_.clear();
      clr_v_.clear();
      for (auto& it : mp) {
        clr_cl_.push_back(it.first);
        clr_v_.push_back(it.second);
      }
    }

    using TS = typename M::template OpCatT<Scal>; 
    using TVS = typename M::template OpCatVT<Scal>; 
    m.Bcast(std::make_shared<TS>(&clr_cl_));
    m.Bcast(std::make_shared<TVS>(&clr_v_));
  }
  if (sa_ && sem("sphavg-sph")) {
    const Scal shrr = var.Double["shell_rr"]; // shell inner radius relative 
                                              // to equivalent radius
    const Scal shr = var.Double["shell_r"]; // shell inner radius absolute
    // shell total radius: rr * req + r
    const Scal shh = var.Double["shell_h"]; // shell thickness relative to h
    auto h = m.GetCellSize();

    auto& ss = sa_ss_;
    ss.clear();
    for (size_t i = 0; i < clr_v_.size(); ++i) {
      auto& s = clr_v_[i];
      // XXX: adhoc, assume vf,r,x,y,z in clr_v_
      Vect x(s[2], s[3], s[4]);
      Scal r = s[1] * shrr + shr;
      ss.emplace_back(x, r, h[0] * shh);
    }
  }
  if (sa_ && sem.Nested("sphavg-update")) {
    auto& vf = as_->GetField();
    auto& vel = fs_->GetVelocity();
    auto& velm = fcvm_;
    auto& p = fs_->GetPressure();
    sa_->Update(vf, vel, velm, st_.dt, p, sa_ss_);
  }
  if (sem("color-dump") && dm) {
    if (m.IsRoot()) {
      std::string s = GetDumpName("traj", ".csv", dmptraj_.GetN());
      std::cout << std::fixed << std::setprecision(8)
          << "dump" 
          << " t=" << st_.t
          << " to " << s << std::endl;
      std::ofstream o;
      o.open(s);
      o.precision(20);
      // header
      {
        o << "cl";
        auto& nm = clr_nm_;
        for (size_t i = 0; i < nm.size(); ++i) {
          o << "," << nm[i];
        }
        o << std::endl;
      }
      // content
      for (size_t i = 0; i < clr_cl_.size(); ++i) {
        auto cl = clr_cl_[i];
        auto& v = clr_v_[i];
        o << cl;
        for (size_t j = 0; j < v.size(); ++j) {
          o << "," << v[j];
        }
        o << "\n";
      }
    }
    if (sa_) {
      if (m.IsRoot()) {
        std::string s = GetDumpName("trajsh", ".csv", dmptraj_.GetN());
        std::cout << std::fixed << std::setprecision(8)
            << "dump" 
            << " t=" << st_.t
            << " to " << s << std::endl;
        std::ofstream o;
        o.open(s);
        o.precision(20);
        // header
        {
          auto nn = sa_->GetNames();
          o << "cl";
          for (auto n : nn) {
            o << "," << n;
          }
          o << std::endl;
        }
        // content
        auto cl = clr_cl_;
        auto av = sa_->GetAvg();
        if (cl.size() != av.size()) {
          throw std::runtime_error(
              "trajsh: cl.size()=" + std::to_string(cl.size()) +
              " != av.size()=" + std::to_string(av.size()));
        }
        for (size_t i = 0; i < cl.size(); ++i) {
          o << cl[i];
          for (auto& a : av[i].SerOut()) {
            o << "," << a;
          }
          o << "\n";
        }
      }
    }
  }
}

template <class M>
void Hydro<M>::Run() {
  auto sem = m.GetSem("run");

  if (sem.Nested("init")) {
    Init();
  }

  if (var.Int["dumpinit"]) {
    Dump(sem);
  }

  sem.LoopBegin();

  if (sem("events")) {
    events_.Exec(st_.t);
  }
  if (sem("loop-check")) {
    if (st_.t + st_.dt * 0.25 > var.Double["tmax"] || 
        int(st_.step) >= var.Int["max_step"]) {
      sem.LoopBreak();
    } else {
      if (m.IsRoot()) {
        ReportStep();
      }
    }
  }

  CheckAbort(sem);

  if (sem("updatepar")) {
    Parse<M>(fs_->GetPar(), var);
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
    if (var.Int["clip_vf"]) {
      if (sem("as-clip")) {
        Clip(as_->GetField(), 0., 1.);
      }
    }
    if (var.Int["enable_bubgen"]) {
      if (sem.Nested("bubgen")) {
        StepBubgen();
      }
    }
  }
  if (var.Int["enable_color"] && as_) {
    if (sem.Nested("color")) {
      tr_->Update(as_->GetField());
    }
  }

  if (sem.Nested("dt")) {
    CalcDt(); // must be after CalcStat to keep dt for moving mesh velocity 
  }

  Dump(sem);

  if (sem("inc")) {
    ++st_.step;
  }

  sem.LoopEnd();
}

template <class M>
void Hydro<M>::ReportStep() {
  std::cout << std::fixed << std::setprecision(8)
      << "STEP=" << st_.step 
      << " t=" << st_.t
      << " dt=" << st_.dt
      << " ta=" << as_->GetTime()
      << " dta=" << as_->GetTimeStep()
      << std::endl;
}

template <class M>
void Hydro<M>::CheckAbort(Sem& sem) {
  if (sem("abort-local")) {
    nabort_ = 0.;
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
    }
    catch (const std::runtime_error& e) {
      std::cout << e.what() << std::endl;
      nabort_ += 1.;
    }
    m.Reduce(&nabort_, "sum");
  }

  if (sem("abort-reduce")) {
    if (nabort_ != 0.) {
      if (m.IsRoot()) {
        std::cout << "nabort_ = " << nabort_ << std::endl;
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
    var.Int["iter"] = st_.iter;
    if (m.IsRoot()) {
      std::cout << std::scientific << std::setprecision(16)
          << ".....iter=" << fs_->GetIter()
          << ", diff=" << diff_ << std::endl;
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
      std::cout << std::fixed << std::setprecision(8)
          << ".....adv: t=" << as_->GetTime() 
          << " dt=" << as_->GetTimeStep()
          << std::endl;
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
  Scal t0 = var.Double["bubgen_t0"];
  Scal tper = var.Double["bubgen_per"];
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
      fc_vf_.Reinit(m, 0);
      bgf_(fc_vf_, m);
      m.Comm(&fc_vf_);
    }
    if (sem("as-bubgen-apply")) {
      auto& u = const_cast<FieldCell<Scal>&>(as_->GetField());
      for (auto c : m.AllCells()) {
        u[c] = std::max(u[c], fc_vf_[c]);
      }
      bgt_ = st_.t;
    }
  }
}
