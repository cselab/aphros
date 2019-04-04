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
  using P::IsRoot;
  using P::IsLead;

 private:
  void Init();
  static FieldCell<Vect> GetVort(
      const FieldCell<Vect>& fcv, 
      const MapFace<std::shared_ptr<solver::CondFace>>& mf, M& m);
  void InitVort();
  // zero-gradient bc vect
  MapFace<std::shared_ptr<solver::CondFace>> GetBcVz() const;
  // zero-gradient bc scal
  MapFace<std::shared_ptr<solver::CondFace>> GetBcSz() const;
  void Dump(Sem& sem);
  void DumpTraj(bool dm);
  void CalcMixture(const FieldCell<Scal>& vf);
  // Clips v to given range, uses const_cast
  void Clip(const FieldCell<Scal>& v, Scal min, Scal max);
  void CalcStat();
  void CalcDt();

  struct Event {
    double t;
    std::string cmd;
    std::string arg;
  };
  std::map<std::string, Event> ev_;
  // Parse events from var.String and put to ev_
  void ParseEvents();
  // Exec events due and remove from ev_
  void ExecEvents();

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
  Vect GetCellSize() const {
    Vect h; // cell size
    // XXX: specific for structured 3D mesh
    IdxCell c0(0);
    h = m.GetNode(m.GetNeighbourNode(c0, 7)) - 
        m.GetNode(m.GetNeighbourNode(c0, 0));
    assert(std::abs(h.prod() - m.GetVolume(c0)) < 1e-10);
    return h;
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
    fcom_ = GetVort(fcv, GetBcVz(), m);
    fcomm_.Reinit(m);
    for (auto c : m.Cells()) {
      fcomm_[c] = fcom_[c].norm();
    }
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
  // Returns field with the type (index)
  // of boundary conditions in an adjacent face:
  //   0: empty
  //   1: no-slip wall
  //   2: free-slip wall
  //   3: inlet
  //   4: outlet
  //   -1: unknown
  // mf: boundary conditions
  FieldCell<Scal> GetBcField(
      MapFace<std::shared_ptr<solver::CondFaceFluid>>& mf, const M& m) {
    FieldCell<Scal> fc(m, 0);
    for (auto it : mf) {
      IdxFace f = it.GetIdx();
      auto* b = it.GetValue().get();
      size_t nci = b->GetNci();
      IdxCell c = m.GetNeighbourCell(f, nci);
      if (dynamic_cast<solver::fluid_condition::NoSlipWall<M>*>(b)) {
        fc[c] = 1.;
      } else if (dynamic_cast<solver::fluid_condition::SlipWall<M>*>(b)) {
        fc[c] = 2.;
      } else if (dynamic_cast<solver::fluid_condition::Inlet<M>*>(b)) {
        fc[c] = 3.;
      } else if (dynamic_cast<solver::fluid_condition::Outlet<M>*>(b)) {
        fc[c] = 4.;
      } else {
        fc[c] = -1.;
      }
    }
    return fc;
  }

  FieldCell<Scal> fc_mu_; // viscosity
  FieldCell<Scal> fc_rho_; // density
  FieldCell<Scal> fc_src_; // source
  FieldCell<Vect> fc_force_;  // force 
  FieldFace<Scal> ffbp_;  // balanced force projections
  FieldFace<Scal> ffk_;  // curvature on faces
  MapFace<std::shared_ptr<solver::CondFace>> mf_cond_;
  MapFace<std::shared_ptr<solver::CondFaceFluid>> mf_velcond_;
  std::unique_ptr<solver::AdvectionSolver<M>> as_; // advection solver
  std::unique_ptr<FS> fs_; // fluid solver
  std::unique_ptr<TR> tr_; // color tracker
  std::unique_ptr<SA> sa_; // spherical averages
  FieldCell<Scal> fc_vf_; // volume fraction used by constructor 
  FieldCell<Scal> fccl_; // color used by constructor  
  FieldCell<Vect> fc_vel_; // velocity used by constructor
  FieldCell<Scal> fc_smvf_; // smoothed volume fraction used by CalcMixture()
  FieldCell<Scal> fc_smvfst_; // smoothed volume fraction for surface tension
  FieldFace<Scal> ff_smvf_; // interpolated fc_smvf_
  Scal diff_;  // convergence indicator
  Scal pdist_, pdistmin_; // distance to pfixed cell

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

  using Sph = typename SA::Sph;
  std::vector<Sph> sa_ss_;

  Scal nabort_; // number of abort requests, used by Reduce in checknan Run()

  std::function<void(FieldCell<typename M::Scal>&,const M&)> bgf_; // bubgen
  Scal bgt_ = -1.; // bubgen last time 

  struct Stat {
    Scal m1, m2; // mass
    Vect c1, c2;  // center of mass 
    Vect vc1, vc2;  // center of mass velocity
    Vect v1, v2;  // average velocity
    Scal dtt;  // temporary to reduce
    Scal dt;    // dt fluid 
    Scal dta;  // dt advection
    size_t iter; // iter of fluid solver
    Scal dumpt; // last dump time (rounded to nearest dtdump)
    Scal t;
    size_t step;
    size_t dumpn;
    Vect meshpos;  // mesh position
    Vect meshvel;  // mesh velocity
    Scal ekin;  /// kinetic energy
    Vect vlm; // max-norm of v-"vel"
    Vect vl2; // l2-norm of v-"vel"
    Scal p0, p1, pd; // pressure min,max
    Scal boxomm = 0.; // integral of vorticity magnitude over box
    Scal boxomm2 = 0.; // integral of vorticity magnitude over box
    Vect vomm = Vect(0); // velocity weighted by vorticity
    Scal vommw = 0; // integral of vorticity
    Scal enstr = 0; // enstrophy
    Stat()
        : m1(0), m2(0), c1(0), c2(0), vc1(0), vc2(0), v1(0), v2(0)
        , dtt(0), dt(0), dta(0), iter(0), dumpt(-1e10), step(0)
        , dumpn(0), meshpos(0), meshvel(0)
        , ekin(0)
        , vlm(0), vl2(0), p0(0), p1(0), pd(0)
    {}
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
};

// Computes vorticity of vector field.
// fcv: vector field [s]
// mf: boundary conditions for fcv
// Returns:
// fco: vorticity [i]
template <class M>
auto Hydro<M>::GetVort(const FieldCell<Vect>& fcv, 
                       const MapFace<std::shared_ptr<solver::CondFace>>& mf,
                       M& m) -> FieldCell<Vect> {
  auto ffv = solver::Interpolate(fcv, mf, m);

  auto d0 = solver::Gradient(GetComponent(ffv, 0), m);
  auto d1 = solver::Gradient(GetComponent(ffv, 1), m);
  auto d2 = solver::Gradient(GetComponent(ffv, 2), m);

  FieldCell<Vect> r(m);
  for (auto c : m.Cells()) {
    r[c][0] = d2[c][1] - d1[c][2];
    r[c][1] = d0[c][2] - d2[c][0];
    r[c][2] = d1[c][0] - d0[c][1];
  }

  return r;
}

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
  auto& fct = fcomm_;
  auto& fctv = fcom_;
  if (sem("initpois")) {
    ps_ = std::make_shared<solver::PoisSolver<M>>(GetBcSz(), m);
    m.Comm(&fc_vel_);
    fctv.Reinit(m);
  }
  for (size_t d = 0; d < M::dim; ++d) {
    std::string dn = std::to_string(d);
    if (sem("init-" + dn)) {
      fct = GetComponent(fc_vel_, d);
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
    fc_vel_ = GetVort(fctv, GetBcVz(), m);
    m.Comm(&fc_vel_);
  }
}

template <class M>
void Hydro<M>::Init() {
  auto sem = m.GetSem("init");
  if (sem("pfixed-reduce")) {
    if (auto* p = var.Double("pfixed")) {
      Vect x(var.Vect["pfixed_x"]);
      // Find cell nearest to pfixed_x
      IdxCell c = m.FindNearestCell(x);
      // TODO: add reduce minloc and remove .norm()
      pdist_ = m.GetCenter(c).dist(x) + MIdx(bi_.index).norm() * 1e-12;
      pdistmin_ = pdist_;
      m.Reduce(&pdistmin_, "min");
    }
  }

  if (sem("fields")) {
    fc_src_.Reinit(m, 0.);

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
    {
      const std::string vi = var.String["vel_init"];
      if (vi == "taylor-green") {
        for (auto i : m.AllCells()) {
          auto& v = fc_vel_[i];
          auto x = m.GetCenter(i);
          if (var.Int["dim"] == 2) {
            x[2] = 0.;
          }
          v[0] = std::sin(x[0]) * std::cos(x[1]) * std::cos(x[2]);
          v[1] = -std::cos(x[0]) * std::sin(x[1]) * std::cos(x[2]);
          v[2] = 0.;
        }
      } else if (vi == "vortex") {
        Scal pi = M_PI;
        Vect xc1(var.Vect["vort_x1"]);
        Vect xc2(var.Vect["vort_x2"]);
        Scal g1 = var.Double["vort_g1"];
        Scal g2 = var.Double["vort_g2"];
        Scal vmax = var.Double["vort_vmax"];
        std::vector<Vect> xxc = {xc1, xc2};
        std::vector<Scal> gg = {g1, g2};
        for (size_t i = 0; i < xxc.size(); ++i) {
          Vect xc = xxc[i];
          Scal g = gg[i];
          for (auto c : m.AllCells()) {
            Vect x = m.GetCenter(c);
            Scal r = xc.dist(x);
            Scal r2 = r * r;
            Vect v(0);
            v[0] = -(x[1] - xc[1]) * g / r2 * 2. * pi;
            v[1] = (x[0] - xc[0]) * g / r2 * 2. * pi;
            if (v.norm() > vmax) {
              v *= vmax / v.norm();
            }
            fc_vel_[c][0] += v[0];
            fc_vel_[c][1] += v[1];
          }
        }
      } else if (vi == "vortexring") { // XXX: use with initvort=1
        Vect xc(var.Vect["ring_c"]); // center of ring
        Vect n(var.Vect["ring_n"]);  // normal
        n /= n.norm();
        Scal om = var.Double["ring_om"];  // vorticity in crosssection
        Scal r0 = var.Double["ring_r0"];  // inner radius
        Scal r1 = var.Double["ring_r1"];  // outer radius
        Scal qr = (r1 - r0) * 0.5;   // radius
        Vect qc((r1 + r0) * 0.5, 0., 0.);   // center
        const Scal eps = 1e-10;
        for (auto c : m.AllCells()) {
          Vect x = m.GetCenter(c) - xc;
          // along axis
          Scal xn = n.dot(x);
          // along plane
          Scal xt = (x - n * xn).norm();
          // unit along plane
          Vect t = (x - n * xn) / std::max(eps, xt);
          // unit along circle 
          Vect s = n.cross(t);
          Vect q(xt, xn, 0.);
          fc_vel_[c] = ((q - qc).sqrnorm() <= sqr(qr) ? s * om : Vect(0));
        }
      } else if (vi == "vortexgauss") {
        // XXX: Bergdorf 2007, Direct numerical simulations of vortex rings
        Scal pi = M_PI;
        Scal g(var.Double["ring_gamma"]); // circulation
        Scal sig(var.Double["ring_sigma"]); // sigma
        Scal rad(var.Double["ring_r"]); // radius
        Scal nfr(var.Double["ring_noise_freq"]); // noise angular frequency
        Scal namp(var.Double["ring_noise_amp"]); // noise amp relative to r
        Scal nfr2(var.Double["ring_noise2_freq"]); // noise angular frequency
        Scal namp2(var.Double["ring_noise2_amp"]); // noise amp relative to r
        Vect xc(var.Vect["ring_c"]); // center
        Vect n(var.Vect["ring_n"]);  // normal
        n /= n.norm();
        const Scal eps = 1e-10;
        for (auto c : m.AllCells()) {
          Vect x = m.GetCenter(c) - xc;
          // along axis
          Scal xn = n.dot(x);
          // select direction 
          Vect v(0);
          v[n.abs().argmin()] = 1.;
          // unit vectors in plane
          Vect vx = n.cross(v);
          Vect vy = n.cross(vx);
          // angle
          Scal a = std::atan2(vx.dot(x), vy.dot(x));
          // noise
          xn += rad * namp * std::sin(a * nfr);
          // along plane
          Scal xt = (x - n * xn).norm();
          // noise2
          xt += rad * namp2 * std::sin(a * nfr2);
          // unit radial along plane
          Vect et = (x - n * xn) / std::max(eps, xt);
          // unit along circle 
          Vect es = n.cross(et);
          Scal s2 = sqr(xn) + sqr(xt - rad);
          Scal sig2 = sqr(sig);
          Scal om = g / (pi * sig2) * std::exp(-s2 / sig2);
          fc_vel_[c] = es * om;
        }
      } else if (vi == "rot") {
        Vect xc(var.Vect["rot_c"]);
        Scal om = var.Double["rot_om"];
        for (auto c : m.AllCells()) {
          Vect x = m.GetCenter(c);
          Vect v(0);
          v[0] = -(x[1] - xc[1]) * om * 0.5;
          v[1] = (x[0] - xc[0]) * om * 0.5;
          fc_vel_[c] = v;
        }
      } else if (vi == "grad") {
        Vect xc(var.Vect["grad_c"]);
        Vect vx(var.Vect["grad_vx"]); // gradient of vel[0]
        Vect vy(var.Vect["grad_vy"]); // gradient of vel[1]
        Vect vz(var.Vect["grad_vz"]); // gradient of vel[2]
        for (auto c : m.AllCells()) {
          Vect x = m.GetCenter(c) - xc;
          Vect v;
          v[0] = vx.dot(x);
          v[1] = vy.dot(x);
          v[2] = vz.dot(x);
          fc_vel_[c] = v;
        }
      } else if (vi == "pois" || vi == "poisy") {
        // Poiseuille with walls in y
        MIdx gs = m.GetGlobalSize(); // global mesh size
        Scal ext = var.Double["extent"]; // TODO: revise
        Vect gh = Vect(gs) * ext / gs.max(); // global domain length
        Scal pv = var.Double["poisvel"]; // centerline velocity

        for (auto i : m.AllCells()) {
          Scal y = m.GetCenter(i)[1] / gh[1];
          fc_vel_[i][0] = y * (1. - y) * 4. * pv;
        }
      } else if (vi == "poisyz") {
        // Poiseuille with walls in y and z
        // Spiga 1994: Symmetric solution for velocity in rectangular ducts
        MIdx gs = m.GetGlobalSize(); // global mesh size
        Scal ext = var.Double["extent"]; // TODO: revise
        Vect gh = Vect(gs) * ext / gs.max(); // global domain length
        Scal mu = var.Double["poismu"]; // viscosity
        Scal pg = var.Double["poisgrad"]; // pressure gradient
        int im = var.Int["poisiter"]; // depth to evaluate series
        bool wym = var.Int["poiswym"]; // wallym
        bool wyp = var.Int["poiswyp"]; // wallyp
        bool wzm = var.Int["poiswzm"]; // wallzm
        bool wzp = var.Int["poiswzp"]; // wallzp
        Scal pi = M_PI;

        Scal ly = gh[1];
        Scal lz = gh[2];

        if ((!wym && !wyp) || (!wzm && !wzp)) {
          throw std::runtime_error("poisyz: can't remove both walls");
        }

        Scal oy = 0.;
        if (!wym) {
          oy = 0.5;
          ly *= 2;
        } else if (!wyp) {
          oy = 0;
          ly *= 2;
        }

        Scal oz = 0.;
        if (!wzm) {
          oz = 0.5;
          lz *= 2;
        } else if (!wzp) {
          oz = 0;
          lz *= 2;
        }


        Scal p = sqr(ly) * pg / mu;
        Scal b = lz / ly;
        Scal k = 16. * sqr(b) / std::pow(pi, 4);

        // TODO: tests for gh[1] != gh[2]
        for (auto i : m.AllCells()) {
          Scal y = m.GetCenter(i)[1] / ly;
          y = y + oy;
          Scal z = m.GetCenter(i)[2] / lz;
          z = z + oz;
          Scal s = 0.;
          for (int iy = 1; iy < im * 2; iy += 2) {
            for (int iz = 1; iz < im * 2; iz += 2) {
              s += std::sin(iy * pi * y) *  std::sin(iz * pi * z) / 
                  (iy * iz * (sqr(b) * sqr(iy) + sqr(iz)));
            }
          }
          fc_vel_[i][0] = p * s * k;
        }
      } else if (vi == "uniform" ) {
        Vect v(var.Vect["vel"]);
        for (auto i : m.AllCells()) {
          fc_vel_[i] = v;
        }
      } else if (vi == "zero" ) {
        // nop
      } else  {
        throw std::runtime_error("Init(): unknown vel_init=" + vi);
      }
    }
    m.Comm(&fc_vel_);

    // global mesh size
    MIdx gs = m.GetGlobalSize();

    if (m.IsRoot()) {
      std::cout << "global mesh=" << gs << std::endl;
      std::cout << "surface tension dt=" << GetStDt() << std::endl;
    }

    using Dir = typename M::Dir;
    // boundary xm of global mesh
    auto gxm = [this](IdxFace i) -> bool {
      return m.GetDir(i) == Dir::i &&
          m.GetIndexFaces().GetMIdx(i)[0] == 0;
    };
    auto gxp = [this,gs](IdxFace i) -> bool {
      return m.GetDir(i) == Dir::i &&
          m.GetIndexFaces().GetMIdx(i)[0] == gs[0];
    };
    auto gym = [this](IdxFace i) -> bool {
      return m.GetDir(i) == Dir::j &&
          m.GetIndexFaces().GetMIdx(i)[1] == 0;
    };
    auto gyp = [this,gs](IdxFace i) -> bool {
      return m.GetDir(i) == Dir::j &&
          m.GetIndexFaces().GetMIdx(i)[1] == gs[1];
    };
    auto gzm = [this](IdxFace i) -> bool {
      return dim >= 3 && m.GetDir(i) == Dir::k &&
          m.GetIndexFaces().GetMIdx(i)[2] == 0;
    };
    auto gzp = [this,gs](IdxFace i) -> bool {
      return dim >= 3 && m.GetDir(i) == Dir::k &&
          m.GetIndexFaces().GetMIdx(i)[2] == gs[2];
    };
    // Set condition bc for face i on global box boundary
    // choosing proper neighbour cell id (nci)
    // Return true if on global boundary
    auto set_bc = [&](IdxFace i, std::string bc) -> bool {
      if (gxm(i) || gym(i) || gzm(i)) {
        mf_velcond_[i] = solver::Parse(bc, i, 1, m);
        return true;
      } else if (gxp(i) || gyp(i) || gzp(i)) {
        mf_velcond_[i] = solver::Parse(bc, i, 0, m);
        return true;
      } 
      return false;
    };

    // Boundary conditions for fluid 
    auto ff = m.AllFaces();
    std::vector<std::pair<std::string, std::function<bool(IdxFace)>>> pp = 
        {{"bc_xm", gxm}, {"bc_xp", gxp},
         {"bc_ym", gym}, {"bc_yp", gyp},
         {"bc_zm", gzm}, {"bc_zp", gzp}};

    for (auto p : pp) {
      if (auto bc = var.String(p.first)) {
        for (auto i : ff) {
          p.second(i) && set_bc(i, *bc);
        }
      } 
    }

    // boundary conditions for advection
    for (auto it : mf_velcond_) {
      IdxFace i = it.GetIdx();
      solver::CondFaceFluid* cb = it.GetValue().get();
      if (dynamic_cast<solver::fluid_condition::SlipWall<M>*>(cb)) {
        mf_cond_[i] = std::make_shared<solver::
            CondFaceReflect>(it.GetValue()->GetNci());
      } else {
        mf_cond_[i] = std::make_shared<solver::
            CondFaceGradFixed<Scal>>(Scal(0), it.GetValue()->GetNci());
      }
    }
    // selection boxes
    // Parameters (N>=0):
    // string boxN -- bc description
    // vect boxN_a -- lower corner
    // vect boxN_b -- upper corner
    // double boxN_vf -- inlet volume fraction
    // Check at least first nmax indices and all contiguous
    {
      int n = 0;
      const int nmax = 100;
      while (true) {
        std::string k = "box" + std::to_string(n);
        if (auto p = var.String(k)) {
          Vect a(var.Vect[k + "_a"]);
          Vect b(var.Vect[k + "_b"]);
          Scal vf = var.Double[k + "_vf"];
          Rect<Vect> r(a, b);
          for (auto i : m.AllFaces()) {
            Vect x = m.GetCenter(i);
            if (r.IsInside(x)) {
              if (set_bc(i, *p)) {
                auto b = mf_velcond_[i];
                mf_cond_[i] = std::make_shared
                    <solver::CondFaceValFixed<Scal>>(vf, b->GetNci());
              }
            }
          }
        } else if (n > nmax) { 
          break;
        }
        ++n;
      }
    }
    // selection faces
    // Parameters (N>=0):
    // string faceN -- bc description
    // vect faceN_a -- lower corner
    // vect faceN_b -- upper corner
    // double faceN_vf -- inlet volume fraction
    // Check at least first nmax indices and all contiguous
    {
      int n = -1;
      const int nmax = 100;
      while (true) {
        ++n;
        std::string k = "face" + std::to_string(n);
        if (auto p = var.String(k)) {
          // set boundary conditions on faces of box (a,b) 
          // normal to d with outer normal towards (b-a)[d]
          Vect a(var.Vect[k + "_a"]); 
          Vect b(var.Vect[k + "_b"]);
          int d(var.Int[k + "_dir"]); // direction: 0:x, 1:y, 2:z
          Scal vf = var.Double[k + "_vf"];
          Rect<Vect> r(a, b);
          Vect h = m.GetCellSize();
          auto& cb = m.GetInBlockCells();
          auto& fi = m.GetIndexFaces();
          // indices of [a,b), [begin,end)
          Vect xd(0);
          xd[d] = 1.;
          // round to faces
          a = Vect(MIdx((a + h * 0.5) / h)) * h;
          b = Vect(MIdx((b + h * 0.5) / h)) * h;
          // indices
          MIdx wa(a / h + xd * 0.5);
          MIdx wb(b / h + xd * 0.5);
          wb[d] = wa[d] + 1; // size 1 in direction d
          // direction
          MIdx wd(0);
          wd[d] = 1;
          // direction of neighbour cell
          int nci = ((b - a)[d] > 0. ? 0 : 1);
          // box of valid indices
          MIdx w0 = cb.GetBegin();
          MIdx w1 = cb.GetEnd() + wd;
          // clip (a,b) to valid indices
          wa = wa.clip(w0, w1);
          wb = wb.clip(w0, w1);
          // size of local block
          MIdx ws = wb - wa;
          if (ws.prod() == 0) {
            continue;
          }
          typename M::BlockCells bb(wa, ws);
          for (auto w : bb) {
            IdxFace f = fi.GetIdx(w, Dir(d));
            mf_velcond_[f] = solver::Parse(*p, f, nci, m);
          }
          (void) vf;
        } else if (n > nmax) { 
          break;
        }
      }
    }
    // selection spheres
    // Parameters (N>=0):
    // string sphN -- bc description ("inlet")
    // vect sphN_c -- center
    // double sphN_r0 -- radius inner
    // double sphN_r1 -- radius outer
    // double sphN_vf -- inlet volume fraction
    // vect sphN_vel -- velocity
    // Check at least first nmax indices and all contiguous
    {
      int n = 0;
      const int nmax = 100;
      while (true) {
        std::string k = "sph" + std::to_string(n);
        if (auto p = var.String(k)) {
          if (*p == "inlet") {
            Vect xc(var.Vect[k + "_c"]);
            Scal r0 = var.Double[k + "_r0"];
            Scal r1 = var.Double[k + "_r1"];
            Scal vf = var.Double[k + "_vf"];
            Vect vel(var.Vect[k + "_vel"]);
            for (auto i : m.AllFaces()) {
              Vect x = m.GetCenter(i);
              if (xc.dist(x) < r1) {
                Scal r = xc.dist(x);
                Vect v = vel * std::min(1., (r1 - r) / (r1 - r0));
                std::string a = "inlet";
                a += " " + std::to_string(v[0]);
                a += " " + std::to_string(v[1]);
                a += " " + std::to_string(v[2]);
                if (set_bc(i, a)) {
                  auto b = mf_velcond_[i];
                  mf_cond_[i] = std::make_shared
                      <solver::CondFaceValFixed<Scal>>(vf, b->GetNci());
                }
              }
            }
          } else {
            throw std::runtime_error("unknown selection sphere cond: " + *p);
          }
        } else if (n > nmax) { 
          break;
        }
        ++n;
      }
    }
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
    // Init rho, mu and force based on volume fraction
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

  if (sem("solv")) {
    // time step
    const Scal dt = var.Double["dt0"];
    st_.dt = dt;
    st_.dta = dt;
    var.Double.Set("dt", st_.dt);
    var.Double.Set("dta", st_.dta);

    // cell conditions for advection
    // (empty)
    MapCell<std::shared_ptr<solver::CondCell>> mc_cond;
    // cell conditions for fluid
    MapCell<std::shared_ptr<solver::CondCellFluid>> mc_velcond;
    {
      // Fix pressure at one cell
      if (auto* p = var.Double("pfixed")) {
        Vect x(var.Vect["pfixed_x"]);
        IdxCell c = m.FindNearestCell(x);
        if (pdist_ == pdistmin_) {
          std::cout << "pfixed bi=" << MIdx(bi_.index) 
              << " dist=" << pdist_ << std::endl;
          mc_velcond[c] = std::make_shared
              <solver::fluid_condition::GivenPressureFixed<M>>(*p);
        }
      }
    }
    // exclude cells
    {
      int n = -1;
      const int nmax = 100;
      while (true) {
        ++n;
        std::string k = "cellbox" + std::to_string(n);
        if (auto p = var.String(k)) {
          // set boundary conditions on faces of box (a,b) 
          // normal to d with outer normal towards (b-a)[d]
          Vect a(var.Vect[k + "_a"]); 
          Vect b(var.Vect[k + "_b"]);
          Rect<Vect> r(a, b);
          Vect h = m.GetCellSize();
          auto& cb = m.GetInBlockCells();
          auto& ci = m.GetIndexCells();
          // round to faces
          a = Vect(MIdx((a + h * 0.5) / h)) * h;
          b = Vect(MIdx((b + h * 0.5) / h)) * h;
          // indices
          MIdx wa(a / h);
          MIdx wb(b / h);
          // box of valid indices
          MIdx w0 = cb.GetBegin();
          MIdx w1 = cb.GetEnd();
          // clip (a,b) to valid indices
          wa = wa.clip(w0, w1);
          wb = wb.clip(w0, w1);
          // size of local block
          MIdx ws = wb - wa;
          if (ws.prod() == 0) {
            continue;
          }
          typename M::BlockCells bb(wa, ws);
          for (auto w : bb) {
            IdxCell c = ci.GetIdx(w);
            mc_velcond[c] = std::make_shared
                <solver::fluid_condition::
                GivenVelocityAndPressureFixed<M>>(Vect(0), 0.);
          }
        } else if (n > nmax) { 
          break;
        }
      }
    }
    // Init fluid solver
    {
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
            m, fc_vel_, mf_velcond_, mc_velcond, 
            &fc_rho_, &fc_mu_, &fc_force_, &ffbp_,
            &fc_src_, &fc_src_, 0., st_.dt, p));
      // TODO: check fc_src_ for Simple()

      fcbc_ = GetBcField(mf_velcond_, m);
    }

    // Init advection solver
    {
      std::string as = var.String["advection_solver"];
      if (as == "tvd") {
        auto p = std::make_shared<typename AST::Par>();
        Parse<M>(p.get(), var);
        as_.reset(new AST(
              m, fc_vf_, mf_cond_, 
              &fs_->GetVolumeFlux(solver::Layers::time_curr),
              &fc_src_, 0., st_.dta, p));
      } else if (as == "vof") {
        auto p = std::make_shared<typename ASV::Par>();
        Parse<M>(p.get(), var);
        p->dmp = std::unique_ptr<Dumper>(new Dumper(var, "dump_part_"));
        as_.reset(new ASV(
              m, fc_vf_, mf_cond_, 
              &fs_->GetVolumeFlux(solver::Layers::time_curr),
              &fc_src_, 0., st_.dta, p));
      } else {
        throw std::runtime_error("Unknown advection_solver=" + as);
      }
    }

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

    // Stat: var.Double[p] with name n
    /*
    auto on = [this](std::string n, std::string p) {
      return std::make_shared<output::OutScalFunc<Scal>>(
          n, [&,p](){ return var.Double[p]; });
    };
    */

    // Stat: *p with name n
    auto op = [this](std::string n,  Scal* p) {
      return std::make_shared<output::OutScalFunc<Scal>>(
          n, [p](){ return *p; });
    };

    if (IsRoot()) {
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
          op("ekin", &s.ekin),
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
        con.push_back(op("p0", &s.p0));
        con.push_back(op("p1", &s.p1));
        con.push_back(op("pd", &s.pd));
      }
      if (var.Int["enstrophy"]) {
        con.push_back(op("enstr", &s.enstr));
      }
      ost_ = std::make_shared<output::SerScalPlain<Scal>>(con, "stat.dat");
    }

    ParseEvents();
  }
}


template <class M>
Hydro<M>::Hydro(Vars& var, const MyBlockInfo& bi, Par& par) 
    : KernelMeshPar<M,Par>(var, bi, par)
    , dumper_(var, "dump_field_")
    , dmptraj_(var, "dump_traj_")
    , dmptrep_(var, "dump_trep_")
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
    }

    m.Reduce(&s.m1, "sum");
    m.Reduce(&s.m2, "sum");
    m.Reduce(&s.ekin, "sum");
    for (size_t d = 0; d < dim; ++d) {
      m.Reduce(&s.c1[d], "sum");
      m.Reduce(&s.c2[d], "sum");
      m.Reduce(&s.v1[d], "sum");
      m.Reduce(&s.v2[d], "sum");
    }

    if (var.Int["statvel"]) {
      s.vlm = Vect(0);
      s.vl2 = Vect(0);
      s.p0 = 1e10;
      s.p1 = -1e10;
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
        s.p0 = std::min(s.p0, fp[c]);
        s.p1 = std::max(s.p1, fp[c]);
      }
      for (size_t d = 0; d < dim; ++d) {
        m.Reduce(&s.vlm[d], "max");
        m.Reduce(&s.vl2[d], "sum");
      }
      m.Reduce(&s.p0, "min");
      m.Reduce(&s.p1, "max");
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
  }

  if (sem("reduce")) {
    Scal im1 = (s.m1 == 0 ? 0. : 1. / s.m1);
    Scal im2 = (s.m2 == 0 ? 0. : 1. / s.m2);
    s.c1 *= im1;
    s.c2 *= im2;
    s.v1 *= im1;
    s.v2 *= im2;

    if (s.vommw != 0) {
      s.vomm /= s.vommw;
    }

    // Moving mesh
    s.c1 += st_.meshpos;
    s.c2 += st_.meshpos;

    Scal dt = fs_->GetTimeStep();
    s.vc1 = (s.c1 - s.vc1) / dt;
    s.vc2 = (s.c2 - s.vc2) / dt;

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
      s.pd = s.p1 - s.p0;
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
void Hydro<M>::ParseEvents() {
  // Check at least first nmax indices and all contiguous
  int n = 0;
  const int nmax = 100;
  while (true) {
    std::string k = "ev" + std::to_string(n);
    if (auto p = var.String(k)) {
      Event e;
      std::stringstream b(*p); // buf
      b >> std::skipws;

      // format: <time> <cmd> <arg>
      b >> e.t;
      b >> e.cmd;

      char c;
      // Read first non-ws character
      b >> c;
      // Read remaining line
      std::string s;
      std::getline(b, s);
      e.arg = c + s;

      ev_.emplace(k, e);
    } else if (n > nmax) { 
      break;
    }
    ++n;
  }

  if (IsRoot()) {
    std::cout << "Found events: \n=====" << std::endl;
    for (auto p : ev_) {
      Event& e = p.second;
      std::cout << p.first << " " 
          << e.t << " " << e.cmd << " " 
          << e.arg << std::endl;
    }
    std::cout << "=====" << std::endl;
  }
}

// events: evN <time> <command>
// comamnds: set, print, setdt, setdta, vf_init
// set <type> <key> <value>
// echo string
// setdt <value>
// setdta <value>
// vf_init zero|list
template <class M>
void Hydro<M>::ExecEvents() {
  for (auto it = ev_.begin(); it != ev_.end();) {
    auto& e = it->second;
    std::string c = e.cmd;
    std::string a = e.arg;

    if (st_.t >= e.t) {
      if (IsRoot()) {
        std::cout << std::fixed << std::setprecision(8)
            << "Event at t=" << e.t << ": " 
            << c << " " << a << std::endl;
      }
      if (c == "echo") {
        if (IsRoot()) {
          std::cout << a << std::endl;
        }
      } else if (c == "set") {
        Parser p(var);
        if (IsLead()) {
          p.Run(c + " " + a);
        }
      } else {
        throw std::runtime_error("ExecEvents(): Unknown command '" + c + "'");
      }
      it = ev_.erase(it);
    } else {
      ++it;
    }
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
void Hydro<M>::CalcMixture(const FieldCell<Scal>& fc_vf0) {
  auto sem = m.GetSem("mixture");

  if (sem("init")) {
    fc_mu_.Reinit(m);
    fc_rho_.Reinit(m);
    fc_force_.Reinit(m, Vect(0));
    ffbp_.Reinit(m, 0);
    fc_smvf_ = fc_vf0;
    fc_smvfst_ = fc_vf0;
  }

  if (sem.Nested("smooth")) {
    solver::Smoothen(fc_smvf_, mf_cond_, m, var.Int["vfsmooth"]);
  }

  if (sem.Nested("smoothst")) {
    solver::Smoothen(fc_smvfst_, mf_cond_, m, var.Int["vfsmoothst"]);
  }

  if (sem("calc")) {
    auto& a = fc_smvf_;
    auto& af = ff_smvf_;  
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
    if (var.Int["enable_surftens"] && as_) { // (skip if as_ is null)
      auto af = solver::Interpolate(a, mf_cond_, m);
      auto gc = solver::Gradient(af, m); // [s]

      // zero-derivative bc for Vect
      MapFace<std::shared_ptr<solver::CondFace>> mfvz;
      for (auto it : mf_velcond_) {
        IdxFace i = it.GetIdx();
        mfvz[i] = std::make_shared<solver::
          CondFaceGradFixed<Vect>>(Vect(0), it.GetValue()->GetNci());
      }

      // zero-derivative bc for Scal
      MapFace<std::shared_ptr<solver::CondFace>> mfz;
      for (auto it : mf_velcond_) {
        IdxFace i = it.GetIdx();
        mfz[i] = std::make_shared<solver::
            CondFaceGradFixed<Scal>>(0, it.GetValue()->GetNci());
      }


      // gradient on faces
      auto gf = solver::Interpolate(gc, mfvz, m); // [i]

      // node-based gradient on faces
      if (var.Int["normalnode"]) {
        FieldNode<Vect> gn(m, Vect(0));
        FieldNode<Vect> l(m, Vect(0));
        for (auto c : m.SuCells()) {
          Vect xc = m.GetCenter(c);
          for (size_t q = 0; q < m.GetNumNeighbourNodes(c); ++q) {
            IdxNode n = m.GetNeighbourNode(c, q);
            Vect xn = m.GetNode(n);
            for (size_t d = 0; d < dim; ++d) {
              gn[n][d] += (xc[d] - xn[d] > 0. ? 1. : -1.) * a[c];
              l[n][d] += std::abs(xc[d] - xn[d]);
            }
          }
        }
        for (auto n : m.Nodes()) {
          gn[n] /= l[n];
        }
        gf.Reinit(m, Vect(0));
        for (auto f : m.Faces()) {
          for (size_t q = 0; q < m.GetNumNeighbourNodes(f); ++q) {
            IdxNode n = m.GetNeighbourNode(f, q);
            gf[f] += gn[n];
          }
          gf[f] /= m.GetNumNeighbourNodes(f);
        }
        // Zero if boundary
        for (auto it : mf_velcond_) {
          IdxFace f = it.GetIdx();
          gf[f] = Vect(0);
        }
        // zero in z if 2D
        if (var.Int["dim"] <= 2) {
          for (auto f : m.Faces()) {
            using Dir = typename M::Dir;
            if (m.GetIndexFaces().GetDir(f) == Dir::k) {
              gf[f] = Vect(0); // XXX: zero in z
            }
            gf[f][2] = 0.;
          }
        }
      }

      auto st = var.String["surftens"];
      // implementation by tensor divergence
      if (st == "div") {
        auto stdiag = var.Double["stdiag"];
        for (auto c : m.Cells()) {
          Vect r(0); 
          for (auto q : m.Nci(c)) {
            IdxFace f = m.GetNeighbourFace(c, q);
            auto g = gf[f];
            auto n = g / (g.norm() + 1e-6);  // inner normal
            // TODO: revise 1e-6
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
        auto& ast = fc_smvfst_;

        ffk_.Reinit(m, 0);
        ff_sig_.Reinit(m, 0);
        ff_sig_ = solver::Interpolate(fc_sig_, mfz, m);
        // interpolate curvature
        for (auto f : m.Faces()) {
          IdxCell cm = m.GetNeighbourCell(f, 0);
          IdxCell cp = m.GetNeighbourCell(f, 1);
          /*
          if (!IsNan(fck[cm]) && !IsNan(fck[cp])) {
            Scal wm = std::abs(ast[cp] - 0.5);
            Scal wp = std::abs(ast[cm] - 0.5);
            ffk_[f] = (fck[cm] * wm + fck[cp] * wp) / (wm + wp);
          } else if (!IsNan(fck[cm])) {
            ffk_[f] = fck[cm];
          } else {
            ffk_[f] = fck[cp];
          }
          */
          if (std::abs(ast[cm] - 0.5) < std::abs(ast[cp] - 0.5)) {
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
          Vect dm = m.GetVectToCell(f, 0);
          Vect dp = m.GetVectToCell(f, 1);
          Scal ga = (ast[cp] - ast[cm]) / (dp - dm).norm();
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
              << " vf[cm]=" << ast[cm]
              << " vf[cp]=" << ast[cp]
              << " fck[cm]=" << fck[cm]
              << " fck[cp]=" << fck[cp];
          std::cout << s.str() << std::endl;
        }
        // zero on boundaries
        for (auto it : mf_velcond_) {
          IdxFace f = it.GetIdx();
          ff_st[f] = 0.;
        }
        // contact angle on boundaries
        // angle between normal to boundary and normal to interface
        //Scal ang = var.Double["contang"] * pi / 180.;
        // TODO: controller for contact angle
        Scal k = var.Double["contangk"];
        // apply correction to curvature on boundaries
        for (auto it : mf_velcond_) {
          IdxFace fb = it.GetIdx();
          auto nci = it.GetValue()->GetNci();
          IdxCell c = m.GetNeighbourCell(fb, nci);
          const Scal th = 1e-2;
          // TODO: move th to conf
          if (ast[c] > th && ast[c] < 1. - th) { 
            for (auto q : m.Nci(c)) {
              auto f = m.GetNeighbourFace(c, q);
              ff_st[f] *= k;
            }
            /*
            // outer normal to interface (if vf=1 inside bubble)
            Vect n = -gc[c]; 
            n /= n.norm();
            // outer normal to boundary
            Vect nb = m.GetNormal(fb) * (nci == 1 ? 1. : -1.); 
            // vector along boundary
            Vect a = n - nb * nb.dot(n);
            a /= a.norm();
            // vector with angle ang to boundary normal
            Vect ns = nb * std::cos(ang) + a * std::sin(ang);
            // surface tension force
            Vect s = ns * (-fck[c] * k * std::abs(gc[c].dot(ns)));
            // set to all upwind neighbour faces
            if(0)
            std::cout 
              << " nb=" << nb 
              << " x=" << m.GetCenter(c)
              << " s=" << s
              << " gc=" << gc[c]
              << std::endl;
            //s = fck[c] * k * gc[c];
            for (auto q : m.Nci(c)) {
              auto f = m.GetNeighbourFace(c, q);
              //ff_st[f] = ffk_[f] * k * gf[f].dot(m.GetNormal(f));
              ff_st[f] *= k;
              if ((m.GetCenter(f) - m.GetCenter(c)).dot(s) < 0. ||
                  m.GetNormal(f) == m.GetNormal(fb)) {
                ff_st[f] = s.dot(m.GetNormal(f));
              }
            }
            */
          }
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
          using R = Reconst<Scal>;
          if (auto as = dynamic_cast<solver::Vof<M>*>(as_.get())) {
            auto fc_gsig = solver::Gradient(ff_sig_, m);
            auto &fcn = as->GetNormal();
            auto &fca = as->GetAlpha();
            Scal th = 1e-8;
            Vect h = m.GetCellSize();
            for (auto c : m.Cells()) {
              if (ast[c] > th && ast[c] < 1. - th) { // contains interface
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

// TODO: Test m.Dump() with pending Comm
template <class M>
void Hydro<M>::Dump(Sem& sem) {
  if (sem("dump")) {
    if (dumper_.Try(st_.t, st_.dt)) {
      if (IsRoot()) {
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
    }
  }
  if (sem("dumpwrite")) {
    // Empty stage for DumpWrite
    // TODO: revise
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
        auto add = [&v,&i,this](Scal a) {
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

  // Dump init
  if (var.Int["dumpinit"]) {
    Dump(sem);
  }

  sem.LoopBegin();

  if (sem("events")) {
    ExecEvents();
  }
  if (sem("loop-check")) {
    if (st_.t + st_.dt * 0.25 > var.Double["tmax"] || 
        int(st_.step) >= var.Int["max_step"]) {
      sem.LoopBreak();
    } else {
      if (m.IsRoot()) { 
        std::cout << std::fixed << std::setprecision(8)
            << "STEP=" << st_.step 
            << " t=" << st_.t
            << " dt=" << st_.dt
            << " ta=" << as_->GetTime()
            << " dta=" << as_->GetTimeStep()
            << std::endl;
      }
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
  }
  if (sem("abort-check")) {
    if (nabort_ != 0.) {
      if (m.IsRoot()) {
        std::cout << "nabort_ = " << nabort_ << std::endl;
      }
      sem.LoopBreak();
    }
  }
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
      auto sn = m.GetSem("iter"); // sem nested
      sn.LoopBegin();
      if (sn.Nested("iter")) {
        fs_->MakeIteration();
      }
      if (sn("reduce")) {
        diff_ = fs_->GetError();
        m.Reduce(&diff_, "max");
      }
      if (sn("report")) {
        ++st_.iter;
        var.Int["iter"] = st_.iter;
        if (IsRoot()) {
          std::cout << std::scientific << std::setprecision(16)
              << ".....iter=" << fs_->GetIter()
              << ", diff=" << diff_ << std::endl;
        }
      }
      if (sn("convcheck")) {
        assert(fs_->GetError() <= diff_);
        auto it = fs_->GetIter();
        if ((diff_ < var.Double["tol"] && (int)it >= var.Int["min_iter"]) ||
            (int)it >= var.Int["max_iter"]) {
          sn.LoopBreak();
        }
      }
      // TODO: Suspender loop hangs if (probably) Nested is last
      sn.LoopEnd();
    }
  }
  if (sem.Nested("fs-finish")) {
    fs_->FinishStep();
  }

  if (var.Int["enable_advection"]) {
    if (sem.Nested("as-steps")) {
      auto sn = m.GetSem("steps"); // sem nested
      sn.LoopBegin();
      if (sn.Nested("start")) {
        as_->StartStep();
      }
      if (sn("iter")) {
        as_->MakeIteration();
      }
      if (sn.Nested("finish")) {
        as_->FinishStep();
      }
      if (sn("report")) {
        if (IsRoot()) {
          std::cout << std::fixed << std::setprecision(8)
              << ".....adv: t=" << as_->GetTime() 
              << " dt=" << as_->GetTimeStep()
              << std::endl;
        }
      }
      if (sn("convcheck")) {
        if (as_->GetTime() >= fs_->GetTime() - 0.5 * as_->GetTimeStep()) {
          sn.LoopBreak();
        }
      }
      sn.LoopEnd();
    }
    if (var.Int["clip_vf"]) {
      if (sem("as-clip")) {
        Clip(as_->GetField(), 0., 1.);
      }
    }
    // XXX: adhoc bubble generation
    if (var.Int["enable_bubgen"]) {
      Scal t0 = var.Double["bubgen_t0"];
      Scal tper = var.Double["bubgen_per"];
      bool bg = (st_.t > t0 && st_.t - bgt_ >= tper);
      if (sem("as-bubgen-init")) {
        if (bg) {
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
      }
      if (sem("as-bubgen-apply") && bg) {
        if (bg) {
          auto& u = const_cast<FieldCell<Scal>&>(as_->GetField());
          for (auto c : m.AllCells()) {
            u[c] = std::max(u[c], fc_vf_[c]);
          }
          bgt_ = st_.t;
        }
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

  if (sem("dumpstat")) {
    if (IsRoot()) {
      ost_->Write(0., "");
    }
  }

  if (sem("inc")) {
    ++st_.step;
  }

  sem.LoopEnd();
}

