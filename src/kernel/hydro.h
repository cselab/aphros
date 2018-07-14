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

#include "geom/mesh.h"
#include "solver/solver.h"
#include "solver/advection.h"
#include "solver/vof.h"
#include "solver/tvd.h"
#include "solver/simple.h"
#include "kernel.h"
#include "kernelmesh.h"
#include "parse/vars.h"
#include "parse/interp.h"
#include "dump/output.h"
#include "dump/dumper.h"
#include "func/init_u.h"
#include "func/init_cl.h"

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
  void Dump(Sem& sem);
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

  using FS = solver::FluidSimple<M>;
  using AST = solver::AdvectionSolverExplicit<M>; // advection TVD
  using ASV = solver::Vof<M>; // advection VOF

  void Update(typename FS::Par* p) {
    p->prelax = var.Double["prelax"];
    p->vrelax = var.Double["vrelax"];
    p->rhie = var.Double["rhie"];
    p->second = var.Int["second_order"];
    p->simpler = var.Int["simpler"];
    // TODO: add check for inletflux_numid
    p->inletflux_numid = var.Int["inletflux_numid"];
    p->convsc = solver::GetConvSc(var.String["convsc"]);
    p->convdf = var.Double["convdf"];
  }
  void Update(typename AST::Par* p) {
    p->sharp = var.Double["sharp"];
    p->sharpo = var.Double["sharpo"];
    p->split = var.Int["split"];
  }
  void Update(typename ASV::Par* p) {
    p->curvgrad = var.Int["curvgrad"];
    p->part = var.Int["part"];
    p->part_verb = var.Int["part_verb"];
    p->part_relax = var.Double["part_relax"];
    p->part_h = var.Double["part_h"];
    p->part_kstr = var.Double["part_kstr"];
    p->part_kattr = var.Double["part_kattr"];
    p->part_kbend = var.Double["part_kbend"];
    p->part_bendmean = var.Int["part_bendmean"];
    p->part_dump_fr = var.Int["part_dump_fr"];
    p->part_report_fr = var.Int["part_report_fr"];
    p->part_n = var.Int["part_n"];
    p->part_k = var.Int["part_k"];
    p->part_intth = var.Double["part_intth"];
    p->poly_intth = var.Double["poly_intth"];
    p->clipth = var.Double["clipth"];
    p->dim = var.Int["dim"];
    p->dumppoly = var.Int["dumppoly"];
    p->bcc_k0 = var.Double["bcc_k0"];
    p->bcc_k1 = var.Double["bcc_k1"];
    p->bcc_t0 = var.Double["bcc_t0"];
    p->bcc_t1 = var.Double["bcc_t1"];
    p->bcc_y0 = var.Double["bcc_y0"];
    p->bcc_y1 = var.Double["bcc_y1"];
    p->part_constr = var.Int["part_constr"];
    p->part_segcirc = var.Double["part_segcirc"];
    p->part_itermax = var.Int["part_itermax"];
    p->part_tol = var.Double["part_tol"];
    p->part_np = var.Int["part_np"];
    p->part_ns = var.Int["part_ns"];
  }
  void UpdateAsPar() {
    if (auto as = dynamic_cast<AST*>(as_.get())) {
      Update(as->GetPar());
    } else if (auto as = dynamic_cast<ASV*>(as_.get())) {
      Update(as->GetPar());
    }
  }
  // Surface tension time step
  Scal GetStDt() {
    Scal sig = var.Double["sigma"];
    Scal* cflst = var.Double("cflst");
    if (cflst && sig > 0.) {
      Scal pi = M_PI;
      Scal h3 = m.GetVolume(IdxCell(0));
      Scal r1 = var.Double["rho1"];
      Scal r2 = var.Double["rho2"];
      return (*cflst) * std::sqrt(h3 * (r1 + r2) / (4. * pi * sig));
    }
    return std::numeric_limits<Scal>::max();
  }

  FieldCell<Scal> fc_mu_; // viscosity
  FieldCell<Scal> fc_rho_; // density
  FieldFace<Scal> ff_rho_; // density
  FieldCell<Scal> fc_src_; // source
  FieldCell<Vect> fc_force_;  // force 
  FieldFace<Scal> ffbp_;  // balanced force projections
  FieldFace<Scal> ff_st_;  // surface tension projections
  FieldFace<Scal> ffk_;  // curvature on faces
  MapFace<std::shared_ptr<solver::CondFace>> mf_cond_;
  MapFace<std::shared_ptr<solver::CondFaceFluid>> mf_velcond_;
  std::unique_ptr<solver::AdvectionSolver<M>> as_; // advection solver
  std::unique_ptr<FS> fs_; // fluid solver
  FieldCell<Scal> fc_velux_; // velocity
  FieldCell<Scal> fc_veluy_; 
  FieldCell<Scal> fc_veluz_; 
  FieldCell<Scal> fc_p_; // pressure used by Dump()
  FieldCell<Scal> fc_vf_; // volume fraction used by constructor and Dump()
  FieldCell<Scal> fccl_; // color 
  FieldCell<Vect> fc_vel_; // velocity used by constructor
  FieldCell<Scal> fc_smvf_; // smoothed volume fraction used by CalcMixture()
  FieldCell<Scal> fc_smvfst_; // smoothed volume fraction for surface tension
  FieldFace<Scal> ff_smvf_; // interpolated fc_smvf_
  Scal diff_;  // convergence indicator
  Scal pdist_, pdistmin_; // distance to pfixed cell

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
    Stat()
        : m1(0), m2(0), c1(0), c2(0), vc1(0), vc2(0), v1(0), v2(0)
        , dtt(0), dt(0), dta(0), iter(0), dumpt(-1e10), step(0)
        , dumpn(0), meshpos(0)
        , ekin(0)
    {}
  };
  Stat st_;
  std::shared_ptr<output::Session> ost_; // output stat
  Dumper dumper_;
};

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

    // initial color
    fccl_.Reinit(m, 0);
    auto icl = CreateInitCl<M>(var, m.IsRoot());
    icl(fccl_, fc_vf_, m);
    m.Comm(&fccl_);

    // initial velocity
    fc_vel_.Reinit(m, Vect(0));
    {
      const std::string vi = var.String["vel_init"];
      if (vi == "taylor-green") {
        for (auto i : m.AllCells()) {
          auto& v = fc_vel_[i];
          auto x = m.GetCenter(i);
          if (dim == 2) {
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
      } else if (vi == "pois") {
        Scal pv = var.Double["poisvel"];
        for (auto i : m.AllCells()) {
          Vect x = m.GetCenter(i);
          fc_vel_[i][0] = x[1] * (1. - x[1]) * 4. * pv;
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
    MIdx gs;
    {
      MIdx p(var.Int["px"], var.Int["py"], var.Int["pz"]);
      MIdx b(var.Int["bx"], var.Int["by"], var.Int["bz"]);
      MIdx bs(var.Int["bsx"], var.Int["bsy"], var.Int["bsz"]);
      gs = p * b * bs;
    }

    if (m.IsRoot()) {
      std::cout << "global mesh=" << gs << std::endl;
      std::cout << "surface tension dt=" << GetStDt() << std::endl;
    }

    using Dir = typename M::Dir;
    // boundary xm of global mesh
    auto gxm = [this](IdxFace i) -> bool {
      return m.GetDir(i) == Dir::i &&
          m.GetBlockFaces().GetMIdx(i)[0] == 0;
    };
    auto gxp = [this,gs](IdxFace i) -> bool {
      return m.GetDir(i) == Dir::i &&
          m.GetBlockFaces().GetMIdx(i)[0] == gs[0];
    };
    auto gym = [this](IdxFace i) -> bool {
      return m.GetDir(i) == Dir::j &&
          m.GetBlockFaces().GetMIdx(i)[1] == 0;
    };
    auto gyp = [this,gs](IdxFace i) -> bool {
      return m.GetDir(i) == Dir::j &&
          m.GetBlockFaces().GetMIdx(i)[1] == gs[1];
    };
    auto gzm = [this](IdxFace i) -> bool {
      return dim >= 3 && m.GetDir(i) == Dir::k &&
          m.GetBlockFaces().GetMIdx(i)[2] == 0;
    };
    auto gzp = [this,gs](IdxFace i) -> bool {
      return dim >= 3 && m.GetDir(i) == Dir::k &&
          m.GetBlockFaces().GetMIdx(i)[2] == gs[2];
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
    // any boundary of global mesh
    auto gb = [&](IdxFace i) -> bool {
      return gxm(i) || gxp(i) || gym(i) || gyp(i) || gzm(i) || gzp(i);
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

    // init with zero-derivative boundary conditions for advection
    for (auto it : mf_velcond_) {
      IdxFace i = it.GetIdx();
      mf_cond_[i] = std::make_shared
          <solver::CondFaceGradFixed<Scal>>(
              Scal(0), it.GetValue()->GetNci());
    }

    // Set bc with selection boxes
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
            if (r.IsInside(m.GetCenter(i))) {
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
  }

  if (sem.Nested("smooth")) {
    solver::Smoothen(fc_vf_, mf_cond_, m, var.Int["vf_init_sm"]);
  }

  if (sem.Nested("mixture")) {
    // Init rho, mu and force based on volume fraction
    CalcMixture(fc_vf_);
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
    // Init fluid solver
    {
      auto p = std::make_shared<typename FS::Par>();
      Update(p.get());

      fs_.reset(new FS(
            m, fc_vel_, mf_velcond_, mc_velcond, 
            &fc_rho_, &fc_mu_, &fc_force_, &ffbp_,
            &fc_src_, &fc_src_, 0., st_.dt, p));
    }

    // Init advection solver
    {
      std::string as = var.String["advection_solver"];
      if (as == "tvd") {
        auto p = std::make_shared<typename AST::Par>();
        Update(p.get());
        as_.reset(new AST(
              m, fc_vf_, mf_cond_, 
              &fs_->GetVolumeFlux(solver::Layers::time_curr),
              &fc_src_, 0., st_.dta, p));
      } else if (as == "vof") {
        auto p = std::make_shared<typename ASV::Par>();
        Update(p.get());
        p->dmp = std::unique_ptr<Dumper>(new Dumper(var, "dump_part_"));
        as_.reset(new ASV(
              m, fc_vf_, mf_cond_, 
              &fs_->GetVolumeFlux(solver::Layers::time_curr),
              &fc_src_, 0., st_.dta, p));
      } else {
        throw std::runtime_error("Unknown advection_solver=" + as);
      }
    }

    st_.iter = 0;
    var.Int.Set("iter", st_.iter);

    st_.t = fs_->GetTime();
    var.Double.Set("t", st_.t);

    // Stat: var.Double[p] with name n
    auto on = [this](std::string n, std::string p) {
      return std::make_shared<output::EntryScalarFunction<Scal>>(
          n, [&,p](){ return var.Double[p]; });
    };

    // Stat: *p with name n
    auto op = [this](std::string n,  Scal* p) {
      return std::make_shared<output::EntryScalarFunction<Scal>>(
          n, [p](){ return *p; });
    };

    if (IsRoot()) {
      auto& s = st_;
      output::Content con = {
          op("t", &s.t),
          std::make_shared<output::EntryScalarFunction<int>>(
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
          std::make_shared<output::EntryScalarFunction<Scal>>(
              "meshposx", [this](){ return st_.meshpos[0]; }),
          std::make_shared<output::EntryScalarFunction<Scal>>(
              "meshvelx", [this](){ return st_.meshvel[0]; }),
      };
      ost_ = std::make_shared<
          output::SessionPlainScalar<Scal>>(con, "stat.dat");
    }

    ParseEvents();
  }
}


template <class M>
Hydro<M>::Hydro(Vars& var, const MyBlockInfo& bi, Par& par) 
  : KernelMeshPar<M,Par>(var, bi, par)
  , dumper_(var, "dump_field_")
{}

template <class M>
void Hydro<M>::CalcStat() {
  auto sem = m.GetSem("stat");

  auto& s = st_;

  if (sem("local")) {
    auto& fa = as_->GetField();
    auto& fv = fs_->GetVelocity();

    // check abort TODO: revise,move
    for (auto c : m.Cells()) {
      if (fv[c].sqrnorm() > sqr(var.Double["abortvel"])) {
        std::stringstream g;
        g << "abortvel exceeded at x=" << m.GetCenter(c);
        throw std::runtime_error(g.str());
      }
    }

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
    for (auto d = 0; d < dim; ++d) {
      m.Reduce(&s.c1[d], "sum");
      m.Reduce(&s.c2[d], "sum");
      m.Reduce(&s.v1[d], "sum");
      m.Reduce(&s.v2[d], "sum");
    }
  }

  if (sem("reduce")) {
    Scal im1 = (s.m1 == 0 ? 0. : 1. / s.m1);
    Scal im2 = (s.m2 == 0 ? 0. : 1. / s.m2);
    s.c1 *= im1;
    s.c2 *= im2;
    s.v1 *= im1;
    s.v2 *= im2;

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
      } else {
        throw std::runtime_error("Unknown meshvel_auto=" + *s);
      }
      Vect mask(var.Vect["meshvel_mask"]); // components 0 or 1
      v *= mask;
      double w = var.Double["meshvel_weight"];
      Vect vp = fs_->par->meshvel;
      fs_->par->meshvel = v * w + vp * (1. - w);

      st_.meshvel = fs_->par->meshvel;
      st_.meshpos += st_.meshvel * st_.dt;
    }
  }
}

template <class M>
void Hydro<M>::CalcDt() {
  auto sem = m.GetSem("dt");
  auto& s = st_;

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
        Interp p(var);
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
    ff_rho_.Reinit(m);
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

    // Init density and viscosity
    for (auto f : m.AllFaces()) {
      const Scal v2 = af[f];
      const Scal v1 = 1. - v2;
      ff_rho_[f] = r1 * v1 + r2 * v2;
    }

    // Append gravity to force
    for (auto f : m.AllFaces()) {
      Vect n = m.GetNormal(f);
      ffbp_[f] += force.dot(n);
      ffbp_[f] += grav.dot(n) * ff_rho_[f];
    }

    // Surface tension
    if (var.Int["enable_surftens"] && as_) { // (skip if as_ is null)
      auto af = solver::Interpolate(a, mf_cond_, m);
      auto gc = solver::Gradient(af, m); // [s]

      // zero-derivative bc for Vect
      MapFace<std::shared_ptr<solver::CondFace>> mfvz;
      for (auto it : mf_velcond_) {
        IdxFace i = it.GetIdx();
        mfvz[i] = std::make_shared
            <solver::CondFaceGradFixed<Vect>>(
                Vect(0), it.GetValue()->GetNci());
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
            if (m.GetBlockFaces().GetDir(f) == Dir::k) {
              gf[f] = Vect(0); // XXX: zero in z
            }
            gf[f][2] = 0.;
          }
        }
      }

      auto sig = var.Double["sigma"];
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
          fc_force_[c] += r * sig;
        }
      } else if (st == "kn") {  // curvature * normal
        auto& fck = as_->GetCurv(); // [a]
        ff_st_.Reinit(m, 0);
        auto& ast = fc_smvfst_;

        ffk_.Reinit(m, 0);
        // interpolate curvature
        for (auto f : m.Faces()) {
          IdxCell cm = m.GetNeighbourCell(f, 0);
          IdxCell cp = m.GetNeighbourCell(f, 1);
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
        for (auto f : m.Faces()) {
          IdxCell cm = m.GetNeighbourCell(f, 0);
          IdxCell cp = m.GetNeighbourCell(f, 1);
          Vect dm = m.GetVectToCell(f, 0);
          Vect dp = m.GetVectToCell(f, 1);
          Scal ga = (ast[cp] - ast[cm]) / (dp - dm).norm();
          if (ga != 0.) {
            if (IsNan(ffk_[f])) {
              std::stringstream s;
              s << "Nan curvature at x=" << m.GetCenter(f)
                  << " vf[cm]=" << ast[cm]
                  << " vf[cp]=" << ast[cp]
                  << " fck[cm]=" << fck[cm]
                  << " fck[cp]=" << fck[cp];
              throw std::runtime_error(s.str());
            }
            ff_st_[f] += ga * ffk_[f] * sig;
          }
        }
        // zero on boundaries
        for (auto it : mf_velcond_) {
          IdxFace f = it.GetIdx();
          ff_st_[f] = 0.;
        }
        // contact angle on boundaries
        Scal pi = M_PI;
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
              ff_st_[f] *= k;
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
              //ff_st_[f] = ffk_[f] * k * gf[f].dot(m.GetNormal(f));
              ff_st_[f] *= k;
              if ((m.GetCenter(f) - m.GetCenter(c)).dot(s) < 0. ||
                  m.GetNormal(f) == m.GetNormal(fb)) {
                ff_st_[f] = s.dot(m.GetNormal(f));
              }
            }
            */
          }
        }

        // Surface tension decay between x0 and x1 
        // XXX: adhoc TODO: revise
        const Scal x0 = var.Double["zerostx0"];
        const Scal x1 = var.Double["zerostx1"];
        // Append to force
        for (auto f : m.Faces()) {
          Scal x = m.GetCenter(f)[0];
          if (x > x0) {
            ff_st_[f] *= std::max(0., (x1 - x) / (x1 - x0));
          }
        }

        // Append to force
        for (auto f : m.Faces()) {
          ffbp_[f] += ff_st_[f];
        }
      } else {
        throw std::runtime_error("Unknown surftens=" + st);
      }
    }

    // zero force in z if 2D
    if (var.Int["dim"] <= 2) {
      for (auto f : m.Faces()) {
        using Dir = typename M::Dir;
        if (m.GetBlockFaces().GetDir(f) == Dir::k) {
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
      fc_velux_ = GetComponent(fs_->GetVelocity(), 0);
      m.Dump(&fc_velux_, "vx");
      fc_veluy_ = GetComponent(fs_->GetVelocity(), 1);
      m.Dump(&fc_veluy_, "vy");
      fc_veluz_ = GetComponent(fs_->GetVelocity(), 2);
      m.Dump(&fc_veluz_, "vz");
      fc_p_ = fs_->GetPressure();
      m.Dump(&fc_p_, "p"); 
      fc_vf_ = as_->GetField();
      m.Dump(&fc_vf_, "vf"); 
      if (IsRoot()) {
        dumper_.Report();
      }
    }
  }
  if (sem("dumpwrite")) {
    // Empty stage for DumpWrite
    // TODO: revise
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
        st_.step >= var.Int["max_step"]) {
      sem.LoopBreak();
    } else {
      if (IsRoot()) { 
        std::cout << std::fixed << std::setprecision(8)
            << "STEP=" << st_.step 
            << " t=" << st_.t
            << " dt=" << st_.dt
            << " ta=" << as_->GetTime()
            << " dta=" << as_->GetTimeStep()
            << std::endl;
      }
    }
  }
  if (sem("updatepar")) {
    Update(fs_->GetPar());
    UpdateAsPar();
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
        if ((diff_ < var.Double["tol"] && it >= var.Int["min_iter"]) ||
            it >= var.Int["max_iter"]) {
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
  }

  if (sem.Nested("dt")) {
    CalcDt(); // must be after CalcStat to keep dt for moving mesh velocity 
  }

  Dump(sem);

  if (sem("dumpstat")) {
    if (IsRoot()) {
      ost_->Write();
    }
  }

  if (sem("inc")) {
    ++st_.step;
  }

  sem.LoopEnd();
}

// Dependencies and interfaces:
// 
// Instances: Mesh=MeshStructured, Kernel=Hydro, Distr=Cubism
//
// Interfaces:
// - Mesh: cell connectivity, geometry.
//   Must be a template argument for performance.
//   Aware of: FieldCell, IdxCell
// - Kernel: Run(), Comm(), Reduce(), ReadBuffer(), WriteBuffer(), GetComm()
//   Aware of: Mesh (template), FieldCell, IdxCell, BlockInfo
// - KernelFactory: Creates Kernel
//   Aware of: Kernel, Mesh (template)
// - Distr: Takes KernelFactory, instantiates Kernel for every block,
//   reads GetComm(), implements communication
//  
// Requirements:
// - Cubism should not depend on Hydro (only on Kernel and KernelFactory)
// - Hydro should not depend on Cubism (only on Distr)
// - Hydro is Kernel, HydroFactory is KernelFactory
// - Ideally, Hydro doesn't depend on implementation of Mesh
// - Cubism is Distr
// - HydroFactory creates Hydro for a block described in Cubism-independent way
// - Hydro accepts Mesh initalized by HydroFactory 
//   and is not aware of Distr,
//   but it (as well as other mesh-related functions like Gradient())
//   aware of Mesh which has some distributed computing primitives 
//   like Comm() and Reduce() (also Solve() to solve linear systems
//   but that should be avoided)
// - Distributed primivites in Mesh are not dependent on Cubism or even Distr.
//   Typically, it is just a list of fields for communication.
// - Why doesn't Mesh depend on Distr?
//   This way Mesh is simpler and has less dependencies.
//   Possible drawbacks: information from MPI is not available.
//   All routines (e.g. solve linear system) rely on rigid interface
//   and need to be performed from the outside
//   (though this is necessary for any blocking setup).
// - Interface of Hydro: 
//     vector<FieldCell<Scal>*> GetComm()
//     M& GetMesh()
// - Interface of Cubism:
//     void Step()
//   Cubism (and any Distr) acts like a framework.
//   But every Distr has a way to access and communicate data in blocks
//   (in Cubism, it is BlockLab, Block_t and BlockInfo)
// - The question is how to organize interaction between Distr and Kernel
//   (in particular, Cubism and Hydro). Options:
//   1) sufficient interface of Kernel
//   How: same WriteBuffer() and ReadBuffer() but with some generic block buffer.
//   Cons: too rigid and may require data copying
//   2) some entity aware of both Hydro (implementation of Kernel)
//   and Cubsm (implementation of Distr), different for every pair.
//   How: visitor?
//   3) make Cubism aware of Hydro or at least Kernel<MeshStructured>
//
