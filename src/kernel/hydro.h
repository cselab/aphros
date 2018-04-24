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

#include "geom/vect.h"
#include "geom/mesh3d.h"
#include "solver/solver.h"
#include "solver/advection.h"
#include "solver/advection_vof.h"
#include "solver/advection_tvd.h"
#include "solver/conv_diff.h"
#include "solver/fluid.h"
#include "kernel.h"
#include "kernelmesh.h"
#include "parse/vars.h"
#include "parse/interp.h"
#include "dump/output.h"
#include "dump/dumper.h"


template <class M>
class Hydro : public KernelMesh<M> {
 public:
  using KM = KernelMesh<M>;
  using Mesh = M;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  using MIdx = typename Mesh::MIdx;
  using Sem = typename Mesh::Sem;
  static constexpr size_t dim = M::dim;

  template <class T>
  using FieldCell = FieldCell<T>;
  template <class T>
  using FieldFace = FieldFace<T>;
  template <class T>
  using FieldNode = FieldNode<T>;

  Hydro(Vars& var, const MyBlockInfo& bi);
  void Run() override;
  M& GetMesh() { return m; }

 protected:
  using KM::var;
  using KM::bi_;
  using KM::m;
  using KM::IsRoot;
  using KM::IsLead;

 private:
  void Init();
  void Dump(Sem& sem);
  void CalcMixture(const FieldCell<Scal>& vf);
  // Clips v to given range, uses const_cast
  void Clip(const FieldCell<Scal>& v, Scal min, Scal max);
  void CalcStat();

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

  using AS = solver::AdvectionSolverExplicit<M>;
  using FS = solver::FluidSimple<M>;

  FieldCell<Scal> fc_mu_; // viscosity
  FieldCell<Scal> fc_rho_; // density
  FieldCell<Scal> fc_src_; // source
  FieldCell<Vect> fc_force_;  // force
  FieldCell<Vect> fc_stforce_;  // stforce cells TODO: what is st
  FieldFace<Vect> ff_stforce_;  // stforce faces
  MapFace<std::shared_ptr<solver::ConditionFace>> mf_cond_;
  MapFace<std::shared_ptr<solver::ConditionFaceFluid>> mf_velcond_;
  MultiTimer<std::string> timer_; 
  std::unique_ptr<AS> as_; // advection solver
  std::unique_ptr<FS> fs_; // fluid solver
  FieldCell<Scal> fc_velux_; // velocity
  FieldCell<Scal> fc_veluy_; 
  FieldCell<Scal> fc_veluz_; 
  FieldCell<Scal> fc_p_; // pressure
  FieldCell<Scal> fc_vf_; // volume fraction used by constructor and Dump()
  FieldCell<Vect> fc_vel_; // velocity used by constructor
  FieldCell<Scal> fc_smvf_; // smoothed volume fraction used by CalcMixture()
  Scal diff_;  // convergence indicator
  size_t step_;
  Scal pdist_, pdistmin_; // distance to pfixed cell

  struct Stat {
    Scal m1, m2; // mass
    Vect c1, c2;  // center of mass 
    Vect vc1, vc2;  // center of mass velocity
    Vect v1, v2;  // average velocity
    Scal dtt;  // temporary to reduce
    Scal dt;    // dt fluid 
    Scal dta;  // dt advection
    size_t iter;
    Scal dumpt; // last dump time (rounded to nearest dtdump)
    Scal t;
    size_t dumpn;
    Vect meshpos;  // mesh position
    Vect meshvel;  // mesh velocity
    Stat()
      : m1(0), m2(0), c1(0), c2(0), vc1(0), vc2(0), v1(0), v2(0)
      , dtt(0), dt(0), dta(0), iter(0), dumpt(-1e10), dumpn(0), meshpos(0)
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
    // Find cell nearest to pfixed_x
    Vect x(var.Vect["pfixed_x"]);
    IdxCell c = m.FindNearestCell(x);
    // TODO: add reduce minloc and remove .norm()
    pdist_ = m.GetCenter(c).dist(x) + MIdx(bi_.index).norm() * 1e-12;
    pdistmin_ = pdist_;
    m.Reduce(&pdistmin_, "min");
  }

  if (sem("fields")) {
    fc_stforce_.Reinit(m, Vect(0));
    ff_stforce_.Reinit(m, Vect(0));
    fc_src_.Reinit(m, 0.);

    // initial volume fraction
    fc_vf_.Reinit(m, 0);
    {
      const std::string vi = var.String["vf_init"];
      if (vi == "sin") {
        Vect k;
        if (auto p = var.Vect("sin_k")) {
          k = Vect(*p);
        } else {
          k = Vect(2. * M_PI);
        }

        for (auto i : m.AllCells()) {
          Vect z = m.GetCenter(i) * k;
          fc_vf_[i] = std::sin(z[0]) * std::sin(z[1]) * std::sin(z[2]);
        }
      } else if (vi == "circle") {
        const Vect c(var.Vect["circle_c"]);
        const Scal r(var.Double["circle_r"]);
        for (auto i : m.AllCells()) {
          Vect x = m.GetCenter(i);
          fc_vf_[i] = (c.dist(x) < r ? 1. : 0.);
        }
      } else if (vi == "list" ) {
        std::string fn = var.String["list_path"];
        struct P {
          Vect c;
          Scal r;
        };
        std::vector<P> pp;
        std::ifstream f(fn);
        if (!f.good()) {
          throw std::runtime_error("Can't open particle list '" + fn + "'");
        }
        // Read until eof
        while (true) {
          P p;
          // Read single particle: x y z r
          f >> p.c[0] >> p.c[1] >> p.c[2] >> p.r;
          if (f.good()) {
            pp.push_back(p);
          } else {
            break;
          }
        }
        if (IsRoot()) {
          std::cerr << "Read " 
            << pp.size() << " particles from " 
            << "'" << fn << "'" << std::endl;
        }
        // Set volume fraction to 1 inside particles
        for (auto p : pp) {
          for (auto i : m.AllCells()) {
            Vect x = m.GetCenter(i);
            if (p.c.dist(x) <= p.r) {
              fc_vf_[i] = 1.;
            }
          }
        }
      } else if (vi == "zero" ) {
        // nop
      } else {
        throw std::runtime_error("Init(): unknown vf_init=" + vi);
      }
    }

    // initial velocity
    fc_vel_.Reinit(m, Vect(0));
    {
      const std::string vi = var.String["vel_init"];
      if (vi == "taylor-green") {
        for (auto i : m.AllCells()) {
          auto& u = fc_vel_[i][0];
          auto& v = fc_vel_[i][1];
          Scal l = (2 * M_PI);
          Scal x = m.GetCenter(i)[0] * l;
          Scal y = m.GetCenter(i)[1] * l;
          u = std::cos(x) * std::sin(y);
          v = -std::sin(x) * std::cos(y);
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

    st_.iter = 0;
    var.Int.Set("iter", st_.iter);

    // TODO: Comm initial

    // global mesh size
    MIdx gs;
    {
      MIdx p(var.Int["px"], var.Int["py"], var.Int["pz"]);
      MIdx b(var.Int["bx"], var.Int["by"], var.Int["bz"]);
      MIdx bs(var.Int["bsx"], var.Int["bsy"], var.Int["bsz"]);
      gs = p * b * bs;
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

    // zero-derivative boundary conditions for advection
    for (auto it : mf_velcond_) {
      IdxFace i = it.GetIdx();
      mf_cond_[i] = std::make_shared
          <solver::ConditionFaceDerivativeFixed<Scal>>(
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
                    <solver::ConditionFaceValueFixed<Scal>>(vf, b->GetNci());
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
    step_ = 0;

    // cell conditions for advection
    // (empty)
    MapCell<std::shared_ptr<solver::ConditionCell>> mc_cond;
    // cell conditions for fluid
    MapCell<std::shared_ptr<solver::ConditionCellFluid>> mc_velcond;
    {
      // Fix pressure at one cell
      Vect x(var.Vect["pfixed_x"]);
      IdxCell c = m.FindNearestCell(x);
      Scal p = var.Double["pfixed"];
      if (pdist_ == pdistmin_) {
        std::cout << "pfixed bi=" << MIdx(bi_.index) 
            << " dist=" << pdist_ << std::endl;
        mc_velcond[c] = std::make_shared
            <solver::fluid_condition::GivenPressureFixed<Mesh>>(p);
      }
    }

    // Init fluid solver
    {
      auto p = std::make_shared<typename FS::Par>();
      p->prelax = var.Double["prelax"];
      p->vrelax = var.Double["vrelax"];
      p->rhie = var.Double["rhie"];
      p->second = var.Int["second_order"];
      // TODO: add check for inletflux_numid
      p->inletflux_numid = var.Int["inletflux_numid"];

      fs_.reset(new FS(
            m, fc_vel_, mf_velcond_, mc_velcond, 
            &fc_rho_, &fc_mu_, &fc_force_, &fc_stforce_, &ff_stforce_, 
            &fc_src_, &fc_src_, 0., dt, &timer_, p));
    }

    // Init advection solver
    {
      auto p = std::make_shared<typename AS::Par>();
      p->sharp = var.Double["sharp"];
      p->sharpo = var.Double["sharpo"];
      p->split = var.Int["split"];

      as_.reset(new AS(
            m, fc_vf_, mf_cond_, 
            &fs_->GetVolumeFlux(solver::Layers::iter_curr),
            &fc_src_, 0., dt, p));
    }

    // Output from var.Double by name
    auto on = [this](std::string n /*output-name*/,  
                          std::string p /*parameter*/) {
        return std::make_shared<output::EntryScalarFunction<Scal>>(
            n, [&,p](){ return var.Double[p]; });
      };

    // Output by pointer
    auto op = [this](std::string n /*output-name*/,  Scal* p /*pointer*/) {
        return std::make_shared<output::EntryScalarFunction<Scal>>(
            n, [p](){ return *p; });
      };


    st_.t = fs_->GetTime();
    var.Double.Set("t", st_.t);

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
Hydro<M>::Hydro(Vars& var, const MyBlockInfo& bi) 
  : KernelMesh<M>(var, bi)
  , dumper_(var)
{}

template <class M>
void Hydro<M>::CalcStat() {
  auto sem = m.GetSem("stat");

  auto& s = st_;

  if (sem("local")) {
    st_.t = fs_->GetTime();
    var.Double.Set("t", st_.t);

    auto& fa = as_->GetField();
    auto& fv = fs_->GetVelocity();

    // Store vc1 and vc2
    s.vc1 = s.c1;
    s.vc2 = s.c2;

    // mass, center, velocity
    s.m1 = 0;
    s.m2 = 0;
    s.c1 = Vect(0);
    s.c2 = Vect(0);
    s.v1 = Vect(0);
    s.v2 = Vect(0);
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
    }

    m.Reduce(&s.m1, "sum");
    m.Reduce(&s.m2, "sum");
    for (auto d = 0; d < dim; ++d) {
      m.Reduce(&s.c1[d], "sum");
      m.Reduce(&s.c2[d], "sum");
      m.Reduce(&s.v1[d], "sum");
      m.Reduce(&s.v2[d], "sum");
    }
  }

  if (sem("reduce")) {
    Scal im1 = (s.m1 == 0 ? 0. : 1. /s.m1);
    Scal im2 = (s.m2 == 0 ? 0. : 1. /s.m2);
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
  if (sem("dt-local")) {
    st_.dtt = fs_->GetAutoTimeStep();
    m.Reduce(&st_.dtt, "min");
  }
  if (sem("dt-reduce")) {
    if (st_.iter) { // TODO: revise skipping first iter
      if (auto* cfl = var.Double("cfl")) {
        st_.dt = st_.dtt * (*cfl);
        st_.dt = std::min<Scal>(st_.dt, var.Double["dtmax"]);
        fs_->SetTimeStep(st_.dt);
        var.Double["dt"] = st_.dt;
      }

      if (auto* cfla = var.Double("cfla")) {
        st_.dta = st_.dtt * (*cfla); 
        st_.dta = std::min<Scal>(st_.dta, var.Double["dtmax"]);
        as_->SetTimeStep(st_.dta);
        var.Double["dta"] = st_.dta;
      }
    }
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
        std::cout << std::defaultfloat << "Event at t=" << e.t << ": " 
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
    fc_force_.Reinit(m);
    fc_smvf_ = fc_vf0;
  }

  if (sem.Nested("smooth")) {
    solver::Smoothen(fc_smvf_, mf_cond_, m, var.Int["vfsmooth"]);
  }

  if (sem("calc")) {
    auto& a = fc_smvf_; // smoothed vf

    const Vect f(var.Vect["force"]);
    const Vect g(var.Vect["gravity"]);
    const Scal r1(var.Double["rho1"]);
    const Scal r2(var.Double["rho2"]);
    const Scal m1(var.Double["mu1"]);
    const Scal m2(var.Double["mu2"]);

    // Init density and viscosity
    for (auto i : m.AllCells()) {
      const Scal v2 = a[i];
      const Scal v1 = 1. - v2;
      fc_rho_[i] = r1 * v1 + r2 * v2;
      fc_mu_[i] = m1 * v1 + m2 * v2;
    }

    // Init force
    for (auto i : m.AllCells()) {
      Vect x = m.GetCenter(i);
      fc_force_[i] = f;
      fc_force_[i] += g * fc_rho_[i];
    }

    // Surface tension
    if (var.Int["enable_surftens"]) {
      auto af = solver::Interpolate(a, mf_cond_, m);
      auto gc = solver::Gradient(af, m);

      // zero-derivative bc for Vect
      MapFace<std::shared_ptr<solver::ConditionFace>> mfvz;
      for (auto it : mf_velcond_) {
        IdxFace i = it.GetIdx();
        mfvz[i] = std::make_shared
            <solver::ConditionFaceDerivativeFixed<Vect>>(
                Vect(0), it.GetValue()->GetNci());
      }

      fc_stforce_.Reinit(m, Vect(0));

      // surface tension in cells
      auto sig = var.Double["sigma"];
      auto stdiag = var.Double["stdiag"];
      auto gf = solver::Interpolate(gc, mfvz, m);
      for (auto c : m.Cells()) {
        Vect r(0); // result
        for (size_t e = 0; e < m.GetNumNeighbourFaces(c); ++e) {
          IdxFace f = m.GetNeighbourFace(c, e);
          auto g = gf[f];
          auto n = g / (g.norm() + 1e-6); // TODO: revise 1e-6
          auto s = m.GetOutwardSurface(c, e);
          r += s * (g.norm() * stdiag) - g * s.dot(n);
        }
        r /= m.GetVolume(c);     // div(gg/|g|) - div(|g|I)
        r[2] = 0.;
        //fc_stforce_[c] = r * sig;
        fc_force_[c] += r * sig;
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
  Dump(sem);

  sem.LoopBegin();

  if (sem("loop-check")) {
    ++step_;
    if (fs_->GetTime() >= var.Double["tmax"] ||
        step_ > var.Int["max_step"]) {
      sem.LoopBreak();
    } else if (IsRoot()) { 
      std::cerr 
          << "STEP=" << step_ 
          << " t=" << st_.t
          << " dt=" << st_.dt
          << std::endl;
    }
  }
  if (sem.Nested("mixture")) {
    CalcMixture(as_->GetField());
  }
  if (sem.Nested("stat")) {
    CalcStat();
  }
  if (sem.Nested("events")) {
    ExecEvents();
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
        if (IsRoot()) {
          std::cout << std::scientific << std::setprecision(16)
              << ".....iter=" << fs_->GetIter()
              << ", diff=" << diff_ << std::endl;
          var.Int["iter"] = st_.iter;
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
      if (sn("convcheck")) {
        if (as_->GetTime() >= fs_->GetTime() - 0.5 * as_->GetTimeStep()) {
          sn.LoopBreak();
        }
      }
      if (sn.Nested("finish")) {
        as_->FinishStep();
      }
      sn.LoopEnd();
    }
    if (var.Int["clip_vf"]) {
      if (sem("as-clip")) {
        Clip(as_->GetField(), 0., 1.);
      }
    }
  }

  Dump(sem);

  if (sem("dumpstat")) {
    if (IsRoot()) {
      ost_->Write();
    }
  }

  sem.LoopEnd();

  if (sem("fluid-timer")) {
    if (IsRoot() && var.Int["verbose_fluid_timer"]) {
      double a = 0.; // total
      auto& mt = timer_;
      for (auto e : mt.GetMap()) {
        a += e.second;
      }

      std::cout << std::fixed;
      for (auto e : mt.GetMap()) {
        auto n = e.first; // name
        auto t = e.second; // time

        std::cout 
            << n << "\n" 
            << std::setprecision(5) << t << " s = "
            << std::setprecision(3) << 100. * t / a << "%\n";
      }
      std::cout << std::endl;
    }
  }
}

template <class _M>
class HydroFactory : public KernelMeshFactory<_M> {
 public:
  using M = _M;
  using K = Hydro<M>;
  K* Make(Vars& var, const MyBlockInfo& bi) override {
    return new Hydro<M>(var, bi);
  }
};


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
