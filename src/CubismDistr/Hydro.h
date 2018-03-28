#pragma once

#include <cassert>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <ctime>
#include <array>
#include <list>
#include <chrono>
#include <thread>
#include <mpi.h>
#include <stdexcept>

#include "hydro/vect.hpp"
#include "hydro/mesh3d.hpp"
#include "hydro/solver.hpp"
#include "hydro/advection.hpp"
#include "hydro/conv_diff.hpp"
#include "hydro/fluid.hpp"
#include "hydro/output.hpp"
#include "Kernel.h"
#include "KernelMesh.h"
#include "Vars.h"


template <class M>
class Hydro : public KernelMesh<M> {
 public:
  using KM = KernelMesh<M>;
  using Mesh = M;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  using MIdx = typename Mesh::MIdx;
  using Rect = geom::Rect<Vect>;
  using IdxCell = geom::IdxCell;
  using IdxFace = geom::IdxFace;
  using IdxNode = geom::IdxNode;
  using Sem = typename Mesh::Sem;
  static constexpr size_t dim = M::dim;

  template <class T>
  using FieldCell = geom::FieldCell<T>;
  template <class T>
  using FieldFace = geom::FieldFace<T>;
  template <class T>
  using FieldNode = geom::FieldNode<T>;

  Hydro(Vars& par, const MyBlockInfo& bi);
  void Run() override;
  M& GetMesh() { return m; }

 protected:
  using KM::par;
  using KM::bi_;
  using KM::m;

 private:
  void Init();
  void Dump(Sem& sem);
  void CalcMixture(const FieldCell<Scal>& vf);
  void CalcStat();

  using AS = solver::AdvectionSolverExplicit<M, FieldFace<Scal>>;
  using FS = solver::FluidSimple<M>;
  FieldCell<Scal> fc_mu_; // viscosity
  FieldCell<Scal> fc_rho_; // density
  FieldCell<Scal> fc_src_; // source
  FieldCell<Vect> fc_force_;  // force
  FieldCell<Vect> fc_stforce_;  // stforce cells TODO: what is st
  FieldFace<Vect> ff_stforce_;  // stforce faces
  geom::MapFace<std::shared_ptr<solver::ConditionFace>> mf_cond_;
  geom::MapFace<std::shared_ptr<solver::ConditionFaceFluid>> mf_velcond_;
  MultiTimer<std::string> timer_; 
  std::shared_ptr<const solver::LinearSolverFactory> p_lsf_; // linear solver factory
  std::unique_ptr<AS> as_; // advection solver
  std::unique_ptr<FS> fs_; // fluid solver
  FieldCell<Scal> fc_velux_; // velocity
  FieldCell<Scal> fc_veluy_; 
  FieldCell<Scal> fc_veluz_; 
  FieldCell<Scal> fc_p_; // pressure
  FieldCell<Scal> fc_vf_; // volume fraction
  Scal diff_;  // convergence indicator
  Scal tol_;  // convergence tolerance
  size_t step_;
  bool broot_;
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
};

template <class M>
void Hydro<M>::Init() {
  auto sem = m.GetSem("init");
  if (sem("a")) {
    broot_ = (m.GetInBlockCells().GetBegin() == MIdx(0));

    st_.iter = 0;
    par.Int.Set("iter", st_.iter);

    // TODO: Comm initial

    // global mesh size
    MIdx gs;
    {
      MIdx p(par.Int["px"], par.Int["py"], par.Int["pz"]);
      MIdx b(par.Int["bx"], par.Int["by"], par.Int["bz"]);
      MIdx bs(par.Int["bsx"], par.Int["bsy"], par.Int["bsz"]);
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
    // any boundary of global mesh
    auto gb = [&](IdxFace i) -> bool {
      return gxm(i) || gxp(i) || gym(i) || gyp(i) || gzm(i) || gzp(i);
    };

    // Boundary conditions for fluid
    if (auto p = par.String("bc_xm")) {
      for (auto i : m.Faces()) {
        gxm(i) && (mf_velcond_[i] = solver::Parse(*p, i, m));
      }
    } 
    if (auto p = par.String("bc_xp")) {
      for (auto i : m.Faces()) {
        gxp(i) && (mf_velcond_[i] = solver::Parse(*p, i, m));
      }
    } 
    if (auto p = par.String("bc_ym")) {
      for (auto i : m.Faces()) {
        gym(i) && (mf_velcond_[i] = solver::Parse(*p, i, m));
      }
    } 
    if (auto p = par.String("bc_yp")) {
      for (auto i : m.Faces()) {
        gyp(i) && (mf_velcond_[i] = solver::Parse(*p, i, m));
      }
    } 
    if (auto p = par.String("bc_zm")) {
      for (auto i : m.Faces()) {
        gzm(i) && (mf_velcond_[i] = solver::Parse(*p, i, m));
      }
    } 
    if (auto p = par.String("bc_zp")) {
      for (auto i : m.Faces()) {
        gzp(i) && (mf_velcond_[i] = solver::Parse(*p, i, m));
      }
    } 

    // zero-derivative boundary conditions for advection
    for (auto it : mf_velcond_) {
      IdxFace i = it.GetIdx();
      mf_cond_[i] = std::make_shared
          <solver::ConditionFaceDerivativeFixed<Scal>>(Scal(0));
    }
    
    // Set boundaries intersecting blocks
    {
      int n = 0;
      while (true) {
        std::string s = "box" + std::to_string(n);
        if (auto p = par.String(s)) {
          Vect a(par.Vect[s + "_a"]);
          Vect b(par.Vect[s + "_b"]);
          Scal vf = par.Double[s + "_vf"];
          geom::Rect<Vect> r(a, b);
          size_t q = 0;
          for (auto i : m.Faces()) {
            if (gb(i) && r.IsInside(m.GetCenter(i))) {
                mf_velcond_[i] = solver::Parse(*p, i, m);
                mf_cond_[i] = std::make_shared
                    <solver::ConditionFaceValueFixed<Scal>>(vf);
                ++q;
            }
          }
          if (broot_) {
            std::cout << "bc: " 
                << s << ": " 
                << (*p) << ", "
                << q << " faces"
                << std::endl;
          }
          ++n;
        } else {
          break;
        }
      }
    }


    {
      // Fix pressure at one cell
      Vect x(par.Vect["pfixed_x"]);
      IdxCell c = m.FindNearestCell(x);
      // TODO: add reduce minlocal and remove 1e-16
      pdist_ = m.GetCenter(c).dist(x) + MIdx(bi_.index).norm() * 1e-12;
      pdistmin_ = pdist_;
      m.Reduce(&pdistmin_, "min");
    }
  }

  if (sem("b")) {
    // cell conditions for advection
    // (empty)
    geom::MapCell<std::shared_ptr<solver::ConditionCell>> mc_cond;
    // cell conditions for fluid
    geom::MapCell<std::shared_ptr<solver::ConditionCellFluid>> mc_velcond;

    {
      // Fix pressure at one cell
      Vect x(par.Vect["pfixed_x"]);
      IdxCell c = m.FindNearestCell(x);
      Scal p = par.Double["pfixed"];
      if (pdist_ == pdistmin_) {
        std::cout << "pfixed bi=" << MIdx(bi_.index) 
            << " dist=" << pdist_ << std::endl;
        mc_velcond[c] = std::make_shared
            <solver::fluid_condition::GivenPressureFixed<Mesh>>(p);
      }
    }

    // time step
    const Scal dt = par.Double["dt0"];
    st_.dt = dt;
    st_.dta = dt;
    par.Double.Set("dt", st_.dt);
    par.Double.Set("dta", st_.dta);
    step_ = 0;

    Scal prelax = par.Double["prelax"];
    Scal vrelax = par.Double["vrelax"];
    Scal rhie = par.Double["rhie"];
    tol_ = par.Double["tol"];
    int max_iter = par.Int["max_iter"];
    bool so = par.Int["second_order"];

    fc_stforce_.Reinit(m, Vect(0));
    ff_stforce_.Reinit(m, Vect(0));
    fc_src_.Reinit(m, 0.);

    p_lsf_ = std::make_shared<const solver::LinearSolverFactory>(
          std::make_shared<const solver::LuDecompositionFactory>());

    // initial volume fraction
    FieldCell<Scal> fc_u(m);
    {
      const std::string vi = par.String["vf_init"];
      if (vi == "sin") {
        Vect k;
        if (auto p = par.Vect("sin_k")) {
          k = Vect(*p);
        } else {
          k = Vect(2. * M_PI);
        }

        for (auto i : m.AllCells()) {
          Vect z = m.GetCenter(i) * k;
          fc_u[i] = std::sin(z[0]) * std::sin(z[1]) * std::sin(z[2]);
        }
      } else if (vi == "circle") {
        const Vect c(par.Vect["circle_c"]);
        const Scal r(par.Double["circle_r"]);
        for (auto i : m.AllCells()) {
          Vect x = m.GetCenter(i);
          fc_u[i] = (c.dist(x) < r ? 1. : 0.);
        }
      } else if (vi == "list" ) {
        std::string fn = par.String["list_path"];
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
        if (broot_) {
          std::cerr << "Read " 
            << pp.size() << " particles from " 
            << "'" << fn << "'" << std::endl;
        }
        // Set volume fraction to 1 inside particles
        for (auto p : pp) {
          for (auto i : m.AllCells()) {
            Vect x = m.GetCenter(i);
            if (p.c.dist(x) <= p.r) {
              fc_u[i] = 1.;
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
    FieldCell<Vect> fc_vel(m, Vect(0));
    {
      const std::string vi = par.String["vel_init"];
      if (vi == "taylor-green") {
        for (auto i : m.AllCells()) {
          auto& u = fc_vel[i][0];
          auto& v = fc_vel[i][1];
          Scal l = (2 * M_PI);
          Scal x = m.GetCenter(i)[0] * l;
          Scal y = m.GetCenter(i)[1] * l;
          u = std::cos(x) * std::sin(y);
          v = -std::sin(x) * std::cos(y);
        }
      } else if (vi == "pois") {
        Scal pv = par.Double["poisvel"];
        for (auto i : m.AllCells()) {
          Vect x = m.GetCenter(i);
          fc_vel[i][0] = x[1] * (1. - x[1]) * 4. * pv;
        }
      } else if (vi == "zero" ) {
        // nop
      } else  {
        throw std::runtime_error("Init(): unknown vel_init=" + vi);
      }
    }

    // Init rho, mu and force based on volume fraction
    CalcMixture(fc_u);

    // Init fluid solver
    fs_.reset(new FS(
          m, fc_vel, 
          mf_velcond_, mc_velcond, 
          vrelax, prelax, rhie,
          &fc_rho_, &fc_mu_, 
          &fc_force_, &fc_stforce_, &ff_stforce_, 
          &fc_src_, &fc_src_,
          0., dt,
          *p_lsf_, *p_lsf_,
          tol_, max_iter, 
          &timer_, 
          so, false, false, 0., Vect(0)
          ));

    // Init advection solver
    as_.reset(new AS(
          m, fc_u, mf_cond_, 
          &fs_->GetVolumeFlux(solver::Layers::iter_curr),
          &fc_src_, 0., dt));

    // Output from par.Double by name
    auto on = [this](std::string n /*output-name*/,  
                          std::string p /*parameter*/) {
        return std::make_shared<output::EntryScalarFunction<Scal>>(
            n, [&,p](){ return par.Double[p]; });
      };

    // Output by pointer
    auto op = [this](std::string n /*output-name*/,  Scal* p /*pointer*/) {
        return std::make_shared<output::EntryScalarFunction<Scal>>(
            n, [p](){ return *p; });
      };


    st_.t = fs_->GetTime();
    par.Double.Set("t", st_.t);

    if (broot_) {
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
  }
}


// TODO: move construction to Run()
template <class M>
Hydro<M>::Hydro(Vars& par, const MyBlockInfo& bi) 
  : KernelMesh<M>(par, bi)
{
}

template <class M>
void Hydro<M>::CalcStat() {
  auto sem = m.GetSem("stat");

  auto& s = st_;

  if (sem("local")) {
    st_.t = fs_->GetTime();
    par.Double.Set("t", st_.t);

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

    if (std::string* s = par.String("meshvel_auto")) {
      Vect v(0);
      if (*s == "v") {
        v = st_.v2;
      } else if (*s == "vc") {
        v = st_.vc2;
      } else {
        throw std::runtime_error("Unknown meshvel_auto=" + *s);
      }
      Vect mask(par.Vect["meshvel_mask"]); // components 0 or 1
      v *= mask;
      double w = par.Double["meshvel_weight"];
      Vect vp = fs_->GetMeshVel();
      fs_->SetMeshVel(v * w + vp * (1. - w));

      st_.meshvel = fs_->GetMeshVel();
      st_.meshpos += st_.meshvel * st_.dt;
    }
  }
  if (sem("dta")) {
    st_.dtt = fs_->GetAutoTimeStep();
    m.Reduce(&st_.dtt, "min");
  }
  if (sem("dta-reduce")) {
    if (st_.iter) { // TODO: revise skipping first iter
      if (auto* cfl = par.Double("cfl")) {
        st_.dt = st_.dtt * (*cfl);
        st_.dt = std::min<Scal>(st_.dt, par.Double["dtmax"]);
        fs_->SetTimeStep(st_.dt);
        par.Double["dt"] = st_.dt;
      }

      if (auto* cfla = par.Double("cfla")) {
        st_.dta = st_.dtt * (*cfla); 
        st_.dta = std::min<Scal>(st_.dta, par.Double["dtmax"]);
        as_->SetTimeStep(st_.dta);
        par.Double["dta"] = st_.dta;
      }
    }
  }
}

template <class M>
void Hydro<M>::CalcMixture(const FieldCell<Scal>& fc_vf) {
  fc_mu_.Reinit(m);
  fc_rho_.Reinit(m);
  fc_force_.Reinit(m);
  
  const Vect f(par.Vect["force"]);
  const Vect g(par.Vect["gravity"]);
  const Scal r1(par.Double["rho1"]);
  const Scal r2(par.Double["rho2"]);
  const Scal m1(par.Double["mu1"]);
  const Scal m2(par.Double["mu2"]);

  // Init density and viscosity
  for (auto i : m.AllCells()) {
    const Scal v2 = fc_vf[i];
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
  if (par.Int["enable_surftens"]) {
    auto a = fc_vf;

    auto af = solver::Interpolate(a, mf_cond_, m);
    auto gc = solver::Gradient(af, m);


    // zero-derivative bc for Vect
    geom::MapFace<std::shared_ptr<solver::ConditionFace>> mfvz;
    for (auto it : mf_velcond_) {
      IdxFace i = it.GetIdx();
      mfvz[i] = std::make_shared
          <solver::ConditionFaceDerivativeFixed<Vect>>(Vect(0));
    }

    // surface tension in cells
    auto sig = par.Double["sigma"];
    auto gf = solver::Interpolate(gc, mfvz, m);
    for (auto c : m.Cells()) {
      Vect r(0); // result
      for (size_t e = 0; e < m.GetNumNeighbourFaces(c); ++e) {
        IdxFace f = m.GetNeighbourFace(c, e);
        auto g = gf[f];
        auto n = g / (g.norm() + 1e-6); // TODO: revise 1e-6
        auto s = m.GetOutwardSurface(c, e);
        r += g * s.dot(n);
        r -= s * g.norm();
      }
      r /= m.GetVolume(c);     // div(gg/|g|) - div(|g|I)
      fc_stforce_[c] = r * (-sig);
    }
  }
}

// TODO: Test m.Dump() with pending Comm
template <class M>
void Hydro<M>::Dump(Sem& sem) {
  if (sem("dump")) {
    // requirements: 
    // * interval between dumps is at least dumpdt + dt
    // * dumpt % dtumpdt <= dt * 0.5
    auto pd = par.Double["dumpdt"]; // dum[p] [d]t
    auto& pt = st_.dumpt;
    const auto& dt = st_.dt;

    // time of next dump
    Scal ptn = int(std::max<Scal>(pt + pd, 0.) / pd + 0.5) * pd;
    assert(ptn > pt);

    if (par.Int["output"] && 
        st_.t >= ptn - dt * 0.5 &&
        st_.dumpn < par.Int["dumpmax"]) {
      st_.dumpt = st_.t;

      fc_velux_ = geom::GetComponent(fs_->GetVelocity(), 0);
      m.Dump(&fc_velux_, "vx");
      fc_veluy_ = geom::GetComponent(fs_->GetVelocity(), 1);
      m.Dump(&fc_veluy_, "vy");
      fc_veluz_ = geom::GetComponent(fs_->GetVelocity(), 2);
      m.Dump(&fc_veluz_, "vz");
      fc_p_ = fs_->GetPressure();
      m.Dump(&fc_p_, "p"); 
      fc_vf_ = as_->GetField();
      m.Dump(&fc_vf_, "vf"); 

      if (broot_) {
        std::cerr << "Dump " 
            << "n=" << st_.dumpn 
            << " t=" << st_.dumpt 
            << " tn=" << ptn
            << std::endl;
      }

      ++st_.dumpn;
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
    if (fs_->GetTime() >= par.Double["tmax"] ||
        step_ > par.Int["max_step"]) {
      sem.LoopBreak();
    } else if (broot_) { 
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
  if (sem.Nested("fs-start")) {
    fs_->StartStep();
  }
  if (par.Int["enable_fluid"]) {
    if (sem.Nested("fs-iters")) {
      auto sn = m.GetSem("iter"); // sem nested
      sn.LoopBegin();
      if (sn.Nested("iter")) {
        fs_->MakeIteration();
      }
      if (sn("reduce")) {
        diff_ = fs_->GetConvergenceIndicator();
        m.Reduce(&diff_, "max");
      }
      if (sn("report")) {
        ++st_.iter;
        if (broot_) {
          std::cout << std::scientific << std::setprecision(16)
              << ".....iter=" << fs_->GetIterationCount()
              << ", diff=" << diff_ << std::endl;
          par.Int["iter"] = st_.iter;
        }
      }
      if (sn("convcheck")) {
        assert(fs_->GetConvergenceIndicator() <= diff_);
        auto it = fs_->GetIterationCount();
        if ((diff_ < tol_ && it >= par.Int["min_iter"]) ||
            it >= par.Int["max_iter"]) {
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
  if (par.Int["enable_advection"]) {
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
  }

  Dump(sem);

  if (sem("dumpstat")) {
    if (broot_) {
      ost_->Write();
    }
  }

  sem.LoopEnd();

  if (sem("fluid-timer")) {
    if (broot_ && par.Int["verbose_fluid_timer"]) {
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
  K* Make(Vars& par, const MyBlockInfo& bi) override {
    return new Hydro<M>(par, bi);
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
