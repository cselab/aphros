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
  using Scal = double;
  using Vect = typename Mesh::Vect;
  using MIdx = typename Mesh::MIdx;
  using Rect = geom::Rect<Vect>;
  using IdxCell = geom::IdxCell;
  using IdxFace = geom::IdxFace;
  using IdxNode = geom::IdxNode;
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
  void CalcMixture(const FieldCell<Scal>& vf);
  void CalcStat();

  using AS = solver::AdvectionSolverExplicit<M, FieldFace<Scal>>;
  using FS = solver::FluidSimple<M>;
  FieldCell<Scal> fc_mu_; // viscosity
  FieldCell<Scal> fc_rho_; // density
  FieldCell<Scal> fc_src_; // source
  FieldFace<Scal> ff_flux_;  // volume flux
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

  struct Stat {
    Scal m1, m2; // mass
    Vect c1, c2;  // center of mass 
    Vect vc1, vc2;  // center of mass velocity
    Vect v1, v2;  // average velocity
    Stat()
      : m1(0), m2(0), c1(0), c2(0), vc1(0), vc2(0), v1(0), v2(0)
    {}
    void Clear() {
      (*this) = Stat();
    }
  };
  Stat st_;
  std::shared_ptr<output::Session> ost_; // output stat
};


// TODO: move construction to Run()
template <class M>
Hydro<M>::Hydro(Vars& par, const MyBlockInfo& bi) 
  : KernelMesh<M>(par, bi)
{
  broot_ = (m.GetInBlockCells().GetBegin() == MIdx(0));

  par.Int.Set("iter", 0);

  // initial field for advection
  FieldCell<Scal> fc_u(m);
  const std::string vi = par.String["vf_init"];
  if (vi == "sin") {
    Vect k;
    if (auto p = par.Vect("sin_k")) {
      k = Vect(*p);
    } else {
      k = Vect(2. * M_PI);
    }

    for (auto i : m.Cells()) {
      Vect z = m.GetCenter(i) * k;
      fc_u[i] = std::sin(z[0]) * std::sin(z[1]) * std::sin(z[2]);
    }
  } else if (vi == "circle") {
    const Vect c(par.Vect["circle_c"]);
    const Scal r(par.Double["circle_r"]);
    for (auto i : m.Cells()) {
      Vect x = m.GetCenter(i);
      fc_u[i] = (c.dist(x) < r ? 1. : 0.);
    }
  } else {
    std::cerr << "Unknown vf_init=" << vi << std::endl;
    assert(false);
    // TODO: add assert for release mode (NDEBUG=1)
  }

  const Vect vel(par.Vect["vel"]);

  // initial velocity
  FieldCell<Vect> fc_vel(m, Vect(0));
  if (par.Int["taylor-green"]) {
    for (auto i : m.AllCells()) {
      auto& u = fc_vel[i][0];
      auto& v = fc_vel[i][1];
      Scal l = (2 * M_PI);
      Scal x = m.GetCenter(i)[0] * l;
      Scal y = m.GetCenter(i)[1] * l;
      u = std::cos(x) * std::sin(y);
      v = -std::sin(x) * std::cos(y);
    }
  } 

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
  auto gxm = [this](IdxFace i) {
    return m.GetDir(i) == Dir::i &&
        m.GetBlockFaces().GetMIdx(i)[0] == 0;
  };
  auto gxp = [this,gs](IdxFace i) {
    return m.GetDir(i) == Dir::i &&
        m.GetBlockFaces().GetMIdx(i)[0] == gs[0];
  };
  auto gym = [this](IdxFace i) {
    return m.GetDir(i) == Dir::j &&
        m.GetBlockFaces().GetMIdx(i)[1] == 0;
  };
  auto gyp = [this,gs](IdxFace i) {
    return m.GetDir(i) == Dir::j &&
        m.GetBlockFaces().GetMIdx(i)[1] == gs[1];
  };
  auto gzm = [this](IdxFace i) {
    return dim >= 3 && m.GetDir(i) == Dir::k &&
        m.GetBlockFaces().GetMIdx(i)[2] == 0;
  };
  auto gzp = [this,gs](IdxFace i) {
    return dim >= 3 && m.GetDir(i) == Dir::k &&
        m.GetBlockFaces().GetMIdx(i)[2] == gs[2];
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
    // TODO: Reduce and choose nearest block
    if (m.GetCenter(c).dist(x) < 0.1) { // XXX: adhoc
      mc_velcond[c] = std::make_shared
          <solver::fluid_condition::GivenPressureFixed<Mesh>>(p);
    }

  }

  // velocity and flux
  ff_flux_.Reinit(m);
  for (auto idxface : m.Faces()) {
    ff_flux_[idxface] = vel.dot(m.GetSurface(idxface));
  }

  // time step
  const Scal dt = par.Double["dt"];
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

  // Output from par.Double
  auto od = [this, &par](std::string n /*output-name*/,  
                        std::string p /*parameter*/) {
      return std::make_shared<output::EntryScalarFunction<Scal>>(
          n, [&par, p](){ return par.Double[p]; });
    };

  // Output by pointer
  auto op = [this](std::string n /*output-name*/,  Scal* p /*pointer*/) {
      return std::make_shared<output::EntryScalarFunction<Scal>>(
          n, [p](){ return *p; });
    };


  par.Double.Set("t", fs_->GetTime());

  if (broot_) {
    auto& s = st_;
    output::Content con = {
        od("t", "t"),
        std::make_shared<output::EntryScalarFunction<Scal>>(
            "iter", [&par](){ return par.Int["iter"]; }),
        op("diff", &diff_),
        op("m1", &s.m1),
        op("m2", &s.m2),
        op("c1x", &s.c1[0]), op("c1y", &s.c1[1]), op("c1z", &s.c1[2]),
        op("c2x", &s.c2[0]), op("c2y", &s.c2[1]), op("c2z", &s.c2[2]),
        op("vc1x", &s.vc1[0]), op("vc1y", &s.vc1[1]), op("vc1z", &s.vc1[2]),
        op("vc2x", &s.vc2[0]), op("vc2y", &s.vc2[1]), op("vc2z", &s.vc2[2]),
        op("v1x", &s.v1[0]), op("v1y", &s.v1[1]), op("v1z", &s.v1[2]),
        op("v2x", &s.v2[0]), op("v2y", &s.v2[1]), op("v2z", &s.v2[2]),
    };
    ost_ = std::make_shared<
        output::SessionPlainScalar<Scal>>(con, "stat.dat");
  }
}

template <class M>
void Hydro<M>::CalcStat() {
  auto sem = m.GetSem("stat");

  auto& s = st_;

  if (sem("local")) {
    par.Double.Set("t", fs_->GetTime());

    auto& fa = as_->GetField();
    auto& fv = fs_->GetVelocity();

    // Save c for vc
    Vect c1p = st_.c1;
    Vect c2p = st_.c2;

    st_.Clear();

    // TODO: revise
    // Restore c to vc
    st_.vc1 = c1p;
    st_.vc2 = c2p;

    // mass, center, velocity
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
      if (broot_) {
        std::cout << "meshvel = " << v << std::endl;
      }
      double w = par.Double["meshvel_weight"];
      Vect vp = fs_->GetMeshVel();
      fs_->SetMeshVel(v * w + vp * (1. - w));
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


template <class M>
void Hydro<M>::Run() {
  auto sem = m.GetSem("run");

  sem.LoopBegin();

  if (sem("loop-check")) {
    ++step_;
    if (step_ > par.Int["max_step"]) {
      sem.LoopBreak();
    } else if (broot_) { 
      std::cerr 
          << "STEP=" << step_ 
          << " t=" << fs_->GetTime() << std::endl;
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
  if (sem.Nested("as-start")) {
    as_->StartStep();
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
        if (broot_) {
          std::cout << std::scientific << std::setprecision(16)
              << ".....iter=" << fs_->GetIterationCount()
              << ", diff=" << diff_ << std::endl;
          ++par.Int["iter"];
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
  if (par.Int["enable_advection"]) {
    if (sem.Nested("as-iter")) {
      as_->MakeIteration();
    }
  }
  if (sem.Nested("fs-finish")) {
    fs_->FinishStep();
  }
  if (sem.Nested("as-finish")) {
    as_->FinishStep();
  }
  // TODO: Test Dump with pending Comm
  if (sem("dump")) {
    if (par.Int["output"] && 
        step_ % (par.Int["max_step"] / par.Int["num_frames"])  == 0) {
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
    }
  }
  if (sem("dumpstat")) {
    if (broot_) {
      ost_->Write();
    }
  }
  if (sem("dumpwrite")) {
    // Empty stage for DumpWrite
    // TODO: revise
  }

  sem.LoopEnd();

  if (sem("fluid-timer")) {
    if (broot_) {
      double a = 0.; // total
      auto& mt = timer_;
      for (auto e : mt.GetMap()) {
        a += e.second;
      }

      std::cout << std::fixed;
      for (auto e : mt.GetMap()) {
        auto n = e.first; // name
        if (n == "") {
          n = "other";
        }
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
