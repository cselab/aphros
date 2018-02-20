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

#include "hydro/vect.hpp"
#include "hydro/mesh3d.hpp"
#include "hydro/solver.hpp"
#include "hydro/advection.hpp"
#include "hydro/conv_diff.hpp"
#include "hydro/fluid.hpp"

#include "ICubism.h"
#include "ILocal.h"
#include "Kernel.h"
#include "Vars.h"

// Basic Rules:
// 
// 1. Program at interface.
// First describe the interface in a separate class, then implement.
// 
// 2. Use templates when:
// - 
// 
// 3. Use inheritance when:
// - 
// 
// 4. Naming conventions:
// - single letter when possible
// - if first letter, put that word in comment: n // name
// - if another letter, put that word with letter in [...]: a // n[a]me
// - up to 4 letters if possible
// - private variables end with underscore: a_;
//   exceptions: m (mesh)
// 
// 5. Dereference pointers to references or values if possible
// But: use pointer arguments to prevent passing an rvalue
// 
// 6. Lower bound for template argument:
//   template <class B /*: A*/>
// means that argument B needs to be a subtype of A
//  
// 7. Comments:
// - capitalized descriptive for a class or function 
//   (e.g. Creates instance)
// - capitalized imperative for expressions in implementation 
//   (e.g. Create instance)
// - non-capitalized for declarations 
//   (e.g. buffer index)
//
// 8. Data from Kernel returned by reference if possible




template <class M>
class Hydro : public Kernel {
 public:
  using Mesh = M;
  using Scal = double;
  using Vect = typename Mesh::Vect;
  using MIdx = typename Mesh::MIdx;
  using Rect = geom::Rect<Vect>;
  using IdxCell = geom::IdxCell;
  using IdxFace = geom::IdxFace;
  using IdxNode = geom::IdxNode;

  template <class T>
  using FieldCell = geom::FieldCell<T>;
  template <class T>
  using FieldFace = geom::FieldFace<T>;
  template <class T>
  using FieldNode = geom::FieldNode<T>;

  Hydro(Vars& par, const MyBlockInfo& bi);
  void Run() override;
  M& GetMesh() { return m; }

 private:
  M CreateMesh(const MyBlockInfo& bi);

  Vars& par;
  std::string name_;
  MyBlockInfo bi_;
  M m;
  //using AS = solver::AdvectionSolverExplicit<M, FieldFace<Scal>>;
  using AS = solver::ConvectionDiffusionScalarImplicit<M>;
  using FS = solver::FluidSimple<M>;
  FieldCell<Scal> fc_sc_; // scaling
  FieldFace<Scal> ff_d_; // diffusion rate
  FieldCell<Scal> fc_d_; // diffusion rate
  FieldCell<Scal> fc_src_; // source
  FieldFace<Scal> ff_flux_;  // volume flux
  FieldCell<Vect> fc_force_;  // force
  FieldCell<Vect> fc_stforce_;  // stforce cells TODO: what is st
  FieldFace<Vect> ff_stforce_;  // stforce faces
  MultiTimer<std::string> timer_; 
  std::shared_ptr<const solver::LinearSolverFactory> p_lsf_; // linear solver factory
  std::unique_ptr<AS> as_; // advection solver
  std::unique_ptr<FS> fs_; // fluid solver
  FieldCell<Scal> fc_velu_; // single velocity component
  FieldCell<Scal> fc_p_; // pressure
  Scal sum_;
};

template <class M /*: Mesh*/>
void Grad(M& m) {
  auto sem = m.GetSem("grad");
  if (sem()) {
    std::cerr << sem.GetName() <<  ":s1" << std::endl;
  }
  if (sem()) {
    std::cerr << sem.GetName() <<  ":s2" << std::endl;
  }
}


template <class M>
M Hydro<M>::CreateMesh(const MyBlockInfo& bi) {
  using B = MyBlock;
  B& b = *(B*)bi.ptrBlock;
  int hl = bi.hl;
  MIdx si(B::sx, B::sy, B::sz); // block size inner
  MIdx s = si + MIdx(hl * 2);   // block size with halos

  Scal h = bi.h_gridpoint;
  auto w = bi.index;   // block index
  auto c = bi.origin; 
  Vect d0(c[0], c[1], c[2]); // origin coord
  Vect d1 = d0 + Vect(si) * h;      // end coord
  Rect d(d0, d1);

  MIdx o(w[0] * si[0], w[1] * si[1], w[2] * si[2]); // origin index
  o -= MIdx(hl);
  std::cout 
      << "o=" << o 
      << " dom=" << d0 << "," << d1 
      << " h=" << h
      << std::endl;
  
  return geom::InitUniformMesh<M>(d, o, s, hl);
}

template <class M>
Hydro<M>::Hydro(Vars& par, const MyBlockInfo& bi) 
  : par(par), bi_(bi), m(CreateMesh(bi))
{
  name_ = 
      "[" + std::to_string(bi.index[0]) +
      "," + std::to_string(bi.index[1]) +
      "," + std::to_string(bi.index[2]) + "]";

  // initial field for advection
  FieldCell<Scal> fc_u(m);
  for (auto i : m.Cells()) {
    const Scal kx = 2. * M_PI;
    const Scal ky = 2. * M_PI;
    const Scal kz = 2. * M_PI;
    Vect c = m.GetCenter(i);
    fc_u[i] = std::sin(kx * c[0]) * std::sin(ky * c[1]) * std::sin(kz * c[2]);
    //fc_u[i] = fc_u[i] > 0. ? 1. : -1.;
  }

  // zero-derivative boundary conditions for advection
  geom::MapFace<std::shared_ptr<solver::ConditionFace>> mf_cond;
  for (auto idxface : m.Faces()) {
    if (!m.IsExcluded(idxface) && !m.IsInner(idxface)) {
      mf_cond[idxface] =
          std::make_shared
          <solver::ConditionFaceDerivativeFixed<Scal>>(Scal(0));
    }
  }

  const Vect vel(par.Vect["vel"]);

  // initial velocity
  FieldCell<Vect> fc_vel(m, Vect(0));

  // zero-derivative boundary conditions for velocity
  geom::MapFace<std::shared_ptr<solver::ConditionFaceFluid>> mf_velcond;
  for (auto idxface : m.Faces()) {
    if (!m.IsExcluded(idxface) && !m.IsInner(idxface)) {
      mf_velcond[idxface] =
          std::make_shared
          <solver::fluid_condition::NoSlipWallFixed<M>>(vel);
    }
  }
  
  // empty cell conditions for advection
  geom::MapCell<std::shared_ptr<solver::ConditionCell>> mc_cond;


  // empty cell conditions for velocity
  geom::MapCell<std::shared_ptr<solver::ConditionCellFluid>> mc_velcond;
  {
    Vect x(par.Vect["pfixed_x"]);
    IdxCell c = m.FindNearestCell(x);
    Scal p = par.Double["pfixed"];
    if (m.GetCenter(c).dist(x) < 0.1) { // TODO: choose nearest block
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

  Scal relax = par.Double["relax"];
  Scal prelax = par.Double["prelax"];
  Scal vrelax = par.Double["vrelax"];
  Scal rhie = par.Double["rhie"];
  Scal tol = par.Double["tol"];
  int num_iter = par.Int["num_iter"];
  bool so = par.Int["second_order"];

  fc_sc_.Reinit(m, 1.); // scaling for as_ and density for fs_
  ff_d_.Reinit(m, par.Double["mu"]);
  fc_d_.Reinit(m, par.Double["mu"]);
  fc_stforce_.Reinit(m, Vect(0));
  ff_stforce_.Reinit(m, Vect(0));
  fc_src_.Reinit(m, 0.);

  fc_force_.Reinit(m);
  const Vect f(par.Vect["force"]);
  const Vect f2(par.Vect["force2"]);
  for (auto i : m.Cells()) {
    Vect x = m.GetCenter(i);
    if ((x[1] - 0.5) * (x[2] - 0.25) < 0.) {
      fc_force_[i] = f;
    } else {
      fc_force_[i] = f2;
    }
  }

  p_lsf_ = std::make_shared<const solver::LinearSolverFactory>(
        std::make_shared<const solver::LuDecompositionFactory>());

  // Init advection solver
  //as_.reset(new AS(m, fc_u, mf_cond, &ff_flux_, &fc_src_, 0., dt));
  /*
  as_.reset(new AS(
        m, fc_u, mf_cond, mc_cond, 
        relax, 
        &fc_sc_, &ff_d_, &fc_src_, &ff_flux_,
        0., dt,
        *p_lsf_, 
        tol, num_iter, 
        so
        ));
        */

  fs_.reset(new FS(
        m, fc_vel, 
        mf_velcond, mc_velcond, 
        vrelax, prelax, rhie,
        &fc_sc_, &fc_d_, &fc_force_, &fc_stforce_, &ff_stforce_, 
        &fc_src_, &fc_src_,
        0., dt,
        *p_lsf_, *p_lsf_,
        tol, num_iter, 
        &timer_, 
        so, false, false, 0., Vect(0)
        ));
}


template <class M>
void Hydro<M>::Run() {
  auto sem = m.GetSem("run");

  /*
  if (sem.Nested("as->StartStep()")) {
    as_->StartStep();
  }
  if (sem.Nested("as->MakeIteration")) {
    as_->MakeIteration();
  }
  if (sem.Nested("as->FinishStep()")) {
    as_->FinishStep();
  }
  */

  if (sem.Nested("fs->StartStep()")) {
    fs_->StartStep();
  }
  if (sem.Nested("fs->MakeIteration")) {
    fs_->MakeIteration();
  }
  if (sem.Nested("fs->FinishStep()")) {
    fs_->FinishStep();
  }
  if (sem("Comm(&fc_velu_)")) {
    // advection
    //auto& u = const_cast<FieldCell<Scal>&>(as_->GetField());
    //m.Comm(&u);
    // fluid velocity single component
    //fc_velu_ = geom::GetComponent(fs_->GetVelocity(), 0);
    //m.Comm(&fc_velu_); // goes to dumper
    fc_p_ = fs_->GetPressure();
    m.Comm(&fc_p_); // goes to dumper
  }
}

template <class MM>
class HydroFactory : public KernelFactory {
 public:
  using M = MM;
  using K = Hydro<M>;
  std::unique_ptr<Kernel> Make(Vars& par, const MyBlockInfo& bi) override {
    return std::unique_ptr<Hydro<M>>(new Hydro<M>(par, bi));
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
