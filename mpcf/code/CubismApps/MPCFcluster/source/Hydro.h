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

#include "../../hydro/vect.hpp"
#include "../../hydro/mesh3d.hpp"
#include "../../hydro/solver.hpp"
#include "../../hydro/advection.hpp"

#include "ICubism.h"
#include "ILocal.h"
#include "Kernel.h"

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

  Hydro(const MyBlockInfo& bi);
  void Run() override;
  //void ReadBuffer(LabMPI& l) override;
  //void WriteBuffer(Block_t& o) override;
  M& GetMesh() { return m; }

 private:
  M CreateMesh(const MyBlockInfo& bi);
  using LS = typename Mesh::LS;

  std::string name_;
  MyBlockInfo bi_;
  M m;
  using AS = solver::AdvectionSolverExplicit<M, FieldFace<Scal>>;
  FieldCell<Scal> fc_src_;
  FieldFace<Scal> ff_flux_;
  FieldCell<Scal> fc_p_;
  std::unique_ptr<AS> as_;
  Scal sum_;
  // LS
  std::vector<Scal> lsa_, lsb_, lsx_;
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
  int hl = 1;
  // TODO: hl from Distr
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
Hydro<M>::Hydro(const MyBlockInfo& bi) 
  : bi_(bi), m(CreateMesh(bi))
  , fc_src_(m, 0.), ff_flux_(m), fc_p_(m)
{
  name_ = 
      "[" + std::to_string(bi.index[0]) +
      "," + std::to_string(bi.index[1]) +
      "," + std::to_string(bi.index[2]) + "]";

  // Initial field for advection
  FieldCell<Scal> fc_u(m);
  for (auto i : m.Cells()) {
    const Scal kx = 2. * M_PI;
    const Scal ky = 2. * M_PI;
    const Scal kz = 2. * M_PI;
    Vect c = m.GetCenter(i);
    fc_u[i] = std::sin(kx * c[0]) * std::sin(ky * c[1]) * std::sin(kz * c[2]);
    //fc_u[i] = fc_u[i] > 0. ? 1. : -1.;
  }

  // zero-derivative boundary conditions
  geom::MapFace<std::shared_ptr<solver::ConditionFace>> mf_cond;
  for (auto idxface : m.Faces()) {
    if (!m.IsExcluded(idxface) && !m.IsInner(idxface)) {
      mf_cond[idxface] =
          std::make_shared<solver::ConditionFaceDerivativeFixed<Scal>>(Scal(0));
    }
  }

  // velocity and flux
  const Vect vel(1, 1, 1);
  for (auto idxface : m.Faces()) {
    ff_flux_[idxface] = vel.dot(m.GetSurface(idxface));
  }

  // time step
  const Scal dt = 0.0025;

  // Init advection solver
  as_.reset(new AS(m, fc_u, mf_cond, &ff_flux_, &fc_src_, 0., dt));
}

template <class M>
void Hydro<M>::Run() {
  auto sem = m.GetSem();

  if (sem()) {
    as_->StartStep();
    as_->MakeIteration();
    as_->FinishStep();
    m.Comm(&const_cast<FieldCell<Scal>&>(as_->GetField()));

    sum_ = 0.;
    for (auto i : m.Cells()) {
      sum_ += as_->GetField()[i];
    }
    m.Reduce(&sum_);

    // linear system
    // Each block computes the coefficients assuming a uniform stencil
    // (requirement of hypre)
    // Then it computes the rhs and allocates space for result.
    // All three are 1D arrays.
    // Then the block issues a request to solve a linear system 
    // passing pointers to these arrays and stencil description.
    // After going through all blocks, 
    // the processor assembles the system and calls hypre.
    LS l;
    l.st.emplace_back(0, 0, 0);
    l.st.emplace_back(-1, 0, 0);
    l.st.emplace_back(1, 0, 0);
    l.st.emplace_back(0, -1, 0);
    l.st.emplace_back(0, 1, 0);
    l.st.emplace_back(0, 0, -1);
    l.st.emplace_back(0, 0, 1);
    int bs = _BLOCKSIZE_;
    int n = bs * bs *bs;
    lsa_.resize(n*l.st.size());
    for (int i = 0; i < lsa_.size();) {
      lsa_[i++] = -6.;
      lsa_[i++] = 1.;
      lsa_[i++] = 1.;
      lsa_[i++] = 1.;
      lsa_[i++] = 1.;
      lsa_[i++] = 1.;
      lsa_[i++] = 1.;
    }

    lsb_.resize(n, 1.);
    size_t j = 0;
    auto& bc = m.GetBlockCells();
    auto& u = const_cast<FieldCell<Scal>&>(as_->GetField());
    for (auto i : m.Cells()) {
      auto d = m.GetBlockCells().GetMIdx(i) - MIdx(1) - bc.GetBegin(); // TODO: 1 -> h
      if (MIdx(0) <= d && d < MIdx(bs)) {
        lsb_[j++] = u[i];
      }
    }
    assert(j == lsb_.size());

    lsx_.resize(n, 0.);
    l.a = &lsa_;
    l.b = &lsb_;
    l.x = &lsx_;
    m.Solve(l);
  }
  if (sem()) {
    int bs = _BLOCKSIZE_;
    size_t j = 0;
    auto& bc = m.GetBlockCells();
    for (auto i : m.Cells()) {
      auto d = m.GetBlockCells().GetMIdx(i) - MIdx(1) - bc.GetBegin(); // TODO: 1 -> h
      if (MIdx(0) <= d && d < MIdx(bs)) {
        fc_p_[i] = lsx_[j++];
      }
    }
    assert(j == lsx_.size());
    for (auto i : m.Cells()) {
      --j;
    }
    assert(j == 0);
    m.Comm(&fc_p_);
  }
}

/*
template <class M>
void Hydro<M>::ReadBuffer(LabMPI& l) {
  int e = 0; // buffer field idx

  for (auto u : vcm_) {
    for (auto i : m.Cells()) {
      auto& bc = m.GetBlockCells();
      auto d = bc.GetMIdx(i) - MIdx(1) - bc.GetBegin(); // TODO: 1 -> h
      (*u)[i] = l(d[0], d[1], d[2]).a[e];
    }
    ++e;
  }

  vcm_.clear();
}

template <class M>
void Hydro<M>::WriteBuffer(Block_t& o) {
  int bs = _BLOCKSIZE_;

  // Check buffer has enough space for all fields
  assert(vcm_.size() <= Elem::s);

  int e = 0; // buffer field idx

  for (auto u : vcm_) {
    for (auto i : m.Cells()) {
      auto& bc = m.GetBlockCells();
      auto d = m.GetBlockCells().GetMIdx(i) - MIdx(1) - bc.GetBegin(); // TODO: 1 -> h
      if (MIdx(0) <= d && d < MIdx(bs)) {
        o.data[d[2]][d[1]][d[0]].a[e] = (*u)[i];
      }
    }
    ++e;
  }
}
*/

// Class with field 'stencil' needed for SynchronizerMPI::sync(Processing)
template <class MM>
class HydroFactory : public KernelFactory {
 public:
  using M = MM;
  using K = Hydro<M>;
  std::unique_ptr<Kernel> Make(const MyBlockInfo& bi) override {
    return std::unique_ptr<Hydro<M>>(new Hydro<M>(bi));
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
