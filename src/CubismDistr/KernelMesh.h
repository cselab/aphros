#pragma once

#include <memory>

#include "hydro/vect.hpp"
#include "hydro/mesh3d.hpp"
#include "Kernel.h"
#include "Vars.h"

template <class M>
M CreateMesh(const MyBlockInfo& bi) {
  using MIdx = typename M::MIdx;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Rect = geom::Rect<Vect>;
  int hl = bi.hl;

  MIdx bs(bi.bs); // block size inner
  Scal h = bi.h_gridpoint;
  MIdx w(bi.index);   // block index
  Vect d0(bi.origin); // origin coord
  Vect d1 = d0 + Vect(bs) * h;      // end coord
  Rect d(d0, d1);

  MIdx o = w * bs; // origin index
  
  return geom::InitUniformMesh<M>(d, o, bs, hl);
}

// Abstract Kernel aware of Mesh. Dependency of DistrMesh.
template <class M>
class KernelMesh : public Kernel {
 public:
  using Mesh = M;
  using Scal = typename Mesh::Scal;
  using Vect = typename Mesh::Vect;
  using MIdx = typename Mesh::MIdx;
  static constexpr size_t dim = M::dim;

  KernelMesh(Vars& par, const MyBlockInfo& bi);
  void Run() override = 0;
  M& GetMesh() { return m; }
  bool IsRoot() { return bi_.isroot; }
  bool IsLead() { return bi_.islead; }

 protected:
  Vars& par; // Shared among all blocks on each PEs
  MyBlockInfo bi_;
  M m;
};

template <class M>
KernelMesh<M>::KernelMesh(Vars& par, const MyBlockInfo& bi) 
  : par(par), bi_(bi), m(CreateMesh<M>(bi))
{}

// Abstract KernelFactory aware of Mesh. Dependency of DistrMesh.
template <class M_>
class KernelMeshFactory : public KernelFactory {
 public:
  using M = M_;
  KernelMesh<M>* Make(Vars&, const MyBlockInfo&) override = 0;
};



