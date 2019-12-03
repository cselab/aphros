#pragma once

#include "geom/mesh.h"
#include "kernel.h"

template <class M>
M CreateMesh(const MyBlockInfo& bi) {
  using MIdx = typename M::MIdx;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  int hl = bi.hl;

  MIdx bs(bi.bs); // block size inner
  Scal h = bi.h_gridpoint;
  MIdx w(bi.index);   // block index
  Vect d0(bi.origin); // origin coord
  Vect d1 = d0 + Vect(bs) * h;      // end coord
  Rect<Vect> d(d0, d1);
  MIdx gs(bi.gs);
  MIdx o = w * bs; // origin index

  MIdx wmax = gs / bs; // maximum block index
  int id = w[0] + w[1] * wmax[0] + w[2] * wmax[0] * wmax[1]; // unique id

  M m = InitUniformMesh<M>(d, o, bs, hl, bi.isroot, bi.islead, gs, id);
  m.SetMaxComm(bi.maxcomm);
  return m;
}

// Abstract Kernel aware of Mesh. Dependency of DistrMesh.
template <class M_>
class KernelMesh : public Kernel {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;
  static constexpr size_t dim = M::dim;

  KernelMesh(Vars& var, const MyBlockInfo& bi)
      : var(var), var_mutable(var), bi_(bi), m(CreateMesh<M>(bi)) {
    m.SetCN(var.Int["CHECKNAN"]); // TODO: revise, avoid optional setters
    m.SetEdim(var.Int["dim"]);
  }
  void Run() override = 0;
  M& GetMesh() { return m; }
  bool IsRoot() { return bi_.isroot; }
  bool IsLead() { return bi_.islead; }

 protected:
  const Vars& var;    // shared among all blocks on each PEs
  Vars& var_mutable; // shared among all blocks on each PEs
  MyBlockInfo bi_;
  M m;
};

// Abstract KernelFactory aware of Mesh. Dependency of DistrMesh.
template <class M_>
class KernelMeshFactory : public KernelFactory {
 public:
  using M = M_;
  using K = KernelMesh<M>;
  K* Make(Vars&, const MyBlockInfo&) const override = 0;
};

