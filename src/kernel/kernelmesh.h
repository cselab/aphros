// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include "geom/mesh.h"
#include "parse/vars.h"

// TODO: remove h_gridpoint from BlockInfoProxy
struct BlockInfoProxy {
  using Idx = std::array<int, 3>;
  Idx index;
  void* ptrBlock;
  double h_gridpoint;
  double origin[3];
  Idx bs;
  int hl; // number of halo cells
  bool isroot; // root block (one among blocks on all PEs)
  bool islead; // lead block (one per each PE)
  Idx gs; // global size
  size_t maxcomm; // maximum number of communication requests
};

template <class M>
M CreateMesh(const BlockInfoProxy& bi) {
  using MIdx = typename M::MIdx;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  int hl = bi.hl;

  const MIdx bs(bi.bs); // block size inner
  const Scal h = bi.h_gridpoint;
  const MIdx w(bi.index); // block index
  const Vect d0(bi.origin); // origin coord
  const Vect d1 = d0 + Vect(bs) * h; // end coord
  const Rect<Vect> d(d0, d1); // domain box
  const MIdx gs(bi.gs); // global size
  const MIdx o = w * bs; // origin index

  const MIdx wmax = gs / bs; // maximum block index

  const int id = M::Flags::GetIdFromBlock(w, wmax);
  M m = InitUniformMesh<M>(d, o, bs, hl, bi.isroot, bi.islead, gs, id);
  m.flags.global_origin = Vect(0);
  m.flags.global_blocks = wmax;
  m.flags.block_length = Vect(bs) * h;
  return m;
}

// Abstract Kernel aware of Mesh. Dependency of DistrMesh.
template <class M_>
class KernelMesh {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;
  static constexpr size_t dim = M::dim;

  KernelMesh(Vars& var_, const BlockInfoProxy& bi)
      : var(var_), var_mutable(var_), bi_(bi), m(CreateMesh<M>(bi)) {
    m.flags.check_nan_ = var.Int["CHECKNAN"];
    m.flags.
    m.SetEdim(var.Int["dim"]);
  }
  virtual ~KernelMesh() = default;
  virtual void Run() = 0;
  M& GetMesh() {
    return m;
  }
  bool IsRoot() {
    return bi_.isroot;
  }
  bool IsLead() {
    return bi_.islead;
  }

 protected:
  const Vars& var; // shared among all blocks on each PEs
  Vars& var_mutable; // shared among all blocks on each PEs
  BlockInfoProxy bi_;
  M m;
};

// Abstract KernelFactory aware of Mesh. Dependency of DistrMesh.
template <class M_>
class KernelMeshFactory {
 public:
  using M = M_;
  using K = KernelMesh<M>;
  virtual K* Make(Vars&, const BlockInfoProxy&) const = 0;
};
