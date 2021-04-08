// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include "geom/mesh.h"
#include "parse/vars.h"

namespace generic {

template <size_t dim_ = 3>
struct BlockInfoProxy {
  static constexpr size_t dim = dim_;
  using Vect = generic::Vect<double, dim>;
  using MIdx = generic::MIdx<dim>;

  MIdx index;
  Vect cellsize;
  MIdx blocksize;
  int halos;
  bool isroot; // root block, one block over all processors
  bool islead; // lead block, one block in each processor
  MIdx globalsize; // global mesh size
};

} // namespace generic

template <class M>
M CreateMesh(const generic::BlockInfoProxy<M::dim>& p) {
  using MIdx = typename M::MIdx;
  using Vect = typename M::Vect;
  const MIdx bs = p.blocksize;
  const MIdx begin = p.index * bs;
  const Vect h(p.cellsize);
  const Rect<Vect> domain(Vect(begin) * h, Vect(begin + bs) * h);

  const MIdx global_blocks = p.globalsize / bs;
  const int id = M::Flags::GetIdFromBlock(p.index, global_blocks);
  M m(begin, bs, domain, p.halos, p.isroot, p.islead, p.globalsize, id);
  m.flags.global_origin = Vect(0);
  m.flags.global_blocks = global_blocks;
  m.flags.block_length = h * Vect(bs);
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

  KernelMesh(Vars& var_, const generic::BlockInfoProxy<dim>& bi)
      : var(var_), var_mutable(var_), bi_(bi), m(CreateMesh<M>(bi)) {
    m.flags.check_nan = var.Int["CHECKNAN"];
    m.flags.edim = var.Int["dim"];
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
  generic::BlockInfoProxy<dim> bi_;
  M m;
};

// Abstract KernelFactory aware of Mesh. Dependency of DistrMesh.
template <class M_>
class KernelMeshFactory {
 public:
  using M = M_;
  using K = KernelMesh<M>;
  virtual K* Make(Vars&, const generic::BlockInfoProxy<M::dim>&) const = 0;
};
