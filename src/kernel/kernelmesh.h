// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include "geom/mesh.h"
#include "parse/vars.h"

namespace generic {

template <size_t dim_>
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

  KernelMesh(Vars& var_, const generic::BlockInfoProxy<M::dim>& bi)
      : var(var_), var_mutable(var_), m(CreateMesh<M>(bi)) {
    m.flags.check_nan = var.Int["CHECKNAN"];
    m.flags.edim = std::min<size_t>(var.Int["dim"], M::dim);
    var_mutable.Vect.Set("cell_length", m.GetCellSize());
    var_mutable.Vect.Set("mesh_size", m.GetGlobalSize());
    var_mutable.Vect.Set("domain_length", m.GetGlobalLength());
  }
  virtual ~KernelMesh() = default;
  virtual void Run() = 0;
  M& GetMesh() {
    return m;
  }

 protected:
  const Vars& var; // read-only configuration, shared by local blocks
  Vars& var_mutable; // mutable configuration, shared by local blocks
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
