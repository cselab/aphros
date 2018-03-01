#pragma once

#include <memory>

#include "hydro/vect.hpp"
#include "hydro/mesh3d.hpp"
#include "Kernel.h"
#include "Vars.h"

template <class M>
M CreateMesh(const MyBlockInfo& bi) {
  using MIdx = typename M::MIdx;
  using B = MyBlock;
  using Vect = typename M::Vect;
  using Rect = geom::Rect<Vect>;
  B& b = *(B*)bi.ptrBlock;
  int hl = bi.hl;
  MIdx s(B::sx, B::sy, B::sz); // block size inner

  Scal h = bi.h_gridpoint;
  auto w = bi.index;   // block index
  auto c = bi.origin; 
  Vect d0(c[0], c[1], c[2]); // origin coord
  Vect d1 = d0 + Vect(s) * h;      // end coord
  Rect d(d0, d1);

  MIdx o(w[0] * s[0], w[1] * s[1], w[2] * s[2]); // origin index
  std::cout 
      << "o=" << o 
      << " dom=" << d0 << "," << d1 
      << " h=" << h
      << std::endl;
  
  return geom::InitUniformMesh<M>(d, o, s, hl);
}

// Abstract Kernel aware of Mesh 
template <class M>
class KernelMesh : public Kernel {
 public:
  using Mesh = M;
  using Scal = double;
  using Vect = typename Mesh::Vect;
  using MIdx = typename Mesh::MIdx;
  static constexpr size_t dim = M::dim;

  KernelMesh(Vars& par, const MyBlockInfo& bi);
  void Run() override = 0;
  M& GetMesh() { return m; }

 protected:
  Vars& par;
  MyBlockInfo bi_;
  M m;
};

template <class M>
KernelMesh<M>::KernelMesh(Vars& par, const MyBlockInfo& bi) 
  : par(par), bi_(bi), m(CreateMesh<M>(bi))
{}

template <class _M>
class KernelMeshFactory : public KernelFactory {
 public:
  using M = _M;
  using K = KernelMesh<M>;
  K* Make(Vars&, const MyBlockInfo&) override = 0;
};

