#pragma once

#include <memory>

#include "hydro/vect.hpp"
#include "hydro/mesh3d.hpp"
#include "hydro/linear.hpp"
#include "Kernel.h"
#include "Vars.h"

template <class M>
M CreateMesh(const MyBlockInfo& bi) {
  using MIdx = typename M::MIdx;
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
  std::cout 
      << "o=" << o 
      << " dom=" << d0 << "," << d1 
      << " h=" << h
      << std::endl;
  
  return geom::InitUniformMesh<M>(d, o, bs, hl);
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

/*
template <class M, class Expr>
void GetLs(const geom::FieldCell<Expr>& s, // field of expressions
           const Mesh& m) {
  using LS = typename M::LS;
  using MIdx = typename M::MIdx;
  using IdxCell = geom::IdxCell;
  LS l;
  // Get stencil from first inner cell
  {
    IdxCell c = *m.Cells().begin(); 
    auto& e = s[c];
    for (size_t j = 0; j < e.size(); ++j) {
      MIdx dm = bc.GetMIdx(e[j].idx) - bc.GetMIdx(c);
      l.st.emplace_back(dm);
    }
  }

  int n = m.Cells().size();
  la.resize(n * l.st.size());
  lt.resize(n, 1.);
  lx.resize(n, 0.);

  // fill matrix coeffs
  {
    size_t i = 0;
    for (auto c : m.Cells()) {
      auto& e = fc_system_[c];
      for (size_t j = 0; j < e.size(); ++j) {
        // Check stencil
        if (e[j].idx != bc.GetIdx(bc.GetMIdx(c) + MIdx(l.st[j]))) {
          std::cerr << "***"
              << " MIdx(c)=" << bc.GetMIdx(c)
              << " MIdx(e[j].idx)=" << bc.GetMIdx(e[j].idx)
              << " l.st[j]=" << MIdx(l.st[j]) 
              << std::endl;
          assert(false);
        }
        lsa_[i] = e[j].coeff;
        ++i;
      }
    }
    assert(i == n * l.st.size());
  }

  // fill rhs and zero solution
  {
    size_t i = 0;
    for (auto c : m.Cells()) {
      auto& e = fc_system_[c];
      lsb_[i] = -e.GetConstant();
      lsx_[i] = 0.;
      ++i;
    }
    assert(i == lsb_.size());
  }

  l.a = &lsa_;
  l.b = &lsb_;
  l.x = &lsx_;
  m.Solve(l);
}
*/
