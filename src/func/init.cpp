// Created by Petr Karnakov on 28.12.2019
// Copyright 2019 ETH Zurich

#include "init.h"
#include "geom/mesh.h"
#include "init_cl.h"
#include "init_sig.h"
#include "init_u.h"
#include "overlap/overlap.h"

template <class M>
void InitVfList(
    FieldCell<typename M::Scal>& fc, std::istream& list, int approx,
    size_t edim, const M& m, bool verbose) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Primitive = typename UPrimList<Vect>::Primitive;
  const std::vector<Primitive> ppa = UPrimList<Vect>::GetPrimitives(list, edim);
  if (verbose && m.IsRoot()) {
    std::cerr << "Read " << ppa.size() << " primitives" << std::endl;
  }

  const Vect h = m.GetCellSize();
  // filter to bounding box
  auto& bc = m.GetSuBlockCells();
  Rect<Vect> rect(Vect(bc.GetBegin()) * h, Vect(bc.GetEnd()) * h);

  std::vector<Primitive> pp;
  for (auto& p : ppa) {
    if (p.inter(rect)) {
      pp.push_back(p);
    }
  }

  if (pp.empty()) {
    fc.Reinit(m, 0.);
  } else {
    auto lsmax = [&pp](Vect x) -> std::pair<Scal, size_t> {
      Scal lmax = -std::numeric_limits<Scal>::max(); // maximum level-set
      size_t imax = 0; // index of maximum
      for (size_t i = 0; i < pp.size(); ++i) {
        auto& p = pp[i];
        Scal li = p.ls(x);
        if (p.mod_minus) {
          li = -li;
        }
        if (p.mod_and) {
          lmax = std::min(lmax, li);
        } else {
          if (li > lmax) {
            lmax = li;
            imax = i;
          }
        }
      }
      return {lmax, imax};
    };
    if (approx == 0) { // stepwise
      for (auto c : m.Cells()) {
        fc[c] = (lsmax(m.GetCenter(c)).first >= 0 ? 1 : 0);
      }
    } else if (approx == 1) { // level set
      for (auto c : m.Cells()) {
        const auto x = m.GetCenter(c);
        const auto& p = pp[lsmax(x).second];
        fc[c] = GetLevelSetVolume<Scal>(p.ls, x, h);
      }
    } else if (approx == 2) { // overlap
#if USEFLAG(OVERLAP)
      for (auto c : m.Cells()) {
        const auto x = m.GetCenter(c);
        const auto& p = pp[lsmax(m.GetCenter(c)).second];
        Vect qx = (x - p.c) / p.r;
        Vect qh = h / p.r;
        if (edim == 2 && M::dim == 3) {
          qx[2] = 0.;
          qh[2] *= 1e-3; // XXX: adhoc, thin cell for 2d
        }
        using Vect3 = generic::Vect<Scal, 3>;
        fc[c] = GetSphereOverlap(Vect3(qx), Vect3(qh), Vect3(0), 1);
      }
#else
      fassert(false, "overlap is disabled");
#endif
    } else if (approx == 3) { // level-set on nodes
      FieldNode<Scal> fnl(m);
      for (auto n : m.Nodes()) {
        fnl[n] = lsmax(m.GetNode(n)).first;
      }
      fc = GetPositiveVolumeFraction(fnl, m);
    } else {
      fassert(false, "unknown approx=" + std::to_string(approx));
    }
  }
}

#define XX(M)                                                                 \
  template void InitVf(                                                       \
      FieldCell<typename M::Scal>& fcu, const Vars& var, M& m, bool verbose); \
  template std::function<void(                                                \
      FieldCell<typename M::Scal>&, const FieldCell<typename M::Scal>&,       \
      const M&)>                                                              \
  CreateInitCl(const Vars& par, bool verb);                                   \
  template std::function<void(FieldCell<typename M::Scal>&, const M&)>        \
  CreateInitU(const Vars& par, bool verb);                                    \
  template std::function<void(FieldCell<typename M::Scal>&, const M&)>        \
  CreateInitSig(const Vars& var);

#define COMMA ,
#define X(dim) XX(MeshCartesian<double COMMA dim>)
MULTIDIMX
