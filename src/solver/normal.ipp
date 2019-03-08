#include <memory>
#include <cmath>

#include "geom/mesh.h"
#include "solver.h"

// normal from exact sphere
#define ADHOC_NORM 0

namespace solver {

template <class M_>
struct UNormal<M_>::Imp {
  static constexpr size_t dim = M::dim;

  static auto Maxmod(Scal a, Scal b) -> Scal {
    return std::abs(b) < std::abs(a) ? a : b;
  }

  // Computes normal by gradient.
  // uc: volume fraction
  // mfc: boundary conditions for volume fraction
  // Output: modified in cells with msk=1
  // fcn: normal [s] 
  static void CalcNormalGrad(M& m, const FieldCell<Scal>& uc,
                             const MapFace<std::shared_ptr<CondFace>>& mfc,
                              FieldCell<Vect>& fcn) {
    auto uf = Interpolate(uc, mfc, m);
    auto gc = Gradient(uf, m);
    for (auto c : m.AllCells()) {
      Vect g = gc[c];
      fcn[c] = g;
    }
  }
  // Computes normal by Young's scheme (interpolation of gradient from nodes).
  // fcu: volume fraction
  // fci: interface mask (1: contains interface)
  // Output: modified in cells with fci=1, resized to m
  // fcn: normal with norm1()=1, antigradient of fcu [s] 
  // XXX: uses static variables, not suspendable
  // TODO: check non-uniform mesh
  static void CalcNormalYoung(M& m, const FieldCell<Scal>& fcu, 
                              const FieldCell<bool>& fci,
                              FieldCell<Vect>& fcn) {
    static FieldNode<Vect> g;  // gradient
    static FieldNode<Vect> l;  // step from cell to node
    g.Reinit(m, Vect(0));
    l.Reinit(m, Vect(0));
    // values from cells to neighbour nodes
    for (auto c : m.AllCells()) {
      Vect xc = m.GetCenter(c);
      for (size_t q = 0; q < m.GetNumNeighbourNodes(c); ++q) {
        IdxNode n = m.GetNeighbourNode(c, q);
        Vect xn = m.GetNode(n);
        for (size_t d = 0; d < dim; ++d) {
          g[n][d] += (xc[d] - xn[d] > 0. ? 1. : -1.) * fcu[c];
          l[n][d] += std::abs(xc[d] - xn[d]);
        }
      } 
    } 
    // gradient on nodes
    for (auto n : m.SuNodes()) {
      g[n] /= l[n];
    }

    // gradient on cells
    fcn.Reinit(m);
    for (auto c : m.SuCells()) {
      if (fci[c]) {
        // sum over neighbour nodes
        auto& v = fcn[c];
        v = Vect(0);
        for (size_t q = 0; q < m.GetNumNeighbourNodes(c); ++q) {
          IdxNode n = m.GetNeighbourNode(c, q);
          v += g[n];
        }
        // normalize
        v /= -v.norm1();
      }
    }
  }
  // Computes normal and curvature from height functions.
  // fcu: volume fraction
  // fci: interface mask (1: contains interface)
  // edim: effective dimension
  // ow: 1: force overwrite, 0: update only if gives steeper profile
  // Output: modified in cells with fci=1, resized to m
  // fcn: normal, antigradient of fcu, if gives steeper profile or ow=1 [s] 
  // fck: curvature [s] 
  static void CalcNormalHeight(
      M& m, const FieldCell<Scal>& fcu, const FieldCell<bool>& fci,
      size_t edim, bool ow, FieldCell<Vect>& fcn, FieldCell<Scal>& fck) {
    using MIdx = typename M::MIdx;
    using Dir = typename M::Dir;
    auto& bc = m.GetIndexCells();

    fcn.Reinit(m); // XXX
    fck.Reinit(m);

    for (auto c : m.SuCells()) {
      if (!fci[c]) {
        continue;
      }
      Vect bn; // best normal
      Scal bhx = 0., bhy = 0.; // best first derivative
      Scal bk = 0.; // best curvature[k]
      Dir bd = Dir::i;  // best direction
      bool fst = true; // first
      std::vector<Dir> dd; // direction of plane normal
      if (edim == 2) {
        dd = {Dir::i, Dir::j};
      } else {
        dd = {Dir::i, Dir::j, Dir::k};
      }
      for (Dir dn : dd) {
        // directions of plane tangents ([d]irection [t]angents)
        Dir dtx((size_t(dn) + 1) % dim); 
        Dir dty((size_t(dn) + 2) % dim); 

        MIdx w = bc.GetMIdx(c);

        // offset in normal direction
        MIdx on = MIdx(dn);
        // offset in dtx,dty
        MIdx otx = MIdx(dtx);
        MIdx oty = MIdx(dty);
        // mesh step
        const Scal lx = m.GetCenter(c).dist(m.GetCenter(bc.GetIdx(w - otx)));
        const Scal ly = m.GetCenter(c).dist(m.GetCenter(bc.GetIdx(w - oty)));
        const Scal ln = m.GetCenter(c).dist(m.GetCenter(bc.GetIdx(w - on)));

        // Evaluates height function
        // o: offset from w
        auto hh = [&](MIdx o) -> Scal {
          return 
            (fcu[bc.GetIdx(w + o - on)] + 
            fcu[bc.GetIdx(w + o)] + 
            fcu[bc.GetIdx(w + o + on)]) * ln;
        };

        // height function
        const Scal h = hh(MIdx(0));
        const Scal hxm = hh(-otx);
        const Scal hxp = hh(otx);
        const Scal hym = hh(-oty);
        const Scal hyp = hh(oty);
        // corners: hxy
        const Scal hmm = hh(-otx - oty); 
        const Scal hmp = hh(-otx + oty);
        const Scal hpm = hh(otx - oty);
        const Scal hpp = hh(otx + oty);

        // first derivative (slope)
        Scal hx = (hxp - hxm) / (2. * lx);  // centered
        Scal hy = (hyp - hym) / (2. * ly); 
        // sign: +1 if u increases in dn
        Scal sg = 
            (fcu[bc.GetIdx(w + on)] - fcu[bc.GetIdx(w - on)] > 0. ? 1. : -1.);
        // second derivative 
        Scal hxx = (hxp - 2. * h + hxm) / (lx * lx);
        Scal hyy = (hyp - 2. * h + hym) / (ly * ly);
        Scal hxy = ((hpp - hmp) - (hpm - hmm)) / (4. * lx * ly);
        // curvature
        Scal k = (2. * hx * hy * hxy 
            -(sqr(hy) + 1.) * hxx -(sqr(hx) + 1.) * hyy) / 
            std::pow(sqr(hx) + sqr(hy) + 1., 3. / 2.);
        // outer normal
        Vect n;
        n[size_t(dtx)] = -hx;
        n[size_t(dty)] = -hy;
        n[size_t(dn)] = -sg;
        // select best with minimal slope
        if (fst || 
            std::abs(hx) + std::abs(hy) < std::abs(bhx) + std::abs(bhy)) {
          bn = n;
          bhx = hx;
          bhy = hy;
          bk = k;
          bd = dn;
          fst = false;
        } 
      }
      bn /= bn.norm1(); // normalize

      // update if ow=1 or gives steeper profile in plane dn
      if (ow || std::abs(bn[size_t(bd)]) < std::abs(fcn[c][size_t(bd)])) {
        fcn[c] = bn;
      }

      // curvature
      fck[c] = bk;

      #if ADHOC_NORM
      auto cr = GetBubble();
      auto q = m.GetCenter(c) - cr.first;
      if (edim == 2) {
        q[2] = 0.;
      }
      fcn[c] = q / q.norm();
      #endif 
    }
  }
  // Computes normal by combined Young's scheme and height-functions
  // fcu: volume fraction
  // fci: interface mask (1: contains interface)
  // edim: effective dimension
  // Output: set to NaN if fci=0
  // fcn: normal with norm1()=1, antigradient of fcu [s] 
  // fck: curvature
  static void CalcNormal(
      M& m, const FieldCell<Scal>& fcu, const FieldCell<bool>& fci,
      size_t edim, FieldCell<Vect>& fcn, FieldCell<Scal>& fck) {
    fcn.Reinit(m, Vect(GetNan<Scal>()));
    fck.Reinit(m, GetNan<Scal>());
    CalcNormalYoung(m,fcu, fci, fcn);
    CalcNormalHeight(m, fcu, fci, edim, false, fcn, fck);
  }
};

template <class M_>
void UNormal<M_>::CalcNormal(
    M& m, const FieldCell<Scal>& fcu, const FieldCell<bool>& fci,
    size_t edim, FieldCell<Vect>& fcn, FieldCell<Scal>& fck) {
  Imp::CalcNormal(m, fcu, fci, edim, fcn, fck);
}

template <class M_>
void UNormal<M_>::CalcNormalYoung(
    M& m, const FieldCell<Scal>& fcu, const FieldCell<bool>& fci,
    FieldCell<Vect>& fcn) {
  Imp::CalcNormalYoung(m, fcu, fci, fcn);
}

} // namespace solver
