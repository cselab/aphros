#pragma once

#include <vector>
#include <memory>
#include <limits>
#include <iomanip>
#include <set>

#include <march.h>

#include "vof.h"
#include "geom/mesh.h"
#include "dump/vtk.h"
#include "solver/reconst.h"
#include "debug/isnan.h"

namespace solver {

template <class M_>
struct UVof<M_>::Imp {
  using R = Reconst<Scal>;
  static constexpr size_t dim = M::dim;
  static constexpr Scal kClNone = -1;

  Imp() = default;
  ~Imp() = default;

  void DumpPoly(
      const GRange<size_t>& layers,
      const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<const FieldCell<Scal>*>& fccl,
      const Multi<const FieldCell<Vect>*>& fcn,
      const Multi<const FieldCell<Scal>*>& fca,
      const Multi<const FieldCell<bool>*>& fci,
      std::string fn, Scal t, Scal th, bool bin, bool merge, M& m) {
    auto sem = m.GetSem("dumppoly");
    if (sem("local")) {
      dl_.clear();
      dlc_.clear();
      dll_.clear();
      dlcl_.clear();
      auto h = m.GetCellSize();
      for (auto i : layers) {
        for (auto c : m.Cells()) {
          Scal u = (*fcu[i])[c];
          if (IsNan(u) || IsNan((*fcn[i])[c]) || IsNan((*fca[i])[c])) {
            continue;
          }
          if ((*fci[i])[c] && u > th && u < 1 - th) {
            dl_.push_back(R::GetCutPoly(
                  m.GetCenter(c), (*fcn[i])[c], (*fca[i])[c], h));
            dlc_.push_back(m.GetHash(c));
            dll_.push_back(i);
            dlcl_.push_back(fccl[i] ? (*fccl[i])[c] : 0);
          }
        }
      }
      using TV = typename M::template OpCatVT<Vect>;
      m.Reduce(std::make_shared<TV>(&dl_));
      using TS = typename M::template OpCatT<Scal>;
      m.Reduce(std::make_shared<TS>(&dlc_));
      m.Reduce(std::make_shared<TS>(&dll_));
      m.Reduce(std::make_shared<TS>(&dlcl_));
    }
    if (sem("write")) {
      if (m.IsRoot()) {
        std::cout << std::fixed << std::setprecision(8)
            << "dump" << " t=" << t << " to " << fn << std::endl;
        WriteVtkPoly(fn, dl_, {&dlc_, &dll_, &dlcl_}, {"c", "l", "cl"},
            "Interface from PLIC", true,
            bin, merge);
      }
    }
  }

  // uu: volume fraction in nodes
  // xc: cell center
  // h: cell size
  // iso: isovalue for surface uu=iso
  std::vector<std::vector<Vect>> GetMarchTriangles(
      const std::array<Scal, 8>& uu, const Vect& xc, 
      const Vect& h, Scal iso) {
    std::array<double, 8> uuz = uu;
    for (auto& u : uuz) {
      u -= iso;
    }
    int n;
    std::array<double, 45> tri;
    march_cube(uuz.data(), &n, tri.data());
    assert(n * 3 * 3 <= tri.size());
    std::vector<std::vector<Vect>> vv;

    vv.resize(n);
    size_t i = 0;
    for (auto& v : vv) {
      v.resize(3);
      for (auto& x : v) {
        x[0] = tri[i++] - 0.5;
        x[1] = tri[i++] - 0.5;
        x[2] = tri[i++] - 0.5;
        x = xc + h * x;
      }
    }
    return vv;
  }
  // Returns values over stencil centered at cell c with color cl.
  // Values for neighbors without color cl are filled with 0.
  // sw: stencil half-width
  template <size_t sw, size_t sn=sw*2+1>
  std::array<Scal, sn*sn*sn> GetStencil(
      const GRange<size_t>& layers,
      const Multi<const FieldCell<Scal>*>& fc,
      const Multi<const FieldCell<Scal>*>& fccl,
      IdxCell c, Scal cl, M& m) {
    using MIdx = typename M::MIdx;
    auto& bc = m.GetIndexCells();
    GBlock<IdxCell, dim> bo(MIdx(-sw), MIdx(sn));
    MIdx w = bc.GetMIdx(c);
    std::array<Scal, sn*sn*sn> uu;
    size_t k = 0;
    for (MIdx wo : bo) {
      IdxCell cn = bc.GetIdx(w + wo);
      Scal u = 0;
      for (auto j : layers) {
        if ((*fccl[j])[cn] == cl) {
          u = (*fc[j])[cn];
          break;
        }
      }
      uu[k++] = u;
    }
    return uu;
  }
  // Interpolates from cells to nodes.
  // stencil half-width
  template <int sw, int sn=sw*2+1, int snn=sw*2>
  std::array<Scal, snn*snn*snn> ToNodes(const std::array<Scal, sn*sn*sn>& uu) {
    std::array<Scal, snn*snn*snn> uun;
    size_t i = 0;
    for (int z = 0; z < snn; ++z) {
      for (int y = 0; y < snn; ++y) {
        for (int x = 0; x < snn; ++x) {
          auto u = [&uu,x,y,z](int dx, int dy, int dz) {
            return uu[(z+dz)*sn*sn + (y+dy)*sn + (x+dx)];
          };
          uun[i++] = (1. / 8.) * (
              u(0,0,0) + u(1,0,0) + u(0,1,0) + u(1,1,0) +
              u(0,0,1) + u(1,0,1) + u(0,1,1) + u(1,1,1));
        }
      }
    }
    return uun;
  }
  void DumpPolyMarch(
      const GRange<size_t>& layers,
      const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<const FieldCell<Scal>*>& fccl,
      const Multi<const FieldCell<Vect>*>& fcn,
      const Multi<const FieldCell<Scal>*>& fca,
      const Multi<const FieldCell<bool>*>& fci,
      std::string fn, Scal t, Scal th, bool bin, bool merge, Scal iso, M& m) {
    (void) fcn;
    (void) fca;
    (void) fci;
    (void) th;
    auto sem = m.GetSem("dumppolymarch");
    if (sem("local")) {
      dl_.clear();
      dlc_.clear();
      dll_.clear();
      dlcl_.clear();
      // visited cells without color
      std::set<std::pair<size_t, Scal>> done;
      auto h = m.GetCellSize();
      FieldCell<bool> in(m, false); // inner cell
      for (auto c : m.Cells()) {
        in[c] = true;
      }
      for (auto i : layers) {
        for (auto c : m.Cells()) {
          Scal cl = (*fccl[i])[c];
          if (cl != kClNone) {
            const size_t sw = 1; // stencil halfwidth
            const size_t sn = sw * 2 + 1;
            auto& bc = m.GetIndexCells();
            using MIdx = typename M::MIdx;
            GBlock<IdxCell, dim> bo(MIdx(-sw), MIdx(sn));
            MIdx w = bc.GetMIdx(c);
            for (MIdx wo : bo) {
              IdxCell cn = bc.GetIdx(w + wo);
              if (!in[cn]) continue;
              bool q = false;
              for (auto j : layers) {
                if ((*fccl[j])[cn] == cl) {
                  q = true;
                  break;
                }
              }
              // Add triangles from c or neighhbor without color
              if (c == cn || !q) {
                auto e = std::make_pair(size_t(cn), cl);
                if (!q && done.count(e)) {
                  continue;
                }
                done.insert(e);
                auto uu = GetStencil<1>(layers, fcu, fccl, cn, cl, m);
                auto uun = ToNodes<1>(uu);
                auto vv = GetMarchTriangles(uun, m.GetCenter(cn), h, iso);
                for (auto& v : vv) {
                  dl_.push_back(v);
                  dlc_.push_back(m.GetHash(cn));
                  dll_.push_back(i);
                  dlcl_.push_back(cl);
                }
              }
            }
          }
        }
      }
      using TV = typename M::template OpCatVT<Vect>;
      m.Reduce(std::make_shared<TV>(&dl_));
      using TS = typename M::template OpCatT<Scal>;
      m.Reduce(std::make_shared<TS>(&dlc_));
      m.Reduce(std::make_shared<TS>(&dll_));
      m.Reduce(std::make_shared<TS>(&dlcl_));
    }
    if (sem("write")) {
      if (m.IsRoot()) {
        std::cout << std::fixed << std::setprecision(8)
            << "dump" << " t=" << t << " to " << fn << std::endl;
        WriteVtkPoly(fn, dl_, {&dlc_, &dll_, &dlcl_}, {"c", "l", "cl"}, 
            "Interface from marching cubes", true, bin, merge);
      }
    }
  }

  std::vector<std::vector<Vect>> dl_; // dump poly
  std::vector<Scal> dlc_; // dump poly cell
  std::vector<Scal> dll_; // dump poly layer
  std::vector<Scal> dlcl_; // dump poly color
};

template <class M_>
UVof<M_>::UVof() : imp(new Imp()) {}

template <class M_>
UVof<M_>::~UVof() = default;

template <class M_>
void UVof<M_>::DumpPoly(
    const GRange<size_t>& layers,
    const Multi<const FieldCell<Scal>*>& fcu,
    const Multi<const FieldCell<Scal>*>& fccl,
    const Multi<const FieldCell<Vect>*>& fcn,
    const Multi<const FieldCell<Scal>*>& fca,
    const Multi<const FieldCell<bool>*>& fci,
    std::string fn, Scal t, Scal th, bool bin, bool merge, M& m) {
  imp->DumpPoly(layers, fcu, fccl, fcn, fca, fci, fn, t, th, bin, merge, m);
}

template <class M_>
void UVof<M_>::DumpPoly(
    const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcn,
    const FieldCell<Scal>& fca, const FieldCell<bool>& fci,
    std::string fn, Scal t, Scal th, bool bin, bool merge, M& m) {
  GRange<size_t> layers(0, 1);
  const FieldCell<Scal>* fccl(nullptr);
  imp->DumpPoly(layers, &fcu, fccl, &fcn, &fca, &fci,
                fn, t, th, bin, merge, m);
}

template <class M_>
void UVof<M_>::DumpPolyMarch(
    const GRange<size_t>& layers,
    const Multi<const FieldCell<Scal>*>& fcu,
    const Multi<const FieldCell<Scal>*>& fccl,
    const Multi<const FieldCell<Vect>*>& fcn,
    const Multi<const FieldCell<Scal>*>& fca,
    const Multi<const FieldCell<bool>*>& fci,
    std::string fn, Scal t, Scal th, bool bin, bool merge, Scal iso, M& m) {
  imp->DumpPolyMarch(
      layers, fcu, fccl, fcn, fca, fci, fn, t, th, bin, merge, iso, m);
}

template <class M>
constexpr typename M::Scal UVof<M>::Imp::kClNone;

} // namespace solver
