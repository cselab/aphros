#pragma once

#include <vector>
#include <memory>
#include <limits>
#include <iomanip>
#include <set>
#include <map>
#include <algorithm>
#include <functional>
#include <stdio.h>

#include <march.h>

#include "vof.h"
#include "geom/mesh.h"
#include "dump/vtk.h"
#include "solver/reconst.h"
#include "debug/isnan.h"

namespace solver {

template <class T>
std::ostream& operator<<(std::ostream& o, const std::vector<T>& v) {
  std::string p = "";
  for (auto a : v) {
    o << p << a;
    p = " ";
  }
  return o;
}

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
    (void) th;
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
          if ((*fci[i])[c]) {
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
        WriteVtkPoly<Vect>(fn, dl_, nullptr, 
            {&dlc_, &dll_, &dlcl_}, {"c", "l", "cl"},
            "Interface from PLIC", true,
            bin, merge);
      }
    }
  }

  // uu: volume fraction in nodes
  // xc: cell center
  // h: cell size
  // iso: isovalue for surface uu=iso
  void GetMarchTriangles(
      const std::array<Scal, 8>& uu, const std::array<Vect, 8>& nn,
      const Vect& xc, 
      const Vect& h, Scal iso,
      std::vector<std::vector<Vect>>& vv,
      std::vector<std::vector<Vect>>& vvn) {
    std::array<double, 8> uuz = uu;
    for (auto& u : uuz) {
      u -= iso;
    }
    int nt;
    const int ms = MARCH_NTRI;
    std::array<double, ms> tri;
    march_cube(uuz.data(), &nt, tri.data());
    std::array<int, ms> vc0;
    std::array<int, ms> vc1;
    std::array<double, ms> vw;
    march_cube_location(vc0.data(), vc1.data(), vw.data());
    assert(nt * 3 * 3 <= tri.size());

    vv.resize(nt);
    {
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
    }
    vvn.resize(nt);
    {
      size_t i = 0;
      for (auto& vn : vvn) {
        vn.resize(3);
        for (auto& n : vn) {
          Scal w = vw[i];
          int c0 = vc0[i];
          int c1 = vc1[i];
          ++i;
          n = nn[c0] * (1 - w) + nn[c1] * w;
        }
      }
    }
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
  // Interpolates from cells to nodes.
  // stencil half-width
  template <int sw, int sn=sw*2+1, int snn=sw*2>
  std::array<Vect, snn*snn*snn> GradientNodes(const std::array<Scal, sn*sn*sn>& uu) {
    std::array<Vect, snn*snn*snn> gg;
    size_t i = 0;
    for (int z = 0; z < snn; ++z) {
      for (int y = 0; y < snn; ++y) {
        for (int x = 0; x < snn; ++x) {
          auto u = [&uu,x,y,z](int dx, int dy, int dz) {
            return uu[(z+dz)*sn*sn + (y+dy)*sn + (x+dx)];
          };
          auto& g = gg[i++];
          g[0] = ((u(1,0,0)+ u(1,1,0)+ u(1,0,1)+ u(1,1,1)) -
                   (u(0,0,0)+ u(0,1,0)+ u(0,0,1)+ u(0,1,1))) * 0.25;
          g[1] = ((u(0,1,0)+ u(1,1,0)+ u(0,1,1)+ u(1,1,1)) -
                   (u(0,0,0)+ u(1,0,0)+ u(0,0,1)+ u(1,0,1))) * 0.25;
          g[2] = ((u(0,0,1)+ u(1,0,1)+ u(0,1,1)+ u(1,1,1)) -
                   (u(0,0,0)+ u(1,0,0)+ u(0,1,0)+ u(1,1,0))) * 0.25;
        }
      }
    }
    return gg;
  }
    // Returns values over stencil centered at cell c with color cl.
    // Values for neighbors without color cl are filled with 0.
    // sw: stencil half-width
    template <size_t sw>
    struct GetStencilPure {
      static constexpr size_t sn = sw * 2 + 1;
      std::array<Scal, sn*sn*sn> operator()(
          const FieldCell<Scal>& fc, IdxCell c, const M& m) {
        using MIdx = typename M::MIdx;
        auto& bc = m.GetIndexCells();
        GBlock<IdxCell, M::dim> bo(MIdx(-sw), MIdx(sn));
        MIdx w = bc.GetMIdx(c);
        std::array<typename M::Scal, sn*sn*sn> uu;
        size_t k = 0;
        for (MIdx wo : bo) {
          IdxCell cn = bc.GetIdx(w + wo);
          uu[k++] = fc[cn];
        }
        return uu;
      }
    };
  // bcfill: if >=0. add triangles from 
  // fcus [a]: sum of volume fractions, add triangles from SuCells if not null
  void DumpPolyMarch(
      const GRange<size_t>& layers,
      const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<const FieldCell<Scal>*>& fccl,
      const Multi<const FieldCell<Vect>*>& fcn,
      const Multi<const FieldCell<Scal>*>& fca,
      const Multi<const FieldCell<bool>*>& fci,
      std::string fn, Scal t, Scal th, bool bin, bool merge, Scal iso, 
      const FieldCell<Scal>* fcus, M& m) {
    (void) fcn;
    (void) fca;
    (void) fci;
    (void) th;
    auto sem = m.GetSem("dumppolymarch");
    if (sem("local")) {
      dl_.clear();
      dln_.clear();
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
                auto uu = GetStencil<M, 1>{}(layers, fcu, fccl, cn, cl, m);
                auto uun = ToNodes<1>(uu);
                auto nn = GradientNodes<1>(uu);
                for (auto& n : nn) {
                  n = -n;
                }
                std::vector<std::vector<Vect>> vv, vvn;
                GetMarchTriangles(uun, nn, m.GetCenter(cn), h, iso, vv, vvn);
                for (size_t j = 0; j < vv.size(); ++j) {
                  dl_.push_back(vv[j]);
                  dln_.push_back(vvn[j]);
                  dlc_.push_back(m.GetHash(cn));
                  dll_.push_back(i);
                  dlcl_.push_back(cl);
                }
              }
            }
          }
        }
      }
      if (fcus) {
        // Append triangles on boundaries
        auto& bc = m.GetIndexCells();
        using MIdx = typename M::MIdx;
        for (auto c : m.SuCells()) {
          MIdx w = bc.GetMIdx(c);
          if (!(MIdx(1) <= w && w < m.GetGlobalSize() - MIdx(1))) {
            auto uu = GetStencilPure<1>{}(*fcus, c, m);
            auto uun = ToNodes<1>(uu);
            auto nn = GradientNodes<1>(uu);
            for (auto& n : nn) { n = -n; }
            std::vector<std::vector<Vect>> vv, vvn;
            GetMarchTriangles(uun, nn, m.GetCenter(c), h, iso, vv, vvn);
            for (size_t j = 0; j < vv.size(); ++j) {
              dl_.push_back(vv[j]);
              dln_.push_back(vvn[j]);
              dlc_.push_back(m.GetHash(c));
              dll_.push_back(-1);
              dlcl_.push_back(-1);
            }
          }
        }
      }
      using TV = typename M::template OpCatVT<Vect>;
      m.Reduce(std::make_shared<TV>(&dl_));
      m.Reduce(std::make_shared<TV>(&dln_));
      using TS = typename M::template OpCatT<Scal>;
      m.Reduce(std::make_shared<TS>(&dlc_));
      m.Reduce(std::make_shared<TS>(&dll_));
      m.Reduce(std::make_shared<TS>(&dlcl_));
    }
    if (sem("write")) {
      if (m.IsRoot()) {
        std::cout << std::fixed << std::setprecision(8)
            << "dump" << " t=" << t << " to " << fn << std::endl;
        WriteVtkPoly<Vect>(fn, dl_, &dln_,
            {&dlc_, &dll_, &dlcl_}, {"c", "l", "cl"}, 
            "Interface from marching cubes", true, bin, merge);
      }
    }
  }

  // Initializes usermap_.
  // fccl0: known colors
  // fccl: colors to reduce
  static void UserMap(const GRange<size_t>& layers,
               const Multi<const FieldCell<Scal>*>& fccl0,
               const Multi<const FieldCell<Scal>*>& fccl,
               std::map<Scal, Scal>& usermap, M& m) {
    auto sem = m.GetSem("usermap");
    struct {
      std::vector<Scal> vcl, vcln;
    }* ctx(sem);
    auto& vcl = ctx->vcl;
    auto& vcln = ctx->vcln;
    if (sem("local")) {
      usermap.clear();
      for (auto c : m.AllCells()) {
        for (auto i : layers) {
          Scal a = (*fccl0[i])[c];
          Scal cl = (*fccl[i])[c];
          if (a != kClNone && cl != kClNone) {
            if (!usermap.count(cl)) {
              usermap[cl] = a;
            }
          }
        }
      }
      vcl.clear();
      vcln.clear();
      for (auto p : usermap) {
        vcl.push_back(p.first);
        vcln.push_back(p.second);
      }
      using T = typename M::template OpCatT<Scal>;
      m.Reduce(std::make_shared<T>(&vcl));
      m.Reduce(std::make_shared<T>(&vcln));
    }
    if (sem("gather")) {
      usermap.clear();
      if (m.IsRoot()) {
        for (size_t i = 0; i < vcl.size(); ++i) {
          usermap[vcl[i]] = vcln[i];
        }
      }
    }
  }

  // Applies grid heuristic 
  static void Grid(const GRange<size_t>& layers,
            const Multi<const FieldCell<Scal>*>& fccl,
            const Multi<FieldCell<Scal>*>& fcclt, M& m) {
    auto sem = m.GetSem("grid");
    struct {
      std::vector<Scal> merge0, merge1; // colors to merge
    }* ctx(sem);
    auto& merge0 = ctx->merge0;
    auto& merge1 = ctx->merge1;
    if (sem("local")) {
      // Collect neighbor colors in corners
      merge0.clear(); // nodes
      merge1.clear(); // parents
      for (auto c : m.Cells()) {
        for (auto l : layers) {
          if ((*fccl[l])[c] != kClNone) {
            for (auto q : {0, 1, 2}) {
              IdxCell cm = m.GetCell(c, q);
              for (auto lm : layers) {
                if ((*fccl[l])[c] == (*fccl[lm])[cm]) {
                  Scal cl = (*fcclt[l])[c];
                  Scal clm = (*fcclt[lm])[cm];
                  if (cl != clm) {
                    merge0.push_back(std::max(cl, clm));
                    merge1.push_back(std::min(cl, clm));
                  }
                }
              }
            }
          }
        }
        break;
      }
      using T = typename M::template OpCatT<Scal>;
      m.Reduce(std::make_shared<T>(&merge0));
      m.Reduce(std::make_shared<T>(&merge1));
    }
    if (sem("reduce")) {
      if (m.IsRoot()) {
        std::map<Scal, Scal> map;

        for (size_t i = 0; i < merge0.size(); ++i) {
          map[merge0[i]] = merge1[i];
        }
        int iter = 0;
        // Find minimal color connected through pairs
        while (true) {
          bool chg = false;
          for (auto& p : map) {
            if (map.count(p.second)) {
              p.second = map[p.second];
              chg = true;
            }
          }
          if (!chg) {
            break;
          }
          ++iter;
        }
        merge0.clear();
        merge1.clear();
        for (auto& p : map) {
          merge0.push_back(p.first);
          merge1.push_back(p.second);
        }
      }
      using T = typename M::template OpCatT<Scal>;
      m.Bcast(std::make_shared<T>(&merge0));
      m.Bcast(std::make_shared<T>(&merge1));
    }
    if (sem("apply")) {
      std::map<Scal, Scal> map;
      for (size_t i = 0; i < merge0.size(); ++i) {
        map[merge0[i]] = merge1[i];
      }
      for (auto f : m.Faces()) {  // FIXME: inner cells traversed twice
        for (size_t q : {0, 1}) {
          auto c = m.GetCell(f, q);
          for (auto l : layers) {
            auto& cl = (*fcclt[l])[c];
            if (map.count(cl)) {
              cl = map[cl];
            }
          }
        }
      }
    }
  }

  // Reduces the color space.
  // usermap: suggested map <old,new>
  //          applies maximum subset of usermap that produces unique colors,
  //          replaces others with new colors.
  static void ReduceColor(const GRange<size_t>& layers,
                   const Multi<FieldCell<Scal>*>& fccl,
                   const std::map<Scal, Scal>& usermap, Scal clfixed, M& m) {
    auto sem = m.GetSem("recolor");
    struct {
      std::vector<std::vector<Scal>> vvcl; // all colors
      std::vector<Scal> vcl, vcln;
    }* ctx(sem);
    auto& vvcl = ctx->vvcl;
    auto& vcl = ctx->vcl;
    // gather all colors from domain
    if (sem("gather")) {
      std::set<Scal> s;
      for (auto l : layers) {
        for (auto c : m.AllCells()) {
          auto& cl = (*fccl[l])[c];
          if (cl != kClNone) {
            s.insert(cl);
          }
        }
      }
      vcl = std::vector<Scal>(s.begin(), s.end());
      vvcl = {vcl};
      using T = typename M::template OpCatVT<Scal>;
      m.Reduce(std::make_shared<T>(&vvcl));
    }
    // replace with reduced set, applying usermap if possible
    auto& vcln = ctx->vcln;
    if (sem("reduce")) {
      if (m.IsRoot()) {
        std::map<Scal, Scal> map;
        std::set<Scal> used;
        // keep some colors
        auto Add = [&](Scal cl, Scal a) {
          map[cl] = a;
          used.insert(a);
        };
        Add(kClNone, kClNone);
        if (clfixed >= 0) {
          Add(clfixed, clfixed);
        }
        Scal an = 0; // next unused color
        for (auto p : usermap) {
          an = std::max(an, p.second);
        }
        an += 1;
        for (auto& v : vvcl) {
          for (auto& cl : v) {
            if (!map.count(cl)) {
              //map[cl] = usermap.count(cl) ? usermap.at(cl) : -2;
              if (!usermap.count(cl) || used.count(usermap.at(cl))) {
                Add(cl, an);
                an += 1;
              } else { // usermap.count(cl) && !used.count(usermap[cl])
                Add(cl, usermap.at(cl));
              }
            }
            cl = map[cl];
          }
        }
        m.Scatter({&vvcl, &vcln});
      } else {
        m.Scatter({nullptr, &vcln});
      }
    }
    // apply the new set from vcln
    if (sem("apply")) {
      std::map<Scal, Scal> map;
      for (size_t i = 0; i < vcl.size(); ++i) {
        map[vcl[i]] = vcln[i];
      }
      for (auto c : m.AllCells()) {
        for (auto l : layers) {
          auto& cl = (*fccl[l])[c];
          if (cl != kClNone) {
            assert(map.count(cl));
            cl = map[cl];
          }
        }
      }
    }
  }

  static void Init(const GRange<size_t>& layers,
      const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<FieldCell<Scal>*>& fccl,
      Multi<FieldCell<Scal>>& fcclt,
      Scal clfixed, Vect clfixed_x, Scal coalth, M& m) {
    auto sem = m.GetSem("recolor_init");
    struct {
      std::pair<typename M::Scal, int> cldist; // color,mesh_id
    }* ctx(sem);
    auto& cldist = ctx->cldist;
    if (sem("clfixed")) {
      // block nearest to clfixed_x
      if (clfixed >= 0) {
        IdxCell c = m.FindNearestCell(clfixed_x);
        cldist.first = m.GetCenter(c).dist(clfixed_x);
        cldist.second = m.GetId();
        m.Reduce(std::make_shared<typename M::OpMinloc>(&cldist));
      } else {
        cldist.second = -1;
      }
    }
    if (sem("init")) {
      fcclt.resize(layers.size());
      fcclt.InitAll(FieldCell<Scal>(m, kClNone));
      // initial unique color
      Scal q = m.GetId() * m.GetInBlockCells().size() * layers.size();
      for (auto i : layers) {
        for (auto c : m.Cells()) {
          if ((*fccl[i])[c] != kClNone) {
            fcclt[i][c] = (q += 1);
          }
        }
        if (cldist.second == m.GetId()) {
          IdxCell c = m.FindNearestCell(clfixed_x);
          if ((*fccl[i])[c] != kClNone) {
            fcclt[i][c] = clfixed;
          }
        }
        m.Comm(&fcclt[i]);
      }

      // detect overlap
      for (auto i : layers) {
        for (auto c : m.Cells()) {
          if ((*fccl[i])[c] != kClNone) {
            for (auto j : layers) {
              if (j != i && (*fccl[j])[c] != kClNone) {
                if ((*fcu[i])[c] + (*fcu[j])[c] > coalth) {
                  Scal cl = std::min(fcclt[i][c], fcclt[j][c]);
                  fcclt[i][c] = cl;
                  fcclt[j][c] = cl;
                }
              }
            }
          }
        }
      }
    }
  }

  static void RecolorDirect(const GRange<size_t>& layers,
      const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<FieldCell<Scal>*>& fccl,
      const Multi<const FieldCell<Scal>*>& fccl0,
      Scal clfixed, Vect clfixed_x, Scal coalth,
      const MapFace<std::shared_ptr<CondFace>>& mfc,
      bool bcc_reflect, bool verb, bool reduce, bool grid, M& m) {
    auto sem = m.GetSem("recolor");
    struct {
      std::map<Scal, Scal> usermap;
      Scal tries;
      Multi<FieldCell<Scal>> fcclt;  // tmp color
    }* ctx(sem);
    auto& fcclt = ctx->fcclt;
    if (sem.Nested()) {
      Init(layers, fcu, fccl, fcclt, clfixed, clfixed_x, coalth, m);
    }
    if (bcc_reflect && sem("reflect")) {
      for (auto i : layers) {
        BcReflect(fcclt[i], mfc, kClNone, false, m);
      }
    }
    sem.LoopBegin();
    if (grid && sem.Nested()) {
      Grid(layers, fccl, fcclt, m);
    }
    if (sem("min")) {
      size_t tries = 0;
      size_t cells = 0;
      using MIdx = typename M::MIdx;
      auto& bc = m.GetIndexCells();
      static constexpr size_t sw = 1;  // stencil half-width
      static constexpr size_t sn = sw * 2 + 1;
      GBlock<IdxCell, M::dim> bo(MIdx(-sw), MIdx(sn));
      while (true) {
        bool chg = false;
        for (auto i : layers) {
          for (auto c : m.Cells()) {
            if ((*fccl[i])[c] != kClNone) {
              // update color with minimum over neighbours
              MIdx w = bc.GetMIdx(c);
              for (MIdx wo : bo) {
                IdxCell cn = bc.GetIdx(w + wo);
                for (auto j : layers) {
                  if ((*fccl[j])[cn] == (*fccl[i])[c]) {
                    if (fcclt[j][cn] < fcclt[i][c]) {
                      chg = true;
                      ++cells;
                      fcclt[i][c] = fcclt[j][cn];
                    }
                  }
                }
              }
            }
          }
        }
        if (!chg) {
          break;
        }
        ++tries;
      }
      for (auto i : layers) {
        m.Comm(&fcclt[i]);
      }
      ctx->tries = tries;
      m.Reduce(&ctx->tries, "max");
    }
    if (bcc_reflect && sem("reflect")) {
      for (auto i : layers) {
        BcReflect(fcclt[i], mfc, kClNone, false, m);
      }
    }
    if (sem("check")) {
      if (verb && m.IsRoot()) {
        std::cerr << "recolor:"
          << " max tries: " << ctx->tries 
          << std::endl;
      }
      if (!ctx->tries) {
        sem.LoopBreak();
      }
    }
    sem.LoopEnd();
    if (reduce && sem.Nested()) {
      UserMap(layers, fccl0, fcclt, ctx->usermap, m);
    }
    if (reduce && sem.Nested()) {
      ReduceColor(layers, fcclt, ctx->usermap, clfixed, m);
    }
    if (sem("copy")) {
      for (auto c : m.AllCells()) {
        for (auto l : layers) {
          (*fccl[l])[c] = fcclt[l][c];
        }
      }
    }
  }

  static void RecolorUnionFind(const GRange<size_t>& layers,
      const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<FieldCell<Scal>*>& fccl,
      const Multi<const FieldCell<Scal>*>& fccl0,
      Scal clfixed, Vect clfixed_x, Scal coalth,
      const MapFace<std::shared_ptr<CondFace>>& mfc,
      bool bcc_reflect, bool verb, bool reduce, bool grid, M& m) {
    auto sem = m.GetSem("recolor");
    struct {
      std::map<Scal, Scal> usermap;
      Scal tries;
      Multi<FieldCell<Scal>> fcclt;  // tmp color
      Multi<FieldCell<IdxCell>> fcc; // root cell
      Multi<FieldCell<char>> fcl;  // root layer
    }* ctx(sem);
    auto& fcclt = ctx->fcclt;
    auto& fcc = ctx->fcc;
    auto& fcl = ctx->fcl;
    if (sem.Nested()) {
      Init(layers, fcu, fccl, fcclt, clfixed, clfixed_x, coalth, m);
    }
    if (bcc_reflect && sem("reflect")) {
      for (auto i : layers) {
        BcReflect(fcclt[i], mfc, kClNone, false, m);
      }
    }
    if (sem("initroot")) {
      fcc.resize(layers.size());
      fcl.resize(layers.size());
      fcc.InitAll(FieldCell<IdxCell>(m));
      fcl.InitAll(FieldCell<char>(m));
      for (auto c : m.AllCells()) {
        for (auto l : layers) {
          fcc[l][c] = c;
          fcl[l][c] = l;
        }
      }
    }
    sem.LoopBegin();
    if (grid && sem.Nested()) {
      Grid(layers, fccl, fcclt, m);
    }
    if (sem("min")) {
      using Pair = std::pair<IdxCell, char>;
      size_t tries = 0;
      // Returns reference to color
      auto Clt = [&](Pair p) -> Scal& {
        return fcclt[p.second][p.first];
      };
      // Returns parent of p.
      auto Get = [&](Pair p) {
        IdxCell c = p.first;
        char l = p.second;
        return Pair(fcc[l][c], fcl[l][c]);
      };
      // Sets parent of p to r.
      auto Set = [&](Pair p, Pair r) {
        IdxCell c = p.first;
        char l = p.second;
        fcc[l][c] = r.first;
        fcl[l][c] = r.second;
        return r;
      };
      // Returns root of p.
      std::function<Pair(Pair)> Find = [&](Pair p) {
        if (Get(p) == p) {
          return p;
        }
        return Set(p, Find(Get(p)));
      };
      // Merges sets with representatives p0 and p1
      auto Union = [&](Pair p0, Pair p1) {
        p0 = Find(p0);
        p1 = Find(p1);
        if (p0 != p1) {
          if (Clt(p0) <= Clt(p1)) {
            Set(p1, p0);
          } else {
            Set(p0, p1);
          }
          ++tries;
        }
      };

      using MIdx = typename M::MIdx;
      auto& bc = m.GetIndexCells();
      static constexpr size_t sw = 1;  // stencil half-width
      static constexpr size_t sn = sw * 2 + 1;
      GBlock<IdxCell, M::dim> bo(MIdx(-sw), MIdx(sn));

      // Merge all neighbors with the same color.
      for (auto c : m.Cells()) {
        MIdx w = bc.GetMIdx(c);
        for (auto l : layers) {
          if ((*fccl[l])[c] != kClNone) {
            for (MIdx wo : bo) {
              IdxCell cn = bc.GetIdx(w + wo);
              for (auto ln : layers) {
                if ((*fccl[l])[c] == (*fccl[ln])[cn]) {
                  Pair p(c, l);
                  Pair pn(cn, ln);
                  Union(pn, p);
                }
              }
            }
          }
        }
      }
      // Set all cells to their root
      for (auto c : m.SuCells()) {
        for (auto l : layers) {
          Pair p(c, l);
          Set(p, Find(p));
        }
      }
      // Update color in root
      // (smaller color can arrive in halo cells)
      for (auto c : m.SuCells()) {
        for (auto l : layers) {
          Pair p(c, l);
          if (Clt(p) < Clt(Get(p))) {
            Clt(Get(p)) = Clt(p);
            ++tries;
          }
        }
      }
      // Update color from root
      for (auto c : m.SuCells()) {
        for (auto l : layers) {
          Pair p(c, l);
          if (Clt(Get(p)) < Clt(p)) {
            Clt(p) = Clt(Get(p));
            ++tries;
          }
        }
      }
      for (auto l : layers) {
        m.Comm(&fcclt[l]);
      }
      ctx->tries = tries;
      m.Reduce(&ctx->tries, "max");
    }
    if (bcc_reflect && sem("reflect")) {
      for (auto i : layers) {
        BcReflect(fcclt[i], mfc, kClNone, false, m);
      }
    }
    if (sem("check")) {
      if (verb && m.IsRoot()) {
        std::cerr << "recolor:"
          << " max tries: " << ctx->tries 
          << std::endl;
      }
      if (!ctx->tries) {
        sem.LoopBreak();
      }
    }
    sem.LoopEnd();
    if (reduce && sem.Nested()) {
      UserMap(layers, fccl0, fcclt, ctx->usermap, m);
    }
    if (reduce && sem.Nested()) {
      ReduceColor(layers, fcclt, ctx->usermap, clfixed, m);
    }
    if (sem("copy")) {
      for (auto c : m.AllCells()) {
        for (auto l : layers) {
          (*fccl[l])[c] = fcclt[l][c];
        }
      }
    }
  }

  static void Recolor(const GRange<size_t>& layers,
      const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<FieldCell<Scal>*>& fccl,
      const Multi<const FieldCell<Scal>*>& fccl0,
      Scal clfixed, Vect clfixed_x, Scal coalth,
      const MapFace<std::shared_ptr<CondFace>>& mfc,
      bool bcc_reflect, bool verb, bool unionfind, bool reduce, bool grid,
      M& m) {
    if (unionfind) {
      return RecolorUnionFind(layers, fcu, fccl, fccl0, clfixed, clfixed_x,
                              coalth, mfc, bcc_reflect, verb, reduce, grid, m);
    }
    return RecolorDirect(layers, fcu, fccl, fccl0, clfixed, clfixed_x,
                        coalth, mfc, bcc_reflect, verb, reduce, grid, m);
  }

  std::vector<std::vector<Vect>> dl_; // dump poly
  std::vector<std::vector<Vect>> dln_; // dump poly normals
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
    std::string fn, Scal t, Scal th, bool bin, bool merge, Scal iso,
    const FieldCell<Scal>* fcus, M& m) {
  imp->DumpPolyMarch(
      layers, fcu, fccl, fcn, fca, fci, fn, t, th, bin, merge, iso, fcus, m);
}

template <class M_>
void UVof<M_>::Recolor(const GRange<size_t>& layers,
    const Multi<const FieldCell<Scal>*>& fcu,
    const Multi<FieldCell<Scal>*>& fccl,
    const Multi<const FieldCell<Scal>*>& fccl0,
    Scal clfixed, Vect clfixed_x, Scal coalth,
    const MapFace<std::shared_ptr<CondFace>>& mfcu,
    bool bcc_reflect, bool verb, bool unionfind, bool reduce, bool grid,
    M& m) {
  imp->Recolor(layers, fcu, fccl, fccl0, clfixed, clfixed_x, coalth, mfcu,
               bcc_reflect, verb, unionfind, reduce, grid, m);
}


template <class M>
constexpr typename M::Scal UVof<M>::Imp::kClNone;

} // namespace solver
