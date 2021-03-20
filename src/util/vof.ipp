// Created by Petr Karnakov on 02.09.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <stdio.h>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <vector>

#include "debug/isnan.h"
#include "dump/vtk.h"
#include "geom/mesh.h"
#include "march/march.h"
#include "solver/approx.h"
#include "solver/reconst.h"
#include "solver/trackerm.h"
#include "vof.h"

template <class T>
std::ostream& operator<<(std::ostream& o, const std::vector<T>& v) {
  std::string p = "";
  for (auto a : v) {
    o << p << a;
    p = " ";
  }
  return o;
}

template <class T>
void Reorder(std::vector<T>& v, const std::vector<size_t> idx) {
  std::vector<T> t = v;
  for (size_t i = 0; i < v.size(); ++i) {
    v[i] = t[idx[i]];
  }
}

template <class M_>
struct UVof<M_>::Imp {
  using R = Reconst<Scal>;
  static constexpr size_t dim = M::dim;
  static constexpr Scal kClNone = -1;
  using Vect2 = generic::Vect<Scal, 2>;
  using Vect3 = generic::Vect<Scal, 3>;
  using Vect4 = generic::Vect<Scal, 4>;

  Imp() = default;
  ~Imp() = default;

  static void DumpPoly(
      const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<const FieldCell<Scal>*>& fccl,
      const Multi<const FieldCell<Vect>*>& fcn,
      const Multi<const FieldCell<Scal>*>& fca,
      const Multi<const FieldCell<bool>*>& fci, std::string fn, Scal t,
      bool bin, bool merge, M& m) {
    auto sem = m.GetSem("dumppoly");
    struct {
      std::vector<std::vector<Vect>> dl; // polygons
      std::vector<Scal> dlc; // cells index
      std::vector<Scal> dll; // layer
      std::vector<Scal> dlcl; // color
    } * ctx(sem);
    auto& dl = ctx->dl;
    auto& dlc = ctx->dlc;
    auto& dll = ctx->dll;
    auto& dlcl = ctx->dlcl;
    if (sem("local")) {
      auto h = m.GetCellSize();
      for (auto i : layers) {
        for (auto c : m.Cells()) {
          Scal u = (*fcu[i])[c];
          if (IsNan(u) || IsNan((*fcn[i])[c]) || IsNan((*fca[i])[c])) {
            continue;
          }
          if ((*fci[i])[c]) {
            dl.push_back(
                R::GetCutPoly(m.GetCenter(c), (*fcn[i])[c], (*fca[i])[c], h));
            dlc.push_back(m.GetHash(c));
            dll.push_back(i);
            dlcl.push_back(fccl[i] ? (*fccl[i])[c] : 0);
          }
        }
      }
      m.Reduce(&dl, Reduction::concat);
      m.Reduce(&dlc, Reduction::concat);
      m.Reduce(&dll, Reduction::concat);
      m.Reduce(&dlcl, Reduction::concat);
    }
    if (sem("write")) {
      if (m.IsRoot()) {
        std::vector<size_t> idx(dlc.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::stable_sort(idx.begin(), idx.end(), [&](size_t i1, size_t i2) {
          return dlc[i1] < dlc[i2];
        });

        Reorder(dl, idx);
        Reorder(dlc, idx);
        Reorder(dll, idx);
        Reorder(dlcl, idx);

        std::cerr << std::fixed << std::setprecision(8) << "dump"
                  << " t=" << t << " to " << fn << std::endl;
        WriteVtkPoly<Vect>(
            fn, dl, nullptr, {&dlc, &dll, &dlcl}, {"c", "l", "cl"},
            "Interface from PLIC", true, bin, merge);
      }
    }
  }

  // Constructs iso-surface triangles with marching cubes.
  // uu: volume fraction in cell nodes, 3D
  // xc: cell center
  // h: cell size
  // iso: isovalue for surface uu=iso
  static void GetMarchTriangles(
      const std::array<Scal, 8>& uu, const std::array<Vect3, 8>& nn,
      const Vect3& xc, const Vect3& h, Scal iso, /*out*/
      std::vector<std::vector<Vect3>>& vv,
      std::vector<std::vector<Vect3>>& vvn) {
    (void)MARCH_O[0][0];
    std::array<double, 8> uuz = uu;
    for (auto& u : uuz) {
      u -= iso;
    }
    int nt;
    constexpr int kMaxNt = MARCH_NTRI;
    std::array<double, kMaxNt> tri;
    std::array<int, kMaxNt> vc0;
    std::array<int, kMaxNt> vc1;
    std::array<double, kMaxNt> vw;
    march_cube_location(
        uuz.data(), &nt, tri.data(), vc0.data(), vc1.data(), vw.data());
    assert(size_t(nt) * 3 * 3 <= tri.size());

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
  // 2D
  static void GetMarchTriangles(
      const std::array<Scal, 4>&, const std::array<Vect2, 4>&, const Vect2&,
      const Vect2&, Scal, std::vector<std::vector<Vect2>>&,
      std::vector<std::vector<Vect2>>&) {}
  // 4D
  static void GetMarchTriangles(
      const std::array<Scal, 16>&, const std::array<Vect4, 16>&, const Vect4&,
      const Vect4&, Scal, std::vector<std::vector<Vect4>>&,
      std::vector<std::vector<Vect4>>&) {}
  // Interpolates from cells to nodes
  // uu: values on a 3x3x3x3 stencil
  template <class T>
  static std::array<T, 16> InterpolateToNodes(const std::array<T, 81>& uu) {
    std::array<T, 16> res;
    size_t i = 0;
    for (int w = 0; w < 2; ++w) {
      for (int z = 0; z < 2; ++z) {
        for (int y = 0; y < 2; ++y) {
          for (int x = 0; x < 2; ++x) {
            auto u = [&uu, x, y, z, w](int dx, int dy, int dz, int dw) {
              return uu[(w + dw) * 27 + (z + dz) * 9 + (y + dy) * 3 + (x + dx)];
            };
            res[i++] =
                (1. / 16) *
                (u(0, 0, 0, 0) + u(1, 0, 0, 0) + u(0, 1, 0, 0) + u(1, 1, 0, 0) +
                 u(0, 0, 1, 0) + u(1, 0, 1, 0) + u(0, 1, 1, 0) + u(1, 1, 1, 0) +
                 u(0, 0, 0, 1) + u(1, 0, 0, 1) + u(0, 1, 0, 1) + u(1, 1, 0, 1) +
                 u(0, 0, 1, 1) + u(1, 0, 1, 1) + u(0, 1, 1, 1) + u(1, 1, 1, 1));
          }
        }
      }
    }
    return res;
  }
  // TODO: move to an util file
  // Interpolates from cells to nodes
  // uu: values on a 3x3x3 stencil
  template <class T>
  static std::array<T, 8> InterpolateToNodes(const std::array<T, 27>& uu) {
    std::array<T, 8> res;
    size_t i = 0;
    for (int z = 0; z < 2; ++z) {
      for (int y = 0; y < 2; ++y) {
        for (int x = 0; x < 2; ++x) {
          auto u = [&uu, x, y, z](int dx, int dy, int dz) {
            return uu[(z + dz) * 9 + (y + dy) * 3 + (x + dx)];
          };
          res[i++] =
              (1. / 8) * (u(0, 0, 0) + u(1, 0, 0) + u(0, 1, 0) + u(1, 1, 0) +
                          u(0, 0, 1) + u(1, 0, 1) + u(0, 1, 1) + u(1, 1, 1));
        }
      }
    }
    return res;
  }
  // Interpolates from cells to nodes
  // uu: values on a 3x3 stencil
  template <class T>
  static std::array<T, 4> InterpolateToNodes(const std::array<T, 9>& uu) {
    std::array<T, 4> res;
    size_t i = 0;
    for (int y = 0; y < 2; ++y) {
      for (int x = 0; x < 2; ++x) {
        auto u = [&uu, x, y](int dx, int dy) {
          return uu[(y + dy) * 3 + (x + dx)];
        };
        res[i++] = (1. / 4) * (u(0, 0) + u(1, 0) + u(0, 1) + u(1, 1));
      }
    }
    return res;
  }
  // TODO: move to an util file
  // Computes sum of cell values adjacent to each node, skips nan values.
  // sw: stencil half-width
  // uu: cell values on 3x3x3x3 stencil
  template <class T>
  static std::array<T, 16> SumToNodesNan(const std::array<T, 81>& uu) {
    std::array<T, 16> res;
    size_t i = 0;
    for (int w = 0; w < 2; ++w) {
      for (int z = 0; z < 2; ++z) {
        for (int y = 0; y < 2; ++y) {
          for (int x = 0; x < 2; ++x) {
            auto u = [&uu, x, y, z, w](int dx, int dy, int dz, int dw) {
              auto q =
                  uu[(w + dw) * 27 + (z + dz) * 9 + (y + dy) * 3 + (x + dx)];
              return IsNan(q) ? T(0) : q;
            };
            res[i++] =
                u(0, 0, 0, 0) + u(1, 0, 0, 0) + u(0, 1, 0, 0) + u(1, 1, 0, 0) +
                u(0, 0, 1, 0) + u(1, 0, 1, 0) + u(0, 1, 1, 0) + u(1, 1, 1, 0) +
                u(0, 0, 0, 1) + u(1, 0, 0, 1) + u(0, 1, 0, 1) + u(1, 1, 0, 1) +
                u(0, 0, 1, 1) + u(1, 0, 1, 1) + u(0, 1, 1, 1) + u(1, 1, 1, 1);
          }
        }
      }
    }
    return res;
  }
  // Computes sum of cell values adjacent to each node, skips nan values.
  // sw: stencil half-width
  // uu: cell values on 3x3x3 stencil
  template <class T>
  static std::array<T, 8> SumToNodesNan(const std::array<T, 27>& uu) {
    std::array<T, 8> res;
    size_t i = 0;
    for (int z = 0; z < 2; ++z) {
      for (int y = 0; y < 2; ++y) {
        for (int x = 0; x < 2; ++x) {
          auto u = [&uu, x, y, z](int dx, int dy, int dz) {
            auto q = uu[(z + dz) * 9 + (y + dy) * 3 + (x + dx)];
            return IsNan(q) ? T(0) : q;
          };
          res[i++] = u(0, 0, 0) + u(1, 0, 0) + u(0, 1, 0) + u(1, 1, 0) +
                     u(0, 0, 1) + u(1, 0, 1) + u(0, 1, 1) + u(1, 1, 1);
        }
      }
    }
    return res;
  }
  // Computes sum of cell values adjacent to each node, skips nan values.
  // sw: stencil half-width
  // uu: cell values on 3x3 stencil
  template <class T>
  static std::array<T, 4> SumToNodesNan(const std::array<T, 9>& uu) {
    std::array<T, 4> res;
    size_t i = 0;
    for (int y = 0; y < 2; ++y) {
      for (int x = 0; x < 2; ++x) {
        auto u = [&uu, x, y](int dx, int dy) {
          auto q = uu[(y + dy) * 3 + (x + dx)];
          return IsNan(q) ? T(0) : q;
        };
        res[i++] = u(0, 0) + u(1, 0) + u(0, 1) + u(1, 1);
      }
    }
    return res;
  }
  // TODO: move to an util file
  // Computes gradients in nodes.
  // uu: volume fractions in cells on a 3x3 stencil
  static std::array<Vect, 4> GradientNodes(const std::array<Scal, 9>& uu) {
    std::array<Vect, 4> res;
    size_t i = 0;
    for (int y = 0; y < 2; ++y) {
      for (int x = 0; x < 2; ++x) {
        auto u = [&uu, x, y](int dx, int dy) {
          return uu[(y + dy) * 3 + (x + dx)];
        };
        auto& g = res[i++];
        g[0] = ((u(1, 0) + u(1, 1)) - (u(0, 0) + u(0, 1))) * 0.5;
        g[1] = ((u(0, 1) + u(1, 1)) - (u(0, 0) + u(1, 0))) * 0.5;
      }
    }
    return res;
  }
  // Computes gradients in nodes.
  // uu: volume fractions in cells on a 3x3x3 stencil
  static std::array<Vect, 8> GradientNodes(const std::array<Scal, 27>& uu) {
    std::array<Vect, 8> res;
    size_t i = 0;
    for (int z = 0; z < 2; ++z) {
      for (int y = 0; y < 2; ++y) {
        for (int x = 0; x < 2; ++x) {
          auto u = [&uu, x, y, z](int dx, int dy, int dz) {
            return uu[(z + dz) * 9 + (y + dy) * 3 + (x + dx)];
          };
          auto& g = res[i++];
          g[0] = ((u(1, 0, 0) + u(1, 1, 0) + u(1, 0, 1) + u(1, 1, 1)) -
                  (u(0, 0, 0) + u(0, 1, 0) + u(0, 0, 1) + u(0, 1, 1))) *
                 0.25;
          g[1] = ((u(0, 1, 0) + u(1, 1, 0) + u(0, 1, 1) + u(1, 1, 1)) -
                  (u(0, 0, 0) + u(1, 0, 0) + u(0, 0, 1) + u(1, 0, 1))) *
                 0.25;
          g[2] = ((u(0, 0, 1) + u(1, 0, 1) + u(0, 1, 1) + u(1, 1, 1)) -
                  (u(0, 0, 0) + u(1, 0, 0) + u(0, 1, 0) + u(1, 1, 0))) *
                 0.25;
        }
      }
    }
    return res;
  }
  // Computes gradients in nodes.
  // uu: volume fractions in cells on a 3x3x3x3 stencil
  static std::array<Vect, 16> GradientNodes(const std::array<Scal, 81>& uu) {
    std::array<Vect, 16> res;
    size_t i = 0;
    for (int w = 0; w < 2; ++w) {
      for (int z = 0; z < 2; ++z) {
        for (int y = 0; y < 2; ++y) {
          for (int x = 0; x < 2; ++x) {
            auto u = [&uu, x, y, z, w](int dx, int dy, int dz, int dw) {
              return uu[(w + dw) * 27 + (z + dz) * 9 + (y + dy) * 3 + (x + dx)];
            };
            auto& g = res[i++];
            g[0] = ((u(1, 0, 0, 0) + u(1, 1, 0, 0) + u(1, 0, 1, 0) +
                     u(1, 1, 1, 0) +
                     (u(1, 0, 0, 1) + u(1, 1, 0, 1) + u(1, 0, 1, 1) +
                      u(1, 1, 1, 1)) -
                     (u(0, 0, 0, 0) + u(0, 1, 0, 0) + u(0, 0, 1, 0) +
                      u(0, 1, 1, 0)) +
                     u(0, 0, 0, 1) + u(0, 1, 0, 1) + u(0, 0, 1, 1) +
                     u(0, 1, 1, 1))) *
                   0.125;
            g[1] = ((u(0, 1, 0, 0) + u(1, 1, 0, 0) + u(0, 1, 1, 0) +
                     u(1, 1, 1, 0) + u(0, 1, 0, 1) + u(1, 1, 0, 1) +
                     u(0, 1, 1, 1) + u(1, 1, 1, 1)) -
                    (u(0, 0, 0, 0) + u(1, 0, 0, 0) + u(0, 0, 1, 0) +
                     u(1, 0, 1, 0) + u(0, 0, 0, 1) + u(1, 0, 0, 1) +
                     u(0, 0, 1, 1) + u(1, 0, 1, 1))) *
                   0.125;
            g[2] = ((u(0, 0, 1, 0) + u(1, 0, 1, 0) + u(0, 1, 1, 0) +
                     u(1, 1, 1, 0) + u(0, 0, 1, 1) + u(1, 0, 1, 1) +
                     u(0, 1, 1, 1) + u(1, 1, 1, 1)) -
                    (u(0, 0, 0, 0) + u(1, 0, 0, 0) + u(0, 1, 0, 0) +
                     u(1, 1, 0, 0) + u(0, 0, 0, 1) + u(1, 0, 0, 1) +
                     u(0, 1, 0, 1) + u(1, 1, 0, 1))) *
                   0.125;
          }
        }
      }
    }
    return res;
  }
  // bcfill: if >=0. add triangles from
  // fcus [a]: sum of volume fractions, add triangles from SuCells if not null
  static void DumpPolyMarch(
      const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<const FieldCell<Scal>*>& fccl,
      const Multi<const FieldCell<Vect>*>& fcn,
      const Multi<const FieldCell<Scal>*>& fca,
      const Multi<const FieldCell<bool>*>& fci, std::string fn, Scal t,
      bool bin, bool merge, Scal iso, const FieldCell<Scal>* fcus, M& m) {
    (void)fca;
    (void)fci;
    auto sem = m.GetSem("dumppolymarch");
    struct {
      std::vector<std::vector<Vect>> dl; // polygons
      std::vector<std::vector<Vect>> dln; // normal
      std::vector<Scal> dlc; // cells index
      std::vector<Scal> dll; // layer
      std::vector<Scal> dlcl; // color
    } * ctx(sem);
    auto& dl = ctx->dl;
    auto& dln = ctx->dln;
    auto& dlc = ctx->dlc;
    auto& dll = ctx->dll;
    auto& dlcl = ctx->dlcl;
    if (sem("local")) {
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
            for (auto cn : m.Stencil(c)) {
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
                // volume fraction in cells on 3x3x3 stencil
                const auto uu =
                    GetStencilValues<Scal>(layers, fcu, fccl, cn, cl, m);
                // volume fraction in nodes on 2x2x2 stencil
                const auto uun = InterpolateToNodes(uu);
                const auto nnc =
                    GetStencilValues<Vect>(layers, fcn, fccl, cn, cl, m);
                const auto nnn = SumToNodesNan(nnc);
                std::vector<std::vector<Vect>> vv, vvn;
                GetMarchTriangles(uun, nnn, m.GetCenter(cn), h, iso, vv, vvn);
                for (size_t j = 0; j < vv.size(); ++j) {
                  dl.push_back(vv[j]);
                  dln.push_back(vvn[j]);
                  dlc.push_back(m.GetHash(cn));
                  dll.push_back(i);
                  dlcl.push_back(cl);
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
          const MIdx w = bc.GetMIdx(c);
          if (!(MIdx(1) <= w && w < m.GetGlobalSize() - MIdx(1))) {
            auto uu = GetStencilValues<Scal>(*fcus, c, m);
            auto uun = InterpolateToNodes(uu);
            auto nn = GradientNodes(uu);
            for (auto& n : nn) {
              n = -n;
            }
            std::vector<std::vector<Vect>> vv, vvn;
            GetMarchTriangles(uun, nn, m.GetCenter(c), h, iso, vv, vvn);
            for (size_t j = 0; j < vv.size(); ++j) {
              dl.push_back(vv[j]);
              dln.push_back(vvn[j]);
              dlc.push_back(m.GetHash(c));
              dll.push_back(-1);
              dlcl.push_back(-1);
            }
          }
        }
      }
      m.Reduce(&dl, Reduction::concat);
      m.Reduce(&dln, Reduction::concat);
      m.Reduce(&dlc, Reduction::concat);
      m.Reduce(&dll, Reduction::concat);
      m.Reduce(&dlcl, Reduction::concat);
    }
    if (sem("write")) {
      if (m.IsRoot()) {
        std::cerr << std::fixed << std::setprecision(8) << "dump"
                  << " t=" << t << " to " << fn << std::endl;
        WriteVtkPoly<Vect>(
            fn, dl, &dln, {&dlc, &dll, &dlcl}, {"c", "l", "cl"},
            "Interface from marching cubes", true, bin, merge);
      }
    }
  }

  // Initializes usermap.
  // fccl0: known colors
  // fccl: colors to reduce
  static void UserMap(
      const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fccl0,
      const Multi<const FieldCell<Scal>*>& fccl, std::map<Scal, Scal>& usermap,
      M& m) {
    auto sem = m.GetSem("usermap");
    struct {
      std::vector<Scal> vcl, vcln;
    } * ctx(sem);
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
      m.Reduce(&vcl, Reduction::concat);
      m.Reduce(&vcln, Reduction::concat);
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
  static void Grid(
      const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fccl,
      const Multi<FieldCell<Scal>*>& fcclt, M& m) {
    auto sem = m.GetSem("grid");
    struct {
      std::vector<Scal> merge0; // colors from corners
      std::vector<Scal> merge1; // colors from parents
    } * ctx(sem);
    auto& merge0 = ctx->merge0;
    auto& merge1 = ctx->merge1;
    if (sem("local")) {
      // Collect neighbor colors in corners
      const IdxCell c = *m.Cells().begin();
      for (auto l : layers) {
        if ((*fccl[l])[c] != kClNone) {
          for (size_t q : {0, 1, 2}) {
            IdxCell cm = m.GetCell(c, IdxNci(q));
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
      m.Reduce(&merge0, Reduction::concat);
      m.Reduce(&merge1, Reduction::concat);
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
      m.Bcast(&merge0);
      m.Bcast(&merge1);
    }
    if (sem("apply")) {
      std::map<Scal, Scal> map;
      for (size_t i = 0; i < merge0.size(); ++i) {
        map[merge0[i]] = merge1[i];
      }
      for (auto f : m.Faces()) { // FIXME: inner cells traversed twice
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
  static void ReduceColor(
      const GRange<size_t>& layers, const Multi<FieldCell<Scal>*>& fccl,
      const std::map<Scal, Scal>& usermap, Scal clfixed, M& m) {
    auto sem = m.GetSem("recolor");
    struct {
      std::vector<std::vector<Scal>> vvcl; // all colors
      std::vector<Scal> vcl, vcln;
    } * ctx(sem);
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
      m.Reduce(&vvcl, Reduction::concat);
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
              // map[cl] = usermap.count(cl) ? usermap.at(cl) : -2;
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

  static void Init(
      const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<FieldCell<Scal>*>& fccl, Multi<FieldCell<Scal>>& fcclt,
      Scal clfixed, Vect clfixed_x, Scal coalth, M& m) {
    auto sem = m.GetSem("recolor_init");
    struct {
      std::pair<typename M::Scal, int> cldist; // color,mesh_id
    } * ctx(sem);
    auto& cldist = ctx->cldist;
    if (sem("clfixed")) {
      // block nearest to clfixed_x
      if (clfixed >= 0) {
        IdxCell c = m.FindNearestCell(clfixed_x);
        cldist.first = m.GetCenter(c).dist(clfixed_x);
        cldist.second = m.GetId();
        m.Reduce(&cldist, Reduction::minloc);
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

  static void RecolorDirect(
      const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<FieldCell<Scal>*>& fccl,
      const Multi<const FieldCell<Scal>*>& fccl0, Scal clfixed, Vect clfixed_x,
      Scal coalth, const MapEmbed<BCond<Scal>>& mfc_cl, bool verb, bool reduce,
      bool grid, M& m) {
    auto sem = m.GetSem("recolor");
    struct {
      std::map<Scal, Scal> usermap;
      Scal tries;
      Multi<FieldCell<Scal>> fcclt; // tmp color
    } * ctx(sem);
    auto& fcclt = ctx->fcclt;
    if (sem.Nested()) {
      Init(layers, fcu, fccl, fcclt, clfixed, clfixed_x, coalth, m);
    }
    if (sem("reflect")) {
      for (auto i : layers) {
        BcApply(fcclt[i], mfc_cl, m);
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
      static constexpr size_t sw = 1; // stencil half-width
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
      m.Reduce(&ctx->tries, Reduction::max);
    }
    if (sem("reflect")) {
      for (auto i : layers) {
        BcApply(fcclt[i], mfc_cl, m);
      }
    }
    if (sem("check")) {
      if (verb && m.IsRoot()) {
        std::cerr << "recolor:"
                  << " max tries: " << ctx->tries << std::endl;
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

  static void RecolorUnionFind(
      const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<FieldCell<Scal>*>& fccl,
      const Multi<const FieldCell<Scal>*>& fccl0, Scal clfixed, Vect clfixed_x,
      Scal coalth, const MapEmbed<BCond<Scal>>& mfc_cl, bool verb, bool reduce,
      bool grid, M& m) {
    auto sem = m.GetSem("recolor");
    struct {
      std::map<Scal, Scal> usermap;
      Scal tries;
      Multi<FieldCell<Scal>> fcclt; // tmp color
      Multi<FieldCell<IdxCell>> fcc; // root cell
      Multi<FieldCell<char>> fcl; // root layer
    } * ctx(sem);
    auto& fcclt = ctx->fcclt;
    auto& fcc = ctx->fcc;
    auto& fcl = ctx->fcl;
    if (sem.Nested()) {
      Init(layers, fcu, fccl, fcclt, clfixed, clfixed_x, coalth, m);
    }
    if (sem("reflect")) {
      for (auto i : layers) {
        BcApply(fcclt[i], mfc_cl, m);
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
      auto Clt = [&](Pair p) -> Scal& { return fcclt[p.second][p.first]; };
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
      static constexpr size_t sw = 1; // stencil half-width
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
      m.Reduce(&ctx->tries, Reduction::max);
    }
    if (sem("reflect")) {
      for (auto i : layers) {
        BcApply(fcclt[i], mfc_cl, m);
      }
    }
    if (sem("check")) {
      if (verb && m.IsRoot()) {
        std::cerr << "recolor:"
                  << " max tries: " << ctx->tries << std::endl;
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

  static void Recolor(
      const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<FieldCell<Scal>*>& fccl,
      const Multi<const FieldCell<Scal>*>& fccl0, Scal clfixed, Vect clfixed_x,
      Scal coalth, const MapEmbed<BCond<Scal>>& mfc, bool verb, bool unionfind,
      bool reduce, bool grid, M& m) {
    if (unionfind) {
      return RecolorUnionFind(
          layers, fcu, fccl, fccl0, clfixed, clfixed_x, coalth, mfc, verb,
          reduce, grid, m);
    }
    return RecolorDirect(
        layers, fcu, fccl, fccl0, clfixed, clfixed_x, coalth, mfc, verb, reduce,
        grid, m);
  }
};

template <class M_>
UVof<M_>::UVof() : imp(new Imp()) {}

template <class M_>
UVof<M_>::~UVof() = default;

template <class M_>
void UVof<M_>::DumpPoly(
    const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fcu,
    const Multi<const FieldCell<Scal>*>& fccl,
    const Multi<const FieldCell<Vect>*>& fcn,
    const Multi<const FieldCell<Scal>*>& fca,
    const Multi<const FieldCell<bool>*>& fci, std::string fn, Scal t, bool bin,
    bool merge, M& m) {
  imp->DumpPoly(layers, fcu, fccl, fcn, fca, fci, fn, t, bin, merge, m);
}

template <class M_>
void UVof<M_>::DumpPoly(
    const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcn,
    const FieldCell<Scal>& fca, const FieldCell<bool>& fci, std::string fn,
    Scal t, bool bin, bool merge, M& m) {
  GRange<size_t> layers(0, 1);
  const FieldCell<Scal>* fccl(nullptr);
  imp->DumpPoly(layers, &fcu, fccl, &fcn, &fca, &fci, fn, t, bin, merge, m);
}

template <class M_>
void UVof<M_>::DumpPolyMarch(
    const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fcu,
    const Multi<const FieldCell<Scal>*>& fccl,
    const Multi<const FieldCell<Vect>*>& fcn,
    const Multi<const FieldCell<Scal>*>& fca,
    const Multi<const FieldCell<bool>*>& fci, std::string fn, Scal t, bool bin,
    bool merge, Scal iso, const FieldCell<Scal>* fcus, M& m) {
  imp->DumpPolyMarch(
      layers, fcu, fccl, fcn, fca, fci, fn, t, bin, merge, iso, fcus, m);
}

template <class M_>
void UVof<M_>::Recolor(
    const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fcu,
    const Multi<FieldCell<Scal>*>& fccl,
    const Multi<const FieldCell<Scal>*>& fccl0, Scal clfixed, Vect clfixed_x,
    Scal coalth, const MapEmbed<BCond<Scal>>& mfcu, bool verb, bool unionfind,
    bool reduce, bool grid, M& m) {
  Imp::Recolor(
      layers, fcu, fccl, fccl0, clfixed, clfixed_x, coalth, mfcu, verb,
      unionfind, reduce, grid, m);
}

template <class M_>
auto UVof<M_>::GetAdvectionBc(
    const M& m, const MapEmbed<BCondAdvection<Scal>>& mfc)
    -> std::tuple<
        MapEmbed<BCond<Scal>>, MapEmbed<BCond<Scal>>, MapEmbed<BCond<Scal>>,
        MapEmbed<BCond<Vect>>, MapEmbed<BCond<Scal>>> {
  using MIdx = typename M::MIdx;
  using TRM = Trackerm<M>;

  MapEmbed<BCond<Scal>> me_vf; // volume fraction
  MapEmbed<BCond<Scal>> me_cl; // color
  MapEmbed<BCond<Scal>> me_im; // image
  MapEmbed<BCond<Vect>> me_n; // normal
  MapEmbed<BCond<Scal>> me_a; // plane constant
  using Halo = typename BCondAdvection<Scal>::Halo;
  for (auto& p : mfc.GetMapFace()) {
    const IdxFace f = p.first;
    const auto& bc = p.second;
    const auto nci = bc.nci;
    switch (bc.halo) {
      case Halo::reflect:
        me_vf[f] = BCond<Scal>(BCondType::reflect, nci);
        me_cl[f] = BCond<Scal>(BCondType::reflect, nci);
        me_im[f] = BCond<Scal>(BCondType::reflect, nci);
        me_n[f] = BCond<Vect>(BCondType::reflect, nci);
        me_a[f] = BCond<Scal>(BCondType::reflect, nci);
        break;
      case Halo::fill:
        me_vf[f] = BCond<Scal>(BCondType::dirichlet, nci, bc.fill_vf);
        me_cl[f] = BCond<Scal>(BCondType::dirichlet, nci, bc.fill_cl);
        MIdx wim(0);
        wim[size_t(m.GetDir(f))] = (nci == 1 ? -1 : 1);
        me_im[f] = BCond<Scal>(BCondType::dirichlet, nci, TRM::Pack(wim));
        me_n[f] = BCond<Vect>(BCondType::dirichlet, nci, GetNan<Vect>());
        me_a[f] = BCond<Scal>(BCondType::dirichlet, nci, GetNan<Scal>());
        break;
    }
  }
  return {me_vf, me_cl, me_im, me_n, me_a};
}

template <class M>
constexpr typename M::Scal UVof<M>::Imp::kClNone;
