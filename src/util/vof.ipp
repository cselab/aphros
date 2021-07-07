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

  static void DumpPoly(const DumpPolyArgs& args, M& m) {
    auto sem = m.GetSem("dumppoly");
    struct {
      std::vector<std::vector<Vect>> polygons;
      std::vector<Scal> d_cell;
      std::map<std::string, std::vector<Scal>> data;
    } * ctx(sem);
    auto& t = *ctx;

    auto& fcu = args.fcu;
    auto& fccl = args.fccl;
    auto& fcn = args.fcn;
    auto& fca = args.fca;
    auto& fci = args.fci;
    auto& layers = args.layers;

    if (sem("local")) {
      fcu.assert_size(layers);
      fccl.assert_size(layers);
      fcn.assert_size(layers);
      fca.assert_size(layers);
      fci.assert_size(layers);
      fassert(args.filename.length());

      auto for_each_cell = [&](auto body) {
        for (auto l : args.layers) {
          for (auto c : m.Cells()) {
            if (!IsNan((*fcu[l])[c]) && !IsNan((*fcn[l])[c]) &&
                !IsNan((*fca[l])[c]) && (*fci[l])[c]) {
              body(l, c);
            }
          }
        }
      };

      for_each_cell([&](size_t l, IdxCell c) { //
        t.polygons.push_back(R::GetCutPoly(
            m.GetCenter(c), (*fcn[l])[c], (*fca[l])[c], m.GetCellSize()));
      });

      for_each_cell([&](size_t, IdxCell c) { //
        t.d_cell.push_back(m.GetHash(c));
      });

      if (args.dump_layer) {
        auto& d = t.data["l"];
        for_each_cell([&](size_t l, IdxCell) { //
          d.push_back(l);
        });
      }
      if (args.dump_color) {
        auto& d = t.data["cl"];
        for_each_cell([&](size_t l, IdxCell c) { //
          d.push_back(fccl[l] ? (*fccl[l])[c] : 0);
        });
      }

      fassert_equal(args.extra_fields.size(), args.extra_names.size());

      for (size_t i = 0; i < args.extra_fields.size(); ++i) {
        auto& d = t.data[args.extra_names[i]];
        for_each_cell([&](size_t l, IdxCell c) { //
          d.push_back((*args.extra_fields[i][l])[c]);
        });
      }

      m.Reduce(&t.polygons, Reduction::concat);
      m.Reduce(&t.d_cell, Reduction::concat);
      for (auto& p : t.data) {
        m.Reduce(&p.second, Reduction::concat);
      }
    }
    if (sem("write")) {
      if (m.IsRoot()) {
        std::vector<size_t> idx(t.d_cell.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::stable_sort(idx.begin(), idx.end(), [&](size_t i0, size_t i1) {
          return t.d_cell[i0] < t.d_cell[i1];
        });

        // Sort output by the global cell index to make the output independent
        // of the domain partitioning
        Reorder(t.polygons, idx);
        Reorder(t.d_cell, idx);
        for (auto& p : t.data) {
          Reorder(p.second, idx);
        }

        if (args.verbose) {
          std::cerr << std::fixed << std::setprecision(8) << "dump"
                    << " t=" << args.time << " to " << args.filename
                    << std::endl;
        }

        std::vector<const std::vector<Scal>*> fields;
        std::vector<std::string> names;

        for (auto& p : t.data) {
          names.push_back(p.first);
          fields.push_back(&p.second);
        }

        WriteVtkPoly<Vect>(
            args.filename, t.polygons, nullptr, fields, names,
            "Interface from PLIC", true, args.binary, args.merge);
      }
    }
  }

  // Constructs iso-surface triangles with marching cubes
  // uu: volume fraction in cell nodes
  // nn: normals in cell nodes
  // xc: cell center
  // h: cell size
  // iso: isovalue for surface uu=iso
  //
  // Output:
  // vv: triangles as arrays of vertices
  // vvn: normals for every vertex
  //
  // 3D
  static void GetMarchTriangles(
      const std::array<Scal, 8>& uu, const std::array<Vect3, 8>& nn,
      const Vect3& xc, const Vect3& h, Scal iso, /*out*/
      std::vector<std::vector<Vect3>>& vv,
      std::vector<std::vector<Vect3>>& vvn) {
    (void)MARCH_O[0][0]; // suppress unused variable warning from march.h
    std::array<double, 8> uuz = uu;
    for (auto& u : uuz) {
      u -= iso;
    }
    int nt; // output number of triangles
    constexpr int kMaxNt = MARCH_NTRI;
    std::array<double, kMaxNt> tri; // triangles as flat array of coordinates
    // Each vertex belongs to a cell edge, the following describes which one.
    std::array<int, kMaxNt> vc0; // index in `uu` of first endpoint of edge
    std::array<int, kMaxNt> vc1; // index in `uu` of second endpoint of edge
    std::array<double, kMaxNt> vw; // position of point on edge
                                   // 0: first endpoint, 1: second endpoint
    march_cube_location(
        uuz.data(), &nt, tri.data(), vc0.data(), vc1.data(), vw.data());
    assert(size_t(nt) * 3 * 3 <= tri.size());

    {
      vv.clear();
      vv.resize(nt, std::vector<Vect3>(3));
      size_t i = 0;
      for (auto& v : vv) {
        for (auto& x : v) {
          x[0] = tri[i++] - 0.5;
          x[1] = tri[i++] - 0.5;
          x[2] = tri[i++] - 0.5;
          x = xc + h * x;
        }
      }
    }
    {
      vvn.clear();
      vvn.resize(nt, std::vector<Vect3>(3));
      size_t i = 0;
      for (auto& vn : vvn) {
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
  // Constructs iso-surface lines with marching cubes
  static void GetMarchTriangles(
      const std::array<Scal, 4>& uu, const std::array<Vect2, 4>& nn,
      const Vect2& xc, const Vect2& h, Scal iso,
      /*out*/ std::vector<std::vector<Vect2>>& vv,
      std::vector<std::vector<Vect2>>& vvn) {
    std::array<double, 8> uuz;
    // Duplicate 2D volume fractions in z.
    for (size_t i = 0; i < 8; ++i) {
      uuz[i] = uu[i % 4] - iso;
    }
    int nt; // output number of triangles
    constexpr int kMaxNt = MARCH_NTRI;
    std::array<double, kMaxNt> tri; // triangles as flat array of coordinates
    // Each vertex belongs to a cell edge, the following describes which one.
    std::array<int, kMaxNt> vc0; // index in `uu` of first endpoint of edge
    std::array<int, kMaxNt> vc1; // index in `uu` of second endpoint of edge
    std::array<double, kMaxNt> vw; // position of point on edge
                                   // 0: first endpoint, 1: second endpoint
    march_cube_location(
        uuz.data(), &nt, tri.data(), vc0.data(), vc1.data(), vw.data());
    assert(size_t(nt) * 3 * 3 <= tri.size());

    // 3D triangles are stored in tri, vc0, vc1, vw.
    // Extract edges with z=0.
    vv.clear();
    vvn.clear();
    for (int t = 0; t < nt; ++t) {
      std::vector<Vect2> v;
      std::vector<Vect2> vn;
      for (int j = 0; j < 3; ++j) {
        const int i = t * 3 + j; // vertex index
        if (tri[3 * i + 2] == 0) {
          // vertex
          const Vect2 x(tri[3 * i] - 0.5, tri[3 * i + 1] - 0.5);
          v.push_back(xc + h * x);
          // normal
          const Scal w = vw[i];
          const int c0 = vc0[i];
          const int c1 = vc1[i];
          const Vect2 n = nn[c0 % 4] * (1 - w) + nn[c1 % 4] * w;
          vn.push_back(n);
        }
      }
      if (v.size() == 2) {
        vv.push_back(v);
        vvn.push_back(vn);
      }
    }
  }
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
  static void DumpPolyMarch(
      const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<const FieldCell<Scal>*>& fccl,
      const Multi<const FieldCell<Vect>*>& fcn, std::string filename, Scal time,
      bool poly, bool bin, bool merge, Scal iso, const FieldCell<Scal>* fcus,
      M& m) {
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
                  << " t=" << time << " to " << filename << std::endl;
        WriteVtkPoly<Vect>(
            filename, dl, &dln, {&dlc, &dll, &dlcl}, {"c", "l", "cl"},
            "Interface from marching cubes", poly, bin, merge);
      }
    }
  }

  // Converts a pair of fields to map, reduced over all blocks.
  // Output:
  // `map`: resulting on root block, and empty on other blocks
  static void PairToMap(
      const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fc_key,
      const Multi<const FieldCell<Scal>*>& fc_value, std::map<Scal, Scal>& map,
      M& m) {
    auto sem = m.GetSem("pair_to_map");
    struct {
      std::vector<Scal> keys, values;
    } * ctx(sem);
    auto& t = *ctx;
    if (sem("local")) {
      map.clear();
      for (auto c : m.AllCells()) {
        for (auto l : layers) {
          const Scal key = (*fc_key[l])[c];
          const Scal value = (*fc_value[l])[c];
          if (key != kClNone && value != kClNone) {
            if (!map.count(key)) {
              map[key] = value;
            }
          }
        }
      }
      for (auto p : map) {
        t.keys.push_back(p.first);
        t.values.push_back(p.second);
      }
      m.Reduce(&t.keys, Reduction::concat);
      m.Reduce(&t.values, Reduction::concat);
    }
    if (sem("gather")) {
      map.clear();
      if (m.IsRoot()) {
        for (size_t i = 0; i < t.keys.size(); ++i) {
          map[t.keys[i]] = t.values[i];
        }
      }
    }
  }

  // Connected component labeling over sparse set of cells
  // including one cell from each block.
  static void RecolorSparse(
      const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fccl,
      const Multi<FieldCell<Scal>*>& fccl_new, M& m) {
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
                Scal cl = (*fccl_new[l])[c];
                Scal clm = (*fccl_new[lm])[cm];
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
        int iters = 0;
        // Find minimal color connected through pairs
        while (true) {
          bool changed = false;
          for (auto& p : map) {
            if (map.count(p.second)) {
              p.second = map[p.second];
              changed = true;
            }
          }
          if (!changed) {
            break;
          }
          ++iters;
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
            auto& cl = (*fccl_new[l])[c];
            if (map.count(cl)) {
              cl = map[cl];
            }
          }
        }
      }
    }
  }

  // Reduces the color space. Applies maximum subset of `stable_map` that
  // produces unique colors and replaces other colors with smallest integers.
  // If `clfixed >= 0`, keeps color `clfixed` unmodified.
  static void ReduceColor(
      const GRange<size_t>& layers, const Multi<FieldCell<Scal>*>& fccl,
      const std::map<Scal, Scal>& stable_map, Scal clfixed, M& m) {
    auto sem = m.GetSem("recolor");
    struct {
      std::vector<std::vector<Scal>> vvcl; // colors from all blocks
      std::vector<Scal> vcl, vcl_new; // colors on current block
    } * ctx(sem);
    auto& t = *ctx;
    if (sem("gather")) {
      // Gather all colors in domain
      std::set<Scal> set;
      for (auto l : layers) {
        for (auto c : m.AllCells()) {
          auto& cl = (*fccl[l])[c];
          if (cl != kClNone) {
            set.insert(cl);
          }
        }
      }
      t.vcl = std::vector<Scal>(set.begin(), set.end());
      t.vvcl = {t.vcl};
      m.Reduce(&t.vvcl, Reduction::concat);
    }
    if (sem("reduce")) {
      // Replace with reduced set, applying stable_map if possible
      if (m.IsRoot()) {
        std::map<Scal, Scal> map;
        std::set<Scal> used;
        auto Add = [&](Scal cl, Scal a) { // add color pair to map
          map[cl] = a;
          used.insert(a);
        };
        Add(kClNone, kClNone);
        if (clfixed >= 0) {
          Add(clfixed, clfixed);
        }
        Scal cl_next = 0; // next unused color
        for (auto p : stable_map) {
          cl_next = std::max(cl_next, p.second);
        }
        cl_next += 1;
        for (auto& vcl : t.vvcl) {
          for (auto& cl : vcl) {
            if (!map.count(cl)) {
              if (!stable_map.count(cl) || used.count(stable_map.at(cl))) {
                Add(cl, cl_next);
                cl_next += 1;
              } else { // stable_map.count(cl) && !used.count(stable_map[cl])
                Add(cl, stable_map.at(cl));
              }
            }
            cl = map[cl];
          }
        }
        m.Scatter({&t.vvcl, &t.vcl_new});
      } else {
        m.Scatter({nullptr, &t.vcl_new});
      }
    }
    if (sem("apply")) {
      // Apply the new set from vcl_new
      std::map<Scal, Scal> map;
      for (size_t i = 0; i < t.vcl.size(); ++i) {
        map[t.vcl[i]] = t.vcl_new[i];
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

  // Initializes color field with unique colors over all cells and layers.
  // fcu: volume fractions
  // fccl: old colors
  // fccl_new: new colors
  // clfixed, clfixed_x: color and position of cell with fixed color
  // coalth: threshold for total volume fraction
  static void InitUniqueColors(
      const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<FieldCell<Scal>*>& fccl, Multi<FieldCell<Scal>>& fccl_new,
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
      fccl_new.Reinit(layers, m, kClNone);
      // initial unique color
      Scal q = m.GetId() * m.GetInBlockCells().size() * layers.size();
      for (auto i : layers) {
        for (auto c : m.Cells()) {
          if ((*fccl[i])[c] != kClNone) {
            fccl_new[i][c] = (q += 1);
          }
        }
        if (cldist.second == m.GetId()) {
          IdxCell c = m.FindNearestCell(clfixed_x);
          if ((*fccl[i])[c] != kClNone) {
            fccl_new[i][c] = clfixed;
          }
        }
        m.Comm(&fccl_new[i]);
      }

      // Merge colors in cells with total volume fraction above threshold.
      for (auto i : layers) {
        for (auto c : m.Cells()) {
          if ((*fccl[i])[c] != kClNone) {
            for (auto j : layers) {
              if (j != i && (*fccl[j])[c] != kClNone) {
                if ((*fcu[i])[c] + (*fcu[j])[c] > coalth) {
                  Scal cl = std::min(fccl_new[i][c], fccl_new[j][c]);
                  fccl_new[i][c] = cl;
                  fccl_new[j][c] = cl;
                }
              }
            }
          }
        }
      }
    }
  }

  // Connected component labeling with iterative propagation
  // of colors over stencil.
  static void RecolorDirect(
      const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<FieldCell<Scal>*>& fccl,
      const Multi<const FieldCell<Scal>*>& fccl_stable, Scal clfixed,
      Vect clfixed_x, Scal coalth, const MapEmbed<BCond<Scal>>& mfc_cl,
      bool verb, bool reduce, bool grid, M& m) {
    auto sem = m.GetSem("recolor");
    struct {
      std::map<Scal, Scal> stable_map;
      Scal iters;
      Multi<FieldCell<Scal>> fccl_new;
    } * ctx(sem);
    auto& t = *ctx;
    auto& fccl_new = ctx->fccl_new;
    if (sem.Nested()) {
      // Initialize `fccl_new` with unique colors over all cells and layers.
      InitUniqueColors(
          layers, fcu, fccl, fccl_new, clfixed, clfixed_x, coalth, m);
    }
    if (sem("reflect")) {
      // Fill colors in halo cells according to boundary conditions.
      for (auto l : layers) { //
        BcApply(fccl_new[l], mfc_cl, m);
      }
    }
    sem.LoopBegin();
    if (grid && sem.Nested()) {
      // Run labeling over a sparse set of points to speedup propagation
      // in large components.
      RecolorSparse(layers, fccl, fccl_new, m);
    }
    if (sem("min")) {
      size_t iters = 0;
      while (true) {
        bool changed = false;
        for (auto l : layers) {
          for (auto c : m.Cells()) {
            if ((*fccl[l])[c] != kClNone) {
              // Update colors with minimum over neighbours
              for (auto cn : m.Stencil(c)) {
                for (auto ln : layers) {
                  if ((*fccl[ln])[cn] == (*fccl[l])[c]) {
                    if (fccl_new[ln][cn] < fccl_new[l][c]) {
                      changed = true;
                      fccl_new[l][c] = fccl_new[ln][cn];
                    }
                  }
                }
              }
            }
          }
        }
        if (!changed) {
          break;
        }
        ++iters;
      }
      for (auto l : layers) {
        m.Comm(&fccl_new[l]);
      }
      t.iters = iters;
      m.Reduce(&t.iters, Reduction::max);
    }
    if (sem("reflect")) {
      // Fill colors in halo cells according to boundary conditions.
      for (auto l : layers) {
        BcApply(fccl_new[l], mfc_cl, m);
      }
    }
    if (sem("check")) {
      if (verb && m.IsRoot()) {
        std::cerr << "recolor:"
                  << " max iters: " << t.iters << std::endl;
      }
      if (!t.iters) {
        sem.LoopBreak();
      }
    }
    sem.LoopEnd();
    if (reduce && sem.Nested()) {
      // Create map `stable_map  from new colors to `fccl_stable`
      PairToMap(layers, fccl_new, fccl_stable, t.stable_map, m);
    }
    if (reduce && sem.Nested()) {
      // Apply `stable_map` to new colors and reduce colors to smaller
      // intergers
      ReduceColor(layers, fccl_new, t.stable_map, clfixed, m);
    }
    if (sem("copy")) {
      for (auto c : m.AllCells()) {
        for (auto l : layers) {
          (*fccl[l])[c] = fccl_new[l][c];
        }
      }
    }
  }

  // Connected component labeling with union-find structure on local block
  static void RecolorUnionFind(
      const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<FieldCell<Scal>*>& fccl,
      const Multi<const FieldCell<Scal>*>& fccl_stable, Scal clfixed,
      Vect clfixed_x, Scal coalth, const MapEmbed<BCond<Scal>>& mfc_cl,
      bool verb, bool reduce, bool grid, M& m) {
    auto sem = m.GetSem("recolor");
    struct {
      std::map<Scal, Scal> stable_map;
      Scal iters;
      Multi<FieldCell<Scal>> fccl_new; // tmp color
      Multi<FieldCell<IdxCell>> fcc; // root cell
      Multi<FieldCell<char>> fcl; // root layer
    } * ctx(sem);
    auto& t = *ctx;
    auto& fccl_new = ctx->fccl_new;
    auto& fcc = ctx->fcc;
    auto& fcl = ctx->fcl;
    if (sem.Nested()) {
      InitUniqueColors(
          layers, fcu, fccl, fccl_new, clfixed, clfixed_x, coalth, m);
    }
    if (sem("reflect")) {
      for (auto i : layers) {
        BcApply(fccl_new[i], mfc_cl, m);
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
      RecolorSparse(layers, fccl, fccl_new, m);
    }
    if (sem("min")) {
      using Pair = std::pair<IdxCell, char>;
      size_t iters = 0;
      // Returns reference to color
      auto Clt = [&](Pair p) -> Scal& { return fccl_new[p.second][p.first]; };
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
          ++iters;
        }
      };

      // Merge all neighbors with the same color.
      for (auto c : m.Cells()) {
        for (auto l : layers) {
          if ((*fccl[l])[c] != kClNone) {
            for (auto cn : m.Stencil(c)) {
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
            ++iters;
          }
        }
      }
      // Update color from root
      for (auto c : m.SuCells()) {
        for (auto l : layers) {
          Pair p(c, l);
          if (Clt(Get(p)) < Clt(p)) {
            Clt(p) = Clt(Get(p));
            ++iters;
          }
        }
      }
      for (auto l : layers) {
        m.Comm(&fccl_new[l]);
      }
      t.iters = iters;
      m.Reduce(&t.iters, Reduction::max);
    }
    if (sem("reflect")) {
      for (auto i : layers) {
        BcApply(fccl_new[i], mfc_cl, m);
      }
    }
    if (sem("check")) {
      if (verb && m.IsRoot()) {
        std::cerr << "recolor:"
                  << " max iters: " << t.iters << std::endl;
      }
      if (!t.iters) {
        sem.LoopBreak();
      }
    }
    sem.LoopEnd();
    if (reduce && sem.Nested()) {
      PairToMap(layers, fccl_new, fccl_stable, t.stable_map, m);
    }
    if (reduce && sem.Nested()) {
      ReduceColor(layers, fccl_new, t.stable_map, clfixed, m);
    }
    if (sem("copy")) {
      for (auto c : m.AllCells()) {
        for (auto l : layers) {
          (*fccl[l])[c] = fccl_new[l][c];
        }
      }
    }
  }

  static void Recolor(
      const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<FieldCell<Scal>*>& fccl,
      const Multi<const FieldCell<Scal>*>& fccl_stable, Scal clfixed,
      Vect clfixed_x, Scal coalth, const MapEmbed<BCond<Scal>>& mfc, bool verb,
      bool unionfind, bool reduce, bool grid, M& m) {
    if (unionfind) {
      return RecolorUnionFind(
          layers, fcu, fccl, fccl_stable, clfixed, clfixed_x, coalth, mfc, verb,
          reduce, grid, m);
    }
    return RecolorDirect(
        layers, fcu, fccl, fccl_stable, clfixed, clfixed_x, coalth, mfc, verb,
        reduce, grid, m);
  }
};

template <class M_>
UVof<M_>::UVof() : imp(new Imp()) {}

template <class M_>
UVof<M_>::~UVof() = default;

template <class M_>
void UVof<M_>::DumpPoly(const DumpPolyArgs& args, M& m) {
  imp->DumpPoly(args, m);
}

template <class M_>
void UVof<M_>::DumpPolyMarch(
    const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fcu,
    const Multi<const FieldCell<Scal>*>& fccl,
    const Multi<const FieldCell<Vect>*>& fcn, std::string filename, Scal t,
    bool poly, bool bin, bool merge, Scal iso, const FieldCell<Scal>* fcus,
    M& m) {
  imp->DumpPolyMarch(
      layers, fcu, fccl, fcn, filename, t, poly, bin, merge, iso, fcus, m);
}

template <class M_>
void UVof<M_>::Recolor(
    const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fcu,
    const Multi<FieldCell<Scal>*>& fccl,
    const Multi<const FieldCell<Scal>*>& fccl_stable, Scal clfixed,
    Vect clfixed_x, Scal coalth, const MapEmbed<BCond<Scal>>& mfcu, bool verb,
    bool unionfind, bool reduce, bool grid, M& m) {
  Imp::Recolor(
      layers, fcu, fccl, fccl_stable, clfixed, clfixed_x, coalth, mfcu, verb,
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
