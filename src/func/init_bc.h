// Created by Petr Karnakov on 20.02.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <cmath>
#include <functional>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <tuple>

#include "debug/isnan.h"
#include "dump/vtk.h"
#include "geom/block.h"
#include "geom/field.h"
#include "geom/vect.h"
#include "parse/codeblocks.h"
#include "parse/vars.h"
#include "primlist.h"
#include "solver/embed.h"
#include "solver/fluid.h"

// desc: discriptor
// nci: neighbour cell id, such that GetCell(f, nci) is an inner cell
template <class Vect>
BCondFluid<Vect> ParseBCondFluid(std::string desc, size_t nci) {
  std::stringstream arg(desc);

  std::string name;
  arg >> name;

  BCondFluid<Vect> bc;
  bc.nci = nci;

  if (name == "wall") {
    bc.type = BCondFluidType::wall;
    arg >> bc.velocity;
  } else if (name == "slipwall") {
    bc.type = BCondFluidType::slipwall;
    bc.velocity = Vect(0);
  } else if (name == "inlet") {
    bc.type = BCondFluidType::inlet;
    arg >> bc.velocity;
  } else if (name == "inletflux") {
    // inletflux <velocity> <id>
    bc.type = BCondFluidType::inletflux;
    arg >> bc.velocity >> bc.inletflux_id;
  } else if (name == "outlet") {
    bc.type = BCondFluidType::outlet;
  } else if (name == "symm") {
    bc.type = BCondFluidType::symm;
  } else {
    throw std::runtime_error(
        std::string() + __func__ + ": unknown name='" + name + "'");
  }
  return bc;
}

template <class M_>
struct UInitEmbedBc {
  using M = M_;
  using EB = Embed<M>;
  static constexpr size_t dim = M::dim;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  // Parses boundary conditions from stream and adds to map.
  // filename: path to file
  // Output:
  // mecond: output map
  // megroup: index of group of primitives
  // Returns descriptors of boundary conditions for each group.
  static std::tuple<
      MapEmbed<BCondFluid<Vect>>, MapEmbed<size_t>, std::vector<std::string>>
  Parse(std::istream& fin, const EB& eb) {
    auto& m = eb.GetMesh();
    const std::vector<CodeBlock> bb = ParseCodeBlocks(fin);
    auto is_boundary = [&m](IdxFace f, size_t& nci) -> bool {
      auto p = m.GetIndexFaces().GetMIdxDir(f);
      const size_t d(p.second);
      const auto w = p.first;
      if (d < m.GetEdim()) {
        if (w[d] == 0) {
          nci = 1;
          return true;
        } else if (w[d] == m.GetGlobalSize()[d]) {
          nci = 0;
          return true;
        }
      }
      return false;
    };
    std::vector<std::string> vdesc;
    MapEmbed<BCondFluid<Vect>> mebc;
    MapEmbed<size_t> megroup;
    for (size_t group = 0; group < bb.size(); ++group) {
      auto& b = bb[group];
      vdesc.push_back(b.name);
      // b.name: descriptor of boundary condition
      // b.content: list of primitives
      std::stringstream s(b.content);
      auto pp = UPrimList<Scal>::Parse(s, false, m.GetEdim());
      auto lsmax = [&pp](Vect x) {
        Scal lmax = -std::numeric_limits<Scal>::max(); // maximum level-set
        for (size_t i = 0; i < pp.size(); ++i) {
          const auto& p = pp[i];
          Scal li = p.ls(x);
          if (p.mod_minus) {
            li = -li;
          }
          if (p.mod_and) {
            lmax = std::min(lmax, li);
          } else {
            if (li > lmax) {
              lmax = li;
            }
          }
        }
        return lmax;
      };
      for (auto f : eb.SuFaces()) {
        size_t nci;
        if (is_boundary(f, nci)) {
          auto ls = lsmax(eb.GetFaceCenter(f));
          if (ls > 0) {
            mebc[f] = ParseBCondFluid<Vect>(b.name, nci);
            megroup[f] = group;
          }
        }
      }
      for (auto c : eb.SuCFaces()) {
        auto ls = lsmax(eb.GetFaceCenter(c));
        if (ls > 0) {
          mebc[c] = ParseBCondFluid<Vect>(b.name, 0);
          megroup[c] = group;
        }
      }
    }
    return {mebc, megroup, vdesc};
  }

  static void DumpPoly(
      const std::string filename, const MapEmbed<size_t>& megroup, const EB& eb,
      M& m) {
    auto sem = m.GetSem("dumppoly");
    struct {
      std::vector<std::vector<Vect>> dpoly;
      std::vector<Scal> dgroup;
    } * ctx(sem);
    auto& dpoly = ctx->dpoly;
    auto& dgroup = ctx->dgroup;
    if (sem("local")) {
      for (auto p : megroup.GetMapCell()) {
        dpoly.push_back(eb.GetCutPoly(p.first));
        dgroup.push_back(p.second);
      }
      for (auto p : megroup.GetMapFace()) {
        dpoly.push_back(eb.GetFacePoly(p.first));
        dgroup.push_back(p.second);
      }
      using TV = typename M::template OpCatVT<Vect>;
      m.Reduce(std::make_shared<TV>(&dpoly));
      using TS = typename M::template OpCatT<Scal>;
      m.Reduce(std::make_shared<TS>(&dgroup));
    }
    if (sem("write")) {
      if (m.IsRoot()) {
        WriteVtkPoly<Vect>(
            filename, dpoly, nullptr, {&dgroup}, {"group"},
            "Boundary conditions", true, true, true);
      }
    }
  }
};
