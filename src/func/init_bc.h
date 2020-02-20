// Created by Petr Karnakov on 20.02.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <cmath>
#include <functional>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "debug/isnan.h"
#include "geom/block.h"
#include "geom/field.h"
#include "geom/unique.h"
#include "geom/vect.h"
#include "parse/codeblocks.h"
#include "parse/vars.h"
#include "primlist.h"
#include "solver/embed.h"
#include "solver/fluid.h"

template <class M_>
struct UInitEmbedBc {
  using M = M_;
  using EB = Embed<M>;
  static constexpr size_t dim = M::dim;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  // Parses boundary conditions from file and writes to map.
  // filename: path to file
  // mec: output map
  static void Parse(
      std::string filename, const EB& eb,
      MapEmbed<UniquePtr<CondFaceFluid>>& mec) {
    auto& m = eb.GetMesh();
    using namespace fluid_condition;
    using Primitive = generic::Primitive<Scal>;
    std::ifstream fin(filename);
    if (!fin.good()) {
      throw std::runtime_error(
          "Can't open list of boundary conditions '" + filename + "'");
    }
    const std::vector<CodeBlock> bb = ParseCodeBlocks(fin);
    fin.close();
    std::vector<std::vector<Primitive>> vpp;
    for (auto& b : bb) {
      // b.name: descriptor of boundary condition
      // b.content: list of primitives
      std::stringstream s(b.content);
      vpp.push_back(UPrimList<Scal>::Parse(s, m.IsRoot(), m.GetEdim()));
    }

    for (auto& pp : vpp) {
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
      for (auto c : eb.SuCFaces()) {
        auto ls = lsmax(eb.GetFaceCenter(c));
        if (ls > 0) {
          mec[c] = UniquePtr<NoSlipWallFixed<M>>(Vect(0), 0);
        }
      }
    }
  }
};
