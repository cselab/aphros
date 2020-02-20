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

// desc: boundary condition descriptor
// nci: neighbour cell id, such that GetCell(f, nci) is an internal cell
template <class M>
UniquePtr<CondFaceFluid> ParseFluidFaceCond(std::string desc, size_t nci) {
  using namespace fluid_condition;
  using Vect = typename M::Vect;
  std::stringstream arg(desc);

  std::string name;
  arg >> name;

  if (name == "wall") {
    // wall <velocity>
    // No-slip wall.
    // zero derivative for pressure, fixed for velocity,
    // fill-conditions for volume fraction.
    Vect vel;
    arg >> vel;
    return UniquePtr<NoSlipWallFixed<M>>(vel, nci);
  } else if (name == "inlet") {
    // inlet <velocity>
    // Fixed velocity inlet.
    Vect vel;
    arg >> vel;
    return UniquePtr<InletFixed<M>>(vel, nci);
  } else if (name == "inletflux") {
    // inletflux <velocity> <id>
    // Fixed flux inlet. Flux defined by given velocity is redistributed
    // over all faces with same id.
    Vect vel;
    int id = 0;
    arg >> vel >> id;
    return UniquePtr<InletFlux<M>>(vel, id, nci);
  } else if (name == "outlet") {
    // Outlet. Velocity is extrapolated from neighbour cells and corrected
    // to yield zero total flux over outlet and inlet faces.
    return UniquePtr<OutletAuto<M>>(nci);
  } else if (name == "slipwall") {
    // Free-slip wall:
    // zero derivative for both pressure, velocity,
    // fill-conditions for volume fraction.
    // TODO: revise, should be non-penetration for velocity
    return UniquePtr<SlipWall<M>>(nci);
  } else if (name == "symm") {
    // Zero derivative for pressure, velocity and volume fraction
    // TODO: revise, should be non-penetration for velocity
    return UniquePtr<Symm<M>>(nci);
  }
  return nullptr;
}

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
    std::ifstream fin(filename);
    if (!fin.good()) {
      throw std::runtime_error(
          "Can't open list of boundary conditions '" + filename + "'");
    }
    const std::vector<CodeBlock> bb = ParseCodeBlocks(fin);
    fin.close();
    for (auto& b : bb) {
      // b.name: descriptor of boundary condition
      // b.content: list of primitives
      std::stringstream s(b.content);
      auto pp = UPrimList<Scal>::Parse(s, m.IsRoot(), m.GetEdim());
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
          mec[c] = ParseFluidFaceCond<M>(b.name, 0);
        }
      }
    }
  }
};
