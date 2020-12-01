// Created by Petr Karnakov on 20.02.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <cmath>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>

#include "debug/isnan.h"
#include "dump/vtk.h"
#include "geom/block.h"
#include "geom/field.h"
#include "geom/vect.h"
#include "parse/codeblocks.h"
#include "parse/util.h"
#include "parse/vars.h"
#include "primlist.h"
#include "solver/embed.h"
#include "solver/fluid.h"

// desc: discriptor
// nci: neighbour cell id, such that GetCell(f, nci) is an inner cell
template <class Vect>
std::pair<bool, BCondFluid<Vect>> ParseBCondFluid(
    std::string desc, size_t nci, Vect face_center, Vect face_normal) {
  using Scal = typename Vect::value_type;
  std::stringstream arg(desc);

  std::string name;
  arg >> name;

  BCondFluid<Vect> bc;
  bc.nci = nci;

  if (name == "wall") {
    bc.type = BCondFluidType::wall;
    arg >> bc.velocity;
    bc.velocity = bc.velocity.orth(face_normal);
  } else if (name == "wall_rotation") {
    // rotation around `center` with angular velocity `omega`
    bc.type = BCondFluidType::wall;
    Vect center, omega;
    arg >> center >> omega;
    bc.velocity = omega.cross(face_center - center);
    bc.velocity = bc.velocity.orth(face_normal);
  } else if (name == "wall_rotation_magn") {
    // rotation around `center` with angular velocity `omega`
    // and velocity magnitude normalized to equal the magnitude of `omega`
    bc.type = BCondFluidType::wall;
    Vect center, omega;
    arg >> center >> omega;
    bc.velocity = omega.cross(face_center - center);
    bc.velocity = bc.velocity.orth(face_normal);
    const Scal norm = bc.velocity.norm();
    if (norm != 0) {
      bc.velocity *= omega.norm() / norm;
    }
  } else if (name == "slipwall") {
    bc.type = BCondFluidType::slipwall;
    bc.velocity = Vect(0);
  } else if (name == "inlet") {
    bc.type = BCondFluidType::inlet;
    arg >> bc.velocity;
  } else if (name == "inletflux") {
    bc.type = BCondFluidType::inletflux;
    arg >> bc.velocity >> bc.inletflux_id;
  } else if (name == "inletpressure") {
    bc.type = BCondFluidType::inletpressure;
    arg >> bc.pressure;
  } else if (name == "outlet") {
    bc.type = BCondFluidType::outlet;
  } else if (name == "outletpressure") {
    bc.type = BCondFluidType::outletpressure;
    arg >> bc.pressure;
  } else if (name == "symm") {
    bc.type = BCondFluidType::symm;
  } else {
    return {false, bc};
  }
  return {true, bc};
}

template <class M_>
struct UInitEmbedBc {
  using M = M_;
  using EB = Embed<M>;
  static constexpr size_t dim = M::dim;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  // Parses boundary conditions set on groups of primitives in stream.
  // Returns:
  // me_group: index of group of primitives
  // me_nci: neighbor cell id
  // vdesc: descriptors of boundary conditions for each group
  template <class MEB>
  static std::tuple<
      MapEmbed<size_t>, MapEmbed<size_t>, std::vector<std::string>>
  ParseGroups(std::istream& fin, const MEB& eb) {
    auto& m = eb.GetMesh();
    const std::vector<CodeBlock> bb = ParseCodeBlocks(fin);
    std::vector<std::string> vdesc;
    MapEmbed<size_t> me_group;
    MapEmbed<size_t> me_nci;
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
        if (m.IsBoundary(f, nci)) {
          auto ls = lsmax(eb.GetFaceCenter(f));
          if (ls > 0) {
            me_group[f] = group;
            me_nci[f] = nci;
          }
        }
      }
      for (auto c : eb.SuCFaces()) {
        auto ls = lsmax(eb.GetFaceCenter(c));
        if (ls > 0) {
          me_group[c] = group;
          me_nci[c] = 0;
        }
      }
    }
    return {me_group, me_nci, vdesc};
  }

  // Returns boundary conditions from descriptors.
  // me_group: index of group of primitives
  // vdesc: descriptors of boundary conditions for each group
  template <class EB>
  static MapEmbed<BCondFluid<Vect>> GetBCondFromGroups(
      const MapEmbed<size_t>& me_group, const MapEmbed<size_t>& me_nci,
      const std::vector<std::string> vdesc, const EB& eb) {
    MapEmbed<BCondFluid<Vect>> mebc;
    for (size_t group = 0; group < vdesc.size(); ++group) {
      me_group.LoopPairs([&](const auto& p) {
        const auto cf = p.first;
        auto bc = ParseBCondFluid<Vect>(
            vdesc[me_group.at(cf)], me_nci.at(cf), eb.GetFaceCenter(cf),
            eb.GetNormal(cf));
        mebc[cf] = bc.second;
      });
    }
    return mebc;
  }

  struct PlainBc {
    MapEmbed<BCond<Scal>> mebc;
    MapEmbed<size_t> me_group;
    std::vector<std::map<std::string, Scal>> vextra;
  };

  template <class MEB>
  static PlainBc GetPlainBc(
      const std::string& bc_path, const MEB& eb,
      const std::set<std::string>& extra_keys) {
    using UI = UInitEmbedBc<M>;

    auto get_key_value = [&](std::string s) {
      std::stringstream buf(s);
      std::string key;
      Scal value;
      buf >> key;
      buf >> value;
      return std::pair<std::string, Scal>({key, value});
    };

    std::set<std::string> known_keys = {"dirichlet", "neumann"};

    std::stringstream buf;
    {
      std::stringstream path(bc_path);
      std::string filename;
      path >> filename;
      if (filename == "inline") {
        buf << path.rdbuf();
      } else {
        filename = path.str();
        std::ifstream fin(filename);
        if (!fin.good()) {
          throw std::runtime_error(
              FILELINE + ": Can't open boundary conditions '" + filename + "'");
        }
        buf << fin.rdbuf();
      }
    }

    PlainBc res;

    std::vector<std::string> vdesc;
    MapEmbed<size_t> me_nci;
    std::tie(res.me_group, me_nci, vdesc) = UI::ParseGroups(buf, eb);
    for (auto desc : vdesc) {
      std::map<std::string, Scal> map;
      for (std::string s : Split(desc, ',')) {
        auto p = get_key_value(s);
        if (extra_keys.count(p.first)) {
          map[p.first] = p.second;
        } else {
          fassert(known_keys.count(p.first), "Unknown bc key: " + p.first);
        }
      }
      res.vextra.push_back(map);
    }

    res.me_group.LoopPairs([&](const auto& cfbc) {
      const auto cf = cfbc.first;
      const auto group = cfbc.second;
      auto& bc = res.mebc[cf];
      bc.nci = me_nci[cf];
      for (std::string s : Split(vdesc[group], ',')) {
        auto p = get_key_value(s);
        const std::string key = p.first;
        const Scal value = p.second;
        if (key == "dirichlet") {
          bc.type = BCondType::dirichlet;
          bc.val = value;
        } else if (key == "neumann") {
          bc.type = BCondType::neumann;
          bc.val = value;
        }
      }
    });
    return res;
  }
};
