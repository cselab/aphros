// Created by Petr Karnakov on 29.06.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>

#include "func/primlist.h"
#include "geom/block.h"
#include "geom/field.h"
#include "geom/vect.h"
#include "parse/vars.h"

// Reads data in plain format.
// dat: data stream
// Output:
// u: field
// size: mesh size
template <class Scal, size_t dim>
void ReadPlain(
    std::istream& dat, FieldCell<Scal>& u, generic::MIdx<dim>& size) {
  dat >> size;
  GIndex<IdxCell, dim> bc(size);
  u.Reinit(bc, 0);
  for (auto c : bc.Range()) {
    dat >> u[c];
  }
}

// Init of color function.
// var: parameters
// M: mesh
// Returns:
// std::function<void(GField<Cell>& fc,const M& m)>
// fc: field to fill [i], resized if needed
// m: mesh
template <class M>
std::function<void(
    FieldCell<typename M::Scal>&, const FieldCell<typename M::Scal>&, const M&)>
CreateInitCl(const Vars& var, bool verb) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  static constexpr Scal kClNone = -1.; // no color

  std::string v = var.String["init_vf"];
  if (v == "list") {
    const size_t edim = var.Int["dim"];
    std::stringstream buf;
    {
      std::stringstream path(var.String["list_path"]);
      std::string filename;
      path >> filename;
      if (filename == "inline") {
        if (verb) {
          std::cerr
              << "InitCl: Reading inline list of primitives from list_path"
              << std::endl;
        }
        buf << path.rdbuf();
      } else {
        filename = path.str();
        std::ifstream fin(filename);
        if (verb) {
          std::cerr << "InitCl: Open list of primitives '" << filename
                    << std::endl;
        }
        fassert(fin.good(), "Can't open list of primitives '" + filename + "'");
        buf << fin.rdbuf();
      }
    }

    // TODO revise with bcast and filter by bounding box,
    // but keep original index for color
    auto pp = UPrimList<Vect>::GetPrimitives(buf, edim);

    return [pp](FieldCell<Scal>& cl, const FieldCell<Scal>& vf, const M& m) {
      cl.Reinit(m, kClNone);
      if (pp.empty()) {
        return;
      }

      for (auto c : m.Cells()) {
        if (vf[c] > 0.) {
          auto x = m.GetCenter(c);
          Scal lmax = -std::numeric_limits<Scal>::max(); // maximum level-set
          Scal clmax = 0; // color of maximum
          for (size_t i = 0; i < pp.size(); ++i) {
            auto& p = pp[i];
            Scal li = p.ls(x);
            if (p.mod_minus) {
              li = -li;
            }
            if (p.mod_and && li < 0 && lmax >= 0) {
              lmax = li;
              clmax = kClNone;
            } else {
              if (li > lmax) {
                lmax = li;
                clmax = i;
              }
            }
          }
          cl[c] = clmax;
        }
      }
    };
  } else if (v == "box") {
    Vect xc(var.Vect["box_c"]);
    Scal s = var.Double["box_s"];
    return [xc, s](FieldCell<Scal>& fc, const FieldCell<Scal>&, const M& m) {
      fc.Reinit(m, kClNone);
      for (auto c : m.Cells()) {
        if ((xc - m.GetCenter(c)).norminf() < s * 0.5) {
          fc[c] = 1.;
        }
      }
    };
  } else if (v == "grid") {
    using MIdx = typename M::MIdx;
    MIdx gridstep;
    for (auto d : M::dirs) {
      gridstep[d] =
          var.Int[std::string("initvf_grid_d") + GDir<M::dim>(d).letter()];
    }
    return [gridstep](FieldCell<Scal>& fc, const FieldCell<Scal>&, const M& m) {
      fc.Reinit(m, kClNone);
      const auto& bc = m.GetIndexCells();
      const MIdx globalsize = m.GetGlobalSize();
      const MIdx gridsize = (globalsize + gridstep - MIdx(1)) / gridstep;
      GIndex<size_t, M::dim> indexer(gridsize);
      for (auto c : m.Cells()) {
        fc[c] = indexer.GetIdx(bc.GetMIdx(c) / gridstep);
      }
    };
  } else if (v == "readplain") {
    const std::string fn = var.String["readplain_path"];
    std::ifstream qdat(fn);

    fassert(qdat.good(), "readplain: can't open '" + fn + "'");
    FieldCell<Scal> qfccl;
    using MIdx = typename M::MIdx;
    MIdx qsize;
    ReadPlain(qdat, qfccl, qsize);
    if (verb) {
      std::cerr << "Read field of size " << qsize << " from '" << fn << "'"
                << std::endl;
    }
    return [qfccl, qsize](
               FieldCell<Scal>& fc, const FieldCell<Scal>&, const M& m) {
      fc.Reinit(m, kClNone);
      const GIndex<IdxCell, M::dim> qbc(qsize);
      auto& bc = m.GetIndexCells();
      const MIdx size = m.GetGlobalSize();
      for (auto c : m.Cells()) {
        const MIdx w = bc.GetMIdx(c);
        const MIdx qw = w * qsize / size;
        const IdxCell qc = qbc.GetIdx(qw);
        fc[c] = qfccl[qc];
      }
    };
  } else if (v == "hdf") {
    return [](FieldCell<Scal>& fc, const FieldCell<Scal>&, const M& m) {
      fc.Reinit(m, 0);
    };
  } else if (v == "zero") {
    return [](FieldCell<Scal>& fc, const FieldCell<Scal>&, const M& m) {
      fc.Reinit(m, kClNone);
    };
  }
  if (verb) {
    std::cerr << "Info: no color function for init_vf=" << v << std::endl;
  }
  return std::function<void(
      FieldCell<typename M::Scal>&, const FieldCell<typename M::Scal>&,
      const M&)>();
}
