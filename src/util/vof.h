// Created by Petr Karnakov on 02.09.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <memory>

#include "parse/vars.h"
#include "solver/advection.h"
#include "solver/cond.h"
#include "solver/multi.h"
#include "util/module.h"

template <class M_>
class UVof {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  UVof();
  ~UVof();

  struct DumpPolyArgs {
    GRange<size_t> layers;
    Multi<const FieldCell<Scal>*> fcu = nullptr; // volume fraction
    Multi<const FieldCell<Scal>*> fccl = nullptr; // color
    Multi<const FieldCell<Vect>*> fcn = nullptr; // normal
    Multi<const FieldCell<Scal>*> fca = nullptr; // plane constant
    Multi<const FieldCell<bool>*> fci = nullptr; // interface mask
    std::string filename; // path to output VTK file
    Scal time = 0; // time used in dump report if verbose
    bool binary = true; // write binary data, else ASCII
    bool merge = true; // merge close points
    bool poly = true; // dump vtk polygons instead of lines
    bool verbose = true; // print dump report to stderr on root
    bool dump_cell_index = false; // dump global cell index `c`
    bool dump_layer = false; // dump layer index `l`
    bool dump_color = true; // dump color `cl`
    std::vector<Multi<const FieldCell<Scal>*>> extra_fields;
    std::vector<std::string> extra_names;
  };

  // Dumps PLIC polygons from multiple layers.
  // filename: filename
  // t: time
  // bin: binary vtk
  // merge: merge close points
  void DumpPoly(const DumpPolyArgs& args, M& m);

  // Dumps marching cube triangles from multiple layers.
  // filename: filename
  // t: time
  // bin: binary vtk
  // merge: merge close points
  // iso: isovalue for surface fcu=iso
  // fcus [a]: sum of volume fractions, add triangles from SuCells if not null
  void DumpPolyMarch(
      const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<const FieldCell<Scal>*>& fccl,
      const Multi<const FieldCell<Vect>*>& fcn, std::string filename, Scal time,
      bool poly, bool bin, bool merge, Scal iso, const FieldCell<Scal>* fcus,
      M& m);

  // Computes unique color for each connected component over all layers.
  // fcu: volume fraction
  // fccl: color to update
  // fccl_stable: known colors to keep (may be same as fccl)
  // clfixed: if >=0, override value for color in cell nearest to clfixed_x
  // coalth: merge two layers i,j if  u_i + u_j > coalth
  // mfcu: boundary conditions for u
  // verb: report color overflow (not enough layers)
  // unionfind: use union-find algorithm (otherwise iterative stencil updates)
  // reduce: reduce color space trying to keep the previous color
  static void Recolor(
      const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<FieldCell<Scal>*>& fccl,
      const Multi<const FieldCell<Scal>*>& fccl_stable, Scal clfixed,
      Vect clfixed_x, Scal coalth, const MapEmbed<BCond<Scal>>& mfcu, bool verb,
      bool unionfind, bool reduce, bool grid, M& m);

  static std::tuple<
      MapEmbed<BCond<Scal>>, MapEmbed<BCond<Scal>>, MapEmbed<BCond<Scal>>,
      MapEmbed<BCond<Vect>>, MapEmbed<BCond<Scal>>>
  GetAdvectionBc(const M& m, const MapEmbed<BCondAdvection<Scal>>& mfc);

  // set volume fraction to 0 or 1 near wall
  static void BcClear(
      FieldCell<Scal>& uc, const MapEmbed<BCondAdvection<Scal>>& mfc,
      const M& m) {
    for (const auto& it : mfc.GetMapFace()) {
      auto& cfa = it.second;
      const IdxCell c = m.GetCell(it.first, cfa.GetNci());
      auto& u = uc[c];
      if (u < cfa.clear0) {
        u = 0;
      } else if (u > cfa.clear1) {
        u = 1;
      }
    }
  }

  // set volume fraction to 0 or 1 near wall
  static void BcClearOverrideColor(
      FieldCell<Scal>& uc, FieldCell<Scal>& clc, Scal cl0,
      const MapEmbed<BCondAdvection<Scal>>& mfc, const M& m) {
    static constexpr Scal kClNone = -1;
    for (const auto& it : mfc.GetMapFace()) {
      auto& cfa = it.second;
      const IdxCell c = m.GetCell(it.first, cfa.GetNci());
      auto& u = uc[c];
      if (u < cfa.clear0) {
        u = 0;
      } else if (u > cfa.clear1) {
        u = 1;
      }
      auto& cl = clc[c];
      if (cfa.clear1 < 1 && cl != cl0) {
        cl = kClNone;
        u = 0;
      }
    }
  }

 public:
  struct Imp;
  std::unique_ptr<Imp> imp;
};

template <class M>
class Labeling {
 public:
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  static constexpr Scal kClNone = -1;

  struct Conf {
    bool verbose = false; // report statistics of algorithm
    bool reduce = true; // reduce color space trying to keep previous colors
    bool grid = true; // use grid heuristic to speed up propagation
    Scal coalth = 1e10; // coalescence is triggered if the sum of volume
                        // fractions in two layers exceeds this value
    Scal clfixed = kClNone; // clfixed: if >=0, override value for color in cell
                            // nearest to clfixed_x
    Vect clfixed_x{0};
  };

  Labeling(const Conf& conf_) : conf(conf_) {}
  virtual ~Labeling() = default;
  virtual void SetConf(const Conf& c) {
    conf = c;
  }
  virtual const Conf& GetConf() {
    return conf;
  }
  // Computes unique color for each connected component over all layers.
  // fcu: volume fraction
  // fccl: color to update
  // fccl_stable: known colors to keep (may be same as fccl)
  // mfcu: boundary conditions for u
  // reduce: reduce color space trying to keep the previous color
  virtual void Recolor(
      const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<FieldCell<Scal>*>& fccl,
      const Multi<const FieldCell<Scal>*>& fccl_stable,
      const MapEmbed<BCond<Scal>>& mfcu, M& m) = 0;

 protected:
  Conf conf;
};

template <class M>
class ModuleLabeling : public Module<ModuleLabeling<M>> {
 public:
  using Vect = typename M::Vect;
  using Conf = typename Labeling<M>::Conf;
  using Module<ModuleLabeling>::Module;
  virtual std::unique_ptr<Labeling<M>> Make(const Conf&, const M& m) = 0;
  static Conf ParseConf(const Vars& var);
};

// Returns values over stencil centered at cell c.
// Values for neighbors without color cl are filled with 0.
// sw: stencil half-width
template <class T, class M>
std::array<T, M::Pow(3, M::dim)> GetStencilValues(
    const FieldCell<T>& fc, IdxCell c, const M& m) {
  std::array<T, M::Pow(3, M::dim)> res;
  size_t i = 0;
  for (auto cn : m.Stencil(c)) {
    res[i++] = fc[cn];
  }
  return res;
}

// Returns values over stencil centered at cell c with color cl.
// Values for neighbors without color cl are filled with 0.
// sw: stencil half-width
template <class T, class M>
std::array<T, M::Pow(3, M::dim)> GetStencilValues(
    const GRange<size_t>& layers, const Multi<const FieldCell<T>*>& fc,
    const Multi<const FieldCell<typename M::Scal>*>& fccl, IdxCell c,
    typename M::Scal cl, const M& m) {
  std::array<T, M::Pow(3, M::dim)> res;
  size_t i = 0;
  for (auto cn : m.Stencil(c)) {
    T u(0);
    for (auto l : layers) {
      if ((*fccl[l])[cn] == cl) {
        u = (*fc[l])[cn];
        break;
      }
    }
    res[i++] = u;
  }
  return res;
}
