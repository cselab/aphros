#pragma once

#include <memory>

#include "solver/multi.h"
#include "solver/cond.h"
#include "solver/advection.h"

namespace solver {

template <class M_>
class UVof {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  UVof();
  ~UVof();

  // Dumps PLIC polygons from multiple layers.
  // fn: filename
  // t: time
  // th: threshold for volume fraction: th < u < 1 - th
  // bin: binary vtk
  // merge: merge close points
  void DumpPoly(
      const GRange<size_t>& layers,
      const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<const FieldCell<Scal>*>& fccl,
      const Multi<const FieldCell<Vect>*>& fcn,
      const Multi<const FieldCell<Scal>*>& fca,
      const Multi<const FieldCell<bool>*>& fci,
      std::string fn, Scal t, Scal th, bool bin, bool merge, M& m);

  // Dumps marching cube triangles from multiple layers.
  // fn: filename
  // t: time
  // th: threshold for volume fraction: th < u < 1 - th
  // bin: binary vtk
  // merge: merge close points
  // iso: isovalue for surface fcu=iso
  void DumpPolyMarch(
      const GRange<size_t>& layers,
      const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<const FieldCell<Scal>*>& fccl,
      const Multi<const FieldCell<Vect>*>& fcn,
      const Multi<const FieldCell<Scal>*>& fca,
      const Multi<const FieldCell<bool>*>& fci,
      std::string fn, Scal t, Scal th, bool bin, bool merge, Scal iso, 
      const FieldCell<Scal>*, M& m);

  // Dumps PLIC polygons from single layer.
  // fn: filename
  // t: time
  // th: threshold for volume fraction: th < u < 1 - th
  // bin: binary vtk
  // merge: merge close points
  void DumpPoly(
      const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcn,
      const FieldCell<Scal>& fca, const FieldCell<bool>& fci,
      std::string fn, Scal t, Scal th, bool bin, bool merge, M& m);

  // Computes unique color for each connected component over all layers.
  // fcu: volume fraction
  // fccl: color to update
  // fccl0: known colors to keep (may be same as fccl)
  // clfixed: if >=0, override value for color in cell nearest to clfixed_x
  // coalth: merge two layers i,j if  u_i + u_j > coalth
  // mfcu: boundary conditions for u
  // verb: report color overflow (not enough layers)
  // unionfind: use union-find algorithm (otherwise iterative stencil updates)
  // reduce: reduce color space trying to keep the previous color
  static void Recolor(
      const GRange<size_t>& layers,
      const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<FieldCell<Scal>*>& fccl,
      const Multi<const FieldCell<Scal>*>& fccl0,
      Scal clfixed, Vect clfixed_x, Scal coalth,
      const MapCondFace& mfcu, bool verb, bool unionfind, bool reduce,
      bool grid, M& m);

  static void GetAdvectionFaceCond(
      const M& m, const MapCondFaceAdvection<Scal>& mfc,
      MapCondFace& mfc_vf, MapCondFace& mfc_cl, MapCondFace& mfc_im,
      MapCondFace& mfc_n, MapCondFace& mfc_a);

  // set volume fraction to 0 or 1 near wall
  static void BcClear(FieldCell<Scal>& uc,
                      const MapCondFace& mfc, const M& m, Scal u0, Scal u1) {
    for (const auto& it : mfc) {
      auto& cb = it.GetValue();
      if (cb.Get<CondFaceVal<Scal>>()) {
        IdxCell c = m.GetNeighbourCell(it.GetIdx(), cb->GetNci());
        auto& u = uc[c];
        if (u < u0) {
          u = 0;
        } else if (u > u1) {
          u = 1;
        }
      }
    }
  }

 public:
  struct Imp;
  std::unique_ptr<Imp> imp;
};

} // namespace solver
