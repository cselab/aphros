#pragma once

#include <memory>
#include "geom/range.h"
#include "solver/multi.h"

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

 public:
  struct Imp;
  std::unique_ptr<Imp> imp;
};

} // namespace solver
