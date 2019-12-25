#pragma once

#include <memory>

#include "cond.h"
#include "dump/dumper.h"
#include "geom/mesh.h"
#include "geom/range.h"
#include "multi.h"
#include "partstr.h"

namespace solver {

template <class M_>
class PartStrMeshM {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using PS = PartStr<Scal>;

  struct Par {
    std::shared_ptr<typename PS::Par> ps;
    size_t dump_fr = -1; // num frames dump
    size_t ns = 2; // number of strings per cell
    Scal tol = 0.01; // tolerance
    size_t itermax = 20;
    int verb = 0; // debug output
    size_t dim = 3;
    Scal maxr = 0; // if input radius of curvature
                   // (e.g. from height functions)
                   // is below than maxr*h,
                   // overwrite with estimate from particles
    bool vtkbin = true; // write binary vtk in DumpPartInter
    bool vtkmerge = true; // merge close points in DumpPartInter
  };

  PartStrMeshM(M& m, std::shared_ptr<Par> par, const GRange<size_t>& layers);
  ~PartStrMeshM();

  // Computes curvature with particles.
  // vfcu: volume fraction
  // vfca: plane constant
  // vfcn: normal
  // vfci: interface mask (1: contains interface)
  // vfccl: color
  void Part(
      const Multi<const FieldCell<Scal>*>& vfcu,
      const Multi<const FieldCell<Scal>*>& vfca,
      const Multi<const FieldCell<Vect>*>& vfcn,
      const Multi<const FieldCell<bool>*>& vfci,
      const Multi<const FieldCell<Scal>*>& vfccl, const FieldCell<Scal>* fck);
  // Dump particles to csv.
  // vfca: plane constant
  // vfcn: normal
  // n: frame index
  // t: time
  void DumpParticles(
      const Multi<const FieldCell<Scal>*>& vfca,
      const Multi<const FieldCell<Vect>*>& vfcn, size_t id, Scal t);
  void DumpPartInter(
      const Multi<const FieldCell<Scal>*>& vfca,
      const Multi<const FieldCell<Vect>*>& vfcn, size_t id, Scal t);
  // Returns curvature field from last call of Part()
  const FieldCell<Scal>& GetCurv(size_t l);

 private:
  struct Imp; // implementation
  // std::unique_ptr<Imp> imp;
  std::unique_ptr<Imp> imp;
};

} // namespace solver
