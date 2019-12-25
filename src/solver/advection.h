#pragma once

#include "geom/mesh.h"
#include "solver/solver.h"

namespace solver {

template <class Scal>
struct CondFaceAdvection {
  enum class Halo { fill, reflect };

  CondFaceAdvection() = default;
  CondFaceAdvection(size_t nci) : nci(nci) {}
  size_t GetNci() const {
    return nci;
  }

  size_t nci;
  Scal clear0 = 0; // snap to 0 if vf<clear0
  Scal clear1 = 1; // snap to 1 if vf>clear1
  Halo halo = Halo::reflect;
  // values for halo cells if halo=fill:
  Scal fill_vf; // volume fraction
  Scal fill_cl; // color
};

template <class Scal>
using MapCondFaceAdvection = MapFace<CondFaceAdvection<Scal>>;

template <class M_>
class AdvectionSolver : public UnsteadyIterativeSolver {
  using M = M_;
  using Scal = typename M::Scal;

 protected:
  M& m;
  const FieldFace<Scal>* ffv_; // volume flux [velocity*area]
  const FieldCell<Scal>* fcs_; // source [value/time]

 public:
  // ffv: volume flux
  // fcs: source
  AdvectionSolver(
      double t, double dt, M& m, const FieldFace<Scal>* ffv,
      const FieldCell<Scal>* fcs)
      : UnsteadyIterativeSolver(t, dt), m(m), ffv_(ffv), fcs_(fcs) {}
  // Postprocessing after time step (dumps)
  virtual void PostStep() {}
  // Volume fraction
  virtual const FieldCell<Scal>& GetField(Layers) const = 0;
  // Volume fraction at last time step
  virtual const FieldCell<Scal>& GetField() const {
    return GetField(Layers::time_curr);
  }
};

} // namespace solver
