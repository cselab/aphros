#pragma once

#include "solver/solver.h"

namespace solver {

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
  AdvectionSolver(double t, double dt, M& m,
                  const FieldFace<Scal>* ffv, const FieldCell<Scal>* fcs)
      : UnsteadyIterativeSolver(t, dt)
      , m(m), ffv_(ffv), fcs_(fcs) {}
  // Postprocessing after time step (curvature, dumps)
  virtual void PostStep() {}
  // Volume fraction
  virtual const FieldCell<Scal>& GetField(Layers) const = 0;
  // Volume fraction at last time step
  virtual const FieldCell<Scal>& GetField() const {
    return GetField(Layers::time_curr);
  }
  // Curvature
  virtual const FieldCell<Scal>& GetCurv() const = 0;
};


} // namespace solver
