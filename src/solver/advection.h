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
  AdvectionSolver(double t, double dt, M& m,
                  const FieldFace<Scal>* ffv /*volume flux*/,
                  const FieldCell<Scal>* fcs /*source*/)
      : UnsteadyIterativeSolver(t, dt)
      , m(m), ffv_(ffv), fcs_(fcs) {}
  virtual const FieldCell<Scal>& GetField(Layers) = 0;
  virtual const FieldCell<Scal>& GetField() {
    return GetField(Layers::time_curr);
  }
  virtual const FieldCell<Scal>& GetCurv() = 0;
};


} // namespace solver
