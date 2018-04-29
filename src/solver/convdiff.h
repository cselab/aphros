#pragma once

#include "solver.h"

namespace solver {

// Solver for convection-diffusion equation
// du/dt + div(v * u) = 1/r(d * grad^2 u + s)
// where:
// u: scalar field solved for
// v: velocity field
// d: diffusion rate
template <class M_>
class ConvectionDiffusionScalar : public UnsteadyIterativeSolver {
  using M = M_;
  using Scal = typename M::Scal;
  static constexpr size_t dim = M::dim;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;

 protected:
  const FieldCell<Scal>* fcr_;  // density
  const FieldFace<Scal>* fcd_;  // diffusion 
  const FieldCell<Scal>* fcs_;  // source
  const FieldFace<Scal>* ffv_;  // volume flux

 public:
  ConvectionDiffusionScalar(
      double t, double dt,
      const FieldCell<Scal>* fcr,
      const FieldFace<Scal>* ffd,
      const FieldCell<Scal>* fcs,
      const FieldFace<Scal>* ffv)
      : UnsteadyIterativeSolver(t, dt)
      , fcr_(fcr) , fcd_(fcd) , fcs_(fcs) , fcv_(fcv)
  {}
  virtual const FieldCell<Scal>& GetField() = 0;
  virtual const FieldCell<Scal>& GetField(Layers layer) = 0;
  virtual void CorrectField(Layers layer,
                            const FieldCell<Scal>& fcc /*correction*/) = 0;
  virtual const FieldCell<Expr>& GetEquations() = 0;
};


} // namespace solver
