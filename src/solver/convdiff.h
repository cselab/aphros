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
class ConvDiffScal : public UnsteadyIterativeSolver {
  using M = M_;
  using Scal = typename M::Scal;
  static constexpr size_t dim = M::dim;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;

 protected:
  M& m;
  const FieldCell<Scal>* fcr_;  // density
  const FieldFace<Scal>* ffd_;  // diffusion 
  const FieldCell<Scal>* fcs_;  // source
  const FieldFace<Scal>* ffv_;  // volume flux

 public:
  ConvDiffScal(
      double t, double dt, M& m,
      const FieldCell<Scal>* fcr /*density*/,
      const FieldFace<Scal>* ffd /*diffusion*/,
      const FieldCell<Scal>* fcs /*source*/,
      const FieldFace<Scal>* ffv /*volume flux*/)
      : UnsteadyIterativeSolver(t, dt)
      , m(m), fcr_(fcr) , ffd_(ffd) , fcs_(fcs) , ffv_(ffv)
  {}
  virtual const FieldCell<Scal>& GetField(Layers) const = 0;
  virtual const FieldCell<Scal>& GetField() const {
    return GetField(Layers::time_curr);
  }
  virtual void CorrectField(Layers, const FieldCell<Scal>&) = 0;
  virtual const FieldCell<Expr>& GetEquations() const = 0;
};


} // namespace solver
