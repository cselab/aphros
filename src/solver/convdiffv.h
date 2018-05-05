#pragma once

#include "convdiff.h"

namespace solver {

template <class M_>
class ConvectionDiffusion : public UnsteadyIterativeSolver {
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  static constexpr size_t dim = M::dim;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;

 protected:
  FieldCell<Scal>* fcr_;  // density
  FieldFace<Scal>* ffd_;  // dynamic viscosity
  FieldCell<Vect>* fcs_;  // force
  FieldFace<Scal>* ffv_;  // volume flux

 public:
  ConvectionDiffusion(double t, double dt,
                      const FieldCell<Scal>* fcr,
                      const FieldFace<Scal>* ffd,
                      const FieldCell<Vect>* fcs,
                      const FieldFace<Scal>* ffv
                      )
      : UnsteadyIterativeSolver(t, dt)
      , fcr_(fcr), ffd_(ffd), fcs_(fcs), ffv_(ffv)
  {}
  virtual const FieldCell<Vect>& GetVelocity() = 0;
  virtual const FieldCell<Vect>& GetVelocity(Layers) = 0;
  virtual void CorrectVelocity(Layers, const FieldCell<Vect>&) = 0;
  virtual const FieldCell<Expr>& GetVelocityEquations(size_t /*comp*/) = 0;
};


} // namespace solver
