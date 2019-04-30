#pragma once

#include "convdiff.h"
#include "solver.h"

namespace solver {

template <class M_>
class ConvDiffVect : public UnsteadyIterativeSolver {
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  static constexpr size_t dim = M::dim;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;

 protected:
  M& m;
  const FieldCell<Scal>* fcr_;  // density
  const FieldFace<Scal>* ffd_;  // dynamic viscosity
  const FieldCell<Vect>* fcs_;  // force
  const FieldFace<Scal>* ffv_;  // volume flux

 public:
  ConvDiffVect(
      double t, double dt, M& m,
      const FieldCell<Scal>* fcr /*density*/,
      const FieldFace<Scal>* ffd /*dynamic viscosity*/,
      const FieldCell<Vect>* fcs /*force*/,
      const FieldFace<Scal>* ffv /*volume flux*/)
      : UnsteadyIterativeSolver(t, dt)
      , m(m), fcr_(fcr), ffd_(ffd), fcs_(fcs), ffv_(ffv)
  {}
  virtual const FieldCell<Vect>& GetVelocity(Layers) const = 0;
  virtual const FieldCell<Vect>& GetVelocity() const {
    return GetVelocity(Layers::time_curr);
  }
  virtual void CorrectVelocity(Layers, const FieldCell<Vect>&) = 0;
  // Returns the diagonal coefficient of the equation in direction d
  virtual FieldCell<Scal> GetDiag(size_t d) const = 0;
  // Returns the constant term of the equation in direction d
  virtual FieldCell<Scal> GetConst(size_t d) const = 0;
  // Assembles linear system
  // fcw: current velocity
  // ffv: volume flux
  // Output:
  // linear system returned by GetVelocityEquations()
  virtual void Assemble(
      const FieldCell<Vect>& fcw, const FieldFace<Scal>& ffv) = 0;
};


} // namespace solver
