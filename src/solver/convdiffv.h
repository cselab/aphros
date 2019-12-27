#pragma once

#include "convdiff.h"
#include "solver.h"

template <class M_>
class ConvDiffVect : public UnsteadyIterativeSolver {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Par = typename ConvDiffScal<M>::Par;
  static constexpr size_t dim = M::dim;

  // par: parameters
  // fcr: density
  // ffd: diffusion coefficient
  // fcs: source
  // ffv: volume flux
  ConvDiffVect(
      double t, double dt, M& m, Par par, const FieldCell<Scal>* fcr,
      const FieldFace<Scal>* ffd, const FieldCell<Vect>* fcs,
      const FieldFace<Scal>* ffv)
      : UnsteadyIterativeSolver(t, dt)
      , m(m)
      , par(par)
      , fcr_(fcr)
      , ffd_(ffd)
      , fcs_(fcs)
      , ffv_(ffv) {}
  virtual const FieldCell<Vect>& GetVelocity(Step) const = 0;
  virtual const FieldCell<Vect>& GetVelocity() const {
    return GetVelocity(Step::time_curr);
  }
  virtual void CorrectVelocity(Step, const FieldCell<Vect>&) = 0;
  // Assembles linear system
  // fcw: current velocity
  // ffv: volume flux
  // Output:
  // linear system returned by GetDiag() and GetConst()
  virtual void Assemble(
      const FieldCell<Vect>& fcw, const FieldFace<Scal>& ffv) = 0;
  // Returns the diagonal coefficient of the equation in direction d
  virtual FieldCell<Scal> GetDiag(size_t d) const = 0;
  // Returns the constant term of the equation in direction d
  virtual FieldCell<Scal> GetConst(size_t d) const = 0;
  virtual const Par& GetPar() const {
    return par;
  }
  virtual void SetPar(Par par0) {
    par = par0;
  }

 protected:
  M& m;
  Par par; // parameters
  const FieldCell<Scal>* fcr_; // density
  const FieldFace<Scal>* ffd_; // dynamic viscosity
  const FieldCell<Vect>* fcs_; // force
  const FieldFace<Scal>* ffv_; // volume flux
};
