#pragma once

#include "solver.h"
#include "geom/mesh.h"
#include "cond.h"

namespace solver {

// Solver for convection-diffusion equation
// du/dt + div(v * u) = 1/r(d * grad^2 u + s)
// where:
// u: scalar field solved for
// v: velocity field
// d: diffusion rate
template <class M_>
class ConvDiffScal : public UnsteadyIterativeSolver {
 public:
  using M = M_;
  static constexpr size_t dim = M::dim;
  using Scal = typename M::Scal;

  struct Par {
    Scal relax = 1.;      // relaxation factor [0,1] (1 -- no relaxation)
    Scal guessextra = 0.; // next iteration guess extrapolation weight [0,1]
    bool second = true; // second order in time
    ConvSc sc = ConvSc::quick; // scheme for convective flux (see convdiffi.h)
    Scal df = 1.; // deferred correction factor
    Scal th = 1e-10; // threshold for flow direction
    bool linreport = false; // report linear solvers
  };

  // par: parameters
  // fcr: density
  // ffd: diffusion coefficient
  // fcs: source
  // ffv: volume flux
  ConvDiffScal(
      double t, double dt, M& m, const Par& par,
      const FieldCell<Scal>* fcr, const FieldFace<Scal>* ffd,
      const FieldCell<Scal>* fcs, const FieldFace<Scal>* ffv)
      : UnsteadyIterativeSolver(t, dt)
      , m(m), par_(par), fcr_(fcr) , ffd_(ffd) , fcs_(fcs) , ffv_(ffv)
  {}
  virtual const FieldCell<Scal>& GetField(Layers) const = 0;
  virtual const FieldCell<Scal>& GetField() const {
    return GetField(Layers::time_curr);
  }
  // Corrects field and exchanges halos.
  // l: layer to correct
  // uc: correction [i]
  virtual void CorrectField(Layers l, const FieldCell<Scal>& uc) = 0;
  // Assembles linear system
  // fcu: field from previous iteration [a]
  // ffv: volume flux
  // Output:
  // fcucs_: linear system, overwritten, returned by GetDiag and GetConst
  virtual void Assemble(
      const FieldCell<Scal>& fcu, const FieldFace<Scal>& ffv) = 0;
  // Returns the diagonal coefficient of the equation
  virtual FieldCell<Scal> GetDiag() const = 0;
  // Returns the constant term of the equation
  virtual FieldCell<Scal> GetConst() const = 0;
  virtual const Par& GetPar() const { return par_; }
  virtual void SetPar(const Par& par) { par_ = par; }

 protected:
  M& m;
  Par par_;  // parameters
  const FieldCell<Scal>* fcr_;  // density
  const FieldFace<Scal>* ffd_;  // diffusion
  const FieldCell<Scal>* fcs_;  // source
  const FieldFace<Scal>* ffv_;  // volume flux
};


} // namespace solver
