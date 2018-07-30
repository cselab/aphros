#pragma once

#include <memory>

#include "convdiffv.h"
#include "convdiffi.h"

namespace solver {

template <class M_>
class ConvDiffVectImp : public ConvDiffVect<M_> {
 public:
  using M = M_;
  using P = ConvDiffVect<M>; // parent
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  static constexpr size_t dim = M::dim;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;
  using CD = ConvDiffScalImp<M>;
  using Par = typename CD::Par;

  // Constructor.
  // fcvel: initial velocity
  // mfc: face conditions
  // mcc: cell conditions
  // fcr: density
  // ffd: dynamic viscosiity
  // fcs: source
  // ffv: volume flux
  // t: initial time
  // dt: time step
  // par: parameters
  ConvDiffVectImp(
      M& m, const FieldCell<Vect>& fcvel,
      const MapFace<std::shared_ptr<CondFace>>& mfc, 
      const MapCell<std::shared_ptr<CondCell>>& mcc,
      const FieldCell<Scal>* fcr, const FieldFace<Scal>* ffd, 
      const FieldCell<Vect>* fcs, const FieldFace<Scal>* ffv, 
      double t, double dt, std::shared_ptr<Par> par);
  ~ConvDiffVectImp();
  // Parameters
  Par* GetPar();
  // Assembles linear system
  // fcw: current velocity
  // ffv: volume flux
  // Output:
  // fcucs_: linear system, overwritten, returned by GetEquations()
  void Assemble(const FieldCell<Vect>& fcw, const FieldFace<Scal>& ffv);
  // Corrects field and comm.
  // fc: correction [i]
  // Output:
  // vel(l) += fc [a]
  void CorrectVelocity(Layers l, const FieldCell<Vect>& fc) override;
  // Equations for velocity
  const FieldCell<Expr>& GetVelocityEquations(size_t d) const override;
  // Velocity conditions.
  // d: component
  MapFace<std::shared_ptr<CondFace>>& GetVelocityCond(size_t d);
  // Instance of ConvDiffScal for component d
  CD& GetSolver(size_t d);
  // ...
  void StartStep() override;
  // ...
  void MakeIteration() override;
  // ...
  void FinishStep() override;
  // ...
  const FieldCell<Vect>& GetVelocity(Layers l) const override;
  // ...
  using P::GetVelocity;
  // ...
  double GetError() const override;

 private:
  struct Imp; // implementation
  std::unique_ptr<Imp> imp;
};

} // namespace solver

