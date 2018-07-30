#pragma once

#include <memory>

#include "convdiff.h"

namespace solver {

template <class M_>
class ConvDiffScalImp : public ConvDiffScal<M_> {
 public:
  using M = M_;
  using P = ConvDiffScal<M>; // parent
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  static constexpr size_t dim = M::dim;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;

  struct Par {
    Scal relax = 1.;      // relaxation factor [0,1] (1 -- no relaxation)
    Scal guessextra = 0.; // next iteration guess extrapolation weight [0,1]
    bool second = true; // second order in time
    ConvSc sc = ConvSc::quick; // scheme for convective flux (see convdiffi.h)
    Scal df = 1.; // deferred correction factor
    Scal th = 1e-10; // threshold for flow direction
  };
  // Constructor.
  // fcu: initial field
  // mfc: face conditions
  // mcc: cell conditions
  // fcr: density
  // ffd: diffusion
  // fcs: source
  // ffv: volume flux
  // t: initial time
  // dt: time step
  // par: parameters
  ConvDiffScalImp(
      M& m, const FieldCell<Scal>& fcu, 
      const MapFace<std::shared_ptr<CondFace>>& mfc, 
      const MapCell<std::shared_ptr<CondCell>>& mcc, 
      const FieldCell<Scal>* fcr, const FieldFace<Scal>* ffd,
      const FieldCell<Scal>* fcs, const FieldFace<Scal>* ffv,
      double t, double dt, std::shared_ptr<Par> par);
  ~ConvDiffScalImp();
  // Parameters
  Par* GetPar();
  // Assembles linear system
  // fcu: field from previous iteration [a]
  // ffv: volume flux
  // Output:
  // fcucs_: linear system, overwritten, returned by GetEquations()
  void Assemble(const FieldCell<Scal>& fcu, const FieldFace<Scal>& ffv);
  // Corrects field and comm.
  // uc: correction [i]
  // Output:
  // u(l) += uc [a]
  void CorrectField(Layers l, const FieldCell<Scal>& uc) override;
  // Equations for velocity
  const FieldCell<Expr>& GetEquations() const;
  // ...
  void StartStep() override;
  // ...
  void MakeIteration() override;
  // ...
  void FinishStep() override;
  // ...
  const FieldCell<Scal>& GetField(Layers layer) const override;
  // ...
  using P::GetField;
  // ...
  double GetError() const override;

 private:
  struct Imp; // implementation
  std::unique_ptr<Imp> imp;
};

} // namespace solver
