#pragma once

#include <memory>

#include "convdiff.h"

namespace solver {

template <class M_>
class ConvDiffScalExp : public ConvDiffScal<M_> {
 public:
  using M = M_;
  using P = ConvDiffScal<M>;
  using Scal = typename M::Scal;
  using Par = typename P::Par;

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
  ConvDiffScalExp(
      M& m, const FieldCell<Scal>& fcu, 
      const MapFace<std::shared_ptr<CondFace>>& mfc, 
      const MapCell<std::shared_ptr<CondCell>>& mcc, 
      const FieldCell<Scal>* fcr, const FieldFace<Scal>* ffd,
      const FieldCell<Scal>* fcs, const FieldFace<Scal>* ffv,
      double t, double dt, const Par& par);
  ~ConvDiffScalExp();
  const FieldCell<Scal>& GetField(Layers) const override;
  using P::GetField;
  void Assemble(const FieldCell<Scal>&, const FieldFace<Scal>&) override;
  void CorrectField(Layers l, const FieldCell<Scal>& uc) override;
  FieldCell<Scal> GetDiag() const override;
  FieldCell<Scal> GetConst() const override;
  void StartStep() override;
  void MakeIteration() override;
  void FinishStep() override;
  double GetError() const override;

 private:
  struct Imp; // implementation
  std::unique_ptr<Imp> imp;
};

} // namespace solver
