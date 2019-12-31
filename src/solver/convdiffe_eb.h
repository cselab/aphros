#pragma once

#include <memory>

#include "convdiff.h"
#include "embed.h"

template <class M_>
class ConvDiffScalExpEmbed final : public ConvDiffScal<M_> {
 public:
  using M = M_;
  using P = ConvDiffScal<M>;
  using Scal = typename M::Scal;
  using Par = typename P::Par;

  // Constructor.
  // eb: embedded boundaries
  // fcu: initial field
  // mfc: face conditions
  // bc: boundary conditions, 0: value, 1: gradient
  // bcu: value or grad.dot.outer_normal
  // fcr: density
  // ffd: diffusion
  // fcs: source
  // ffv: volume flux
  // t: initial time
  // dt: time step
  // par: parameters
  ConvDiffScalExpEmbed(
      M& m, const Embed<M>& eb, const FieldCell<Scal>& fcu,
      const MapCondFace& mfc, size_t bc, Scal bcu, const FieldCell<Scal>* fcr,
      const FieldEmbed<Scal>* fed, const FieldCell<Scal>* fcs,
      const FieldEmbed<Scal>* fev, double t, double dt, Par par);
  ~ConvDiffScalExpEmbed();
  const FieldCell<Scal>& GetField(Step) const override;
  using P::GetField;
  void Assemble(const FieldCell<Scal>&, const FieldFace<Scal>&) override;
  void CorrectField(Step l, const FieldCell<Scal>& uc) override;
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
