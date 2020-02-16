// Created by Petr Karnakov on 27.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <memory>

#include "convdiff.h"

template <class EB_>
class ConvDiffScalImp final : public ConvDiffScal<typename EB_::M> {
 public:
  using EB = EB_;
  using M = typename EB::M;
  using P = ConvDiffScal<M>;
  using Scal = typename M::Scal;
  using Par = typename P::Par;
  using P::dim;

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
      M& m, const EB& eb, const FieldCell<Scal>& fcu, const MapCondFace& mfc,
      const FieldCell<Scal>* fcr, const FieldFace<Scal>* ffd,
      const FieldCell<Scal>* fcs, const FieldFace<Scal>* ffv, double t,
      double dt, Par par);
  ~ConvDiffScalImp();
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
  struct Imp;
  std::unique_ptr<Imp> imp;
};
