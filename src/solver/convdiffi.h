// Created by Petr Karnakov on 27.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <memory>

#include "convdiff.h"

template <class EB_>
class ConvDiffScalImp final : public ConvDiffScal<EB_> {
 public:
  using EB = EB_;
  using M = typename EB::M;
  using Base = ConvDiffScal<EB>;
  using Scal = typename M::Scal;
  using Par = typename ConvDiffScal<M>::Par;
  using Args = ConvDiffArgs<EB>;
  template <class T>
  using FieldFaceb = typename EmbedTraits<EB>::template FieldFaceb<T>;
  using Base::dim;

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
  ConvDiffScalImp(M& m, const EB& eb, const Args& args);
  ~ConvDiffScalImp();
  const FieldCell<Scal>& GetField(Step) const override;
  using Base::GetField;
  void Assemble(const FieldCell<Scal>&, const FieldFaceb<Scal>&) override;
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
