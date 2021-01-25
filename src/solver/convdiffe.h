// Created by Petr Karnakov on 30.04.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <memory>

#include "convdiff.h"

template <class EB_>
class ConvDiffScalExp final : public ConvDiffScal<EB_> {
 public:
  using EB = EB_;
  using M = typename EB::M;
  using Base = ConvDiffScal<EB>;
  using Scal = typename M::Scal;
  using Par = typename ConvDiffScal<M>::Par;
  using Args = ConvDiffArgs<EB>;
  template <class T>
  using FieldFaceb = typename EmbedTraits<EB>::template FieldFaceb<T>;

  // Constructor.
  // fcu: initial field
  // mebc: face conditions
  // mcc: cell conditions
  // fcr: density
  // ffd: diffusion
  // fcs: source
  // ffv: volume flux
  // t: initial time
  // dt: time step
  // par: parameters
  ConvDiffScalExp(M& m, const EB& eb, const Args& args);
  ~ConvDiffScalExp();
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
