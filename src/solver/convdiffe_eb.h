// Created by Petr Karnakov on 31.12.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <memory>

#include "convdiff.h"
#include "embed.h"

template <class EB_>
class ConvDiffScalExpEmbed final : public ConvDiffScal<EB_> {
 public:
  using EB = EB_;
  using M = typename EB::M;
  using Base = ConvDiffScal<EB>;
  using Scal = typename M::Scal;
  using Par = typename Base::Par;
  template <class T>
  using FieldFaceb = typename EmbedTraits<EB>::template FieldFaceb<T>;

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
      M& m, const EB& eb, const FieldCell<Scal>& fcu,
      const MapCondFace& mfc, size_t bc, Scal bcu, const FieldCell<Scal>* fcr,
      const FieldFaceb<Scal>* ffd, const FieldCell<Scal>* fcs,
      const FieldFaceb<Scal>* ffv, double t, double dt, Par par);
  ~ConvDiffScalExpEmbed();
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
