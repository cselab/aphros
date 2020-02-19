// Created by Petr Karnakov on 17.01.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <memory>

#include "convdiffe.h"
#include "convdiffv.h"
#include "embed.h"

template <class MEB_>
class ConvDiffVectEmbed final : public ConvDiffVect<MEB_> {
 public:
  using EB = MEB_;
  using M = typename EB::M;
  using Base = ConvDiffVect<EB>;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  static constexpr size_t dim = M::dim;
  using CD = ConvDiffScalExp<Embed<M>>;
  using Par = typename CD::Par;
  template <class T>
  using FieldFaceb = typename EmbedTraits<EB>::template FieldFaceb<T>;

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
  ConvDiffVectEmbed(
      M& m, const EB& eb, const FieldCell<Vect>& fcvel,
      const MapCondFace& mfc, const FieldCell<Scal>* fcr,
      const FieldFaceb<Scal>* fed, const FieldCell<Vect>* fcs,
      const FieldFaceb<Scal>* ffv, double t, double dt, Par par);
  ~ConvDiffVectEmbed();
  // ...
  void Assemble(
      const FieldCell<Vect>& fcw, const FieldFaceb<Scal>& ffv);
  // Corrects field and comm.
  // fc: correction [i]
  // Output:
  // vel(l) += fc [a]
  void CorrectVelocity(Step l, const FieldCell<Vect>& fc) override;
  // ...
  FieldCell<Scal> GetDiag(size_t d) const override;
  // ...
  FieldCell<Scal> GetConst(size_t d) const override;
  // Velocity conditions.
  // d: component
  MapCondFace& GetVelocityCond(size_t d);
  // ...
  void StartStep() override;
  // ...
  void MakeIteration() override;
  // ...
  void FinishStep() override;
  // ...
  const FieldCell<Vect>& GetVelocity(Step l) const override;
  // ...
  using Base::GetVelocity;
  // ...
  double GetError() const override;

 private:
  struct Imp; // implementation
  std::unique_ptr<Imp> imp;
};
