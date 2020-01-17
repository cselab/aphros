// Created by Petr Karnakov on 17.01.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <memory>

#include "convdiffv.h"
#include "embed.h"

template <class M_>
class ConvDiffVectEmbed final : public ConvDiffVect<M_> {
 public:
  using M = M_;
  using P = ConvDiffVect<M>; // parent
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  static constexpr size_t dim = M::dim;
  using CD = ConvDiffScalExpEmbed<M_>;
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
  ConvDiffVectEmbed(
      M& m, const Embed<M>& eb, const FieldCell<Vect>& fcvel,
      const MapCondFace& mfc, size_t bc, Scal bcu, const FieldCell<Scal>* fcr,
      const FieldEmbed<Scal>* fed, const FieldCell<Vect>* fcs,
      const FieldFace<Scal>* ffv, double t, double dt, Par par);
  ~ConvDiffVectEmbed();
  // ...
  void Assemble(
      const FieldCell<Vect>& fcw, const FieldFace<Scal>& ffv) override;
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
  using P::GetVelocity;
  // ...
  double GetError() const override;

 private:
  struct Imp; // implementation
  std::unique_ptr<Imp> imp;
};
