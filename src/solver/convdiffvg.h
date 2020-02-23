// Created by Petr Karnakov on 20.11.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <memory>

#include "convdiffv.h"
#include "linear/linear.h"

template <class EB_, class CD_>
class ConvDiffVectGeneric final : public ConvDiffVect<EB_> {
 public:
  using EB = EB_;
  using M = typename EB::M;
  using Base = ConvDiffVect<EB>;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  static constexpr size_t dim = M::dim;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;
  using CD = CD_;
  using Par = typename CD::Par;
  template <class T>
  using FieldFaceb = typename EmbedTraits<EB>::template FieldFaceb<T>;

  // Constructor.
  // fcvel: initial velocity
  // mebc: face conditions
  // fcr: density
  // ffd: dynamic viscosiity
  // fcs: source
  // ffv: volume flux
  // t: initial time
  // dt: time step
  // par: parameters
  ConvDiffVectGeneric(
      M& m, const EB& eb, const FieldCell<Vect>& fcvel,
      const MapEmbed<BCond<Vect>>& mebc, const FieldCell<Scal>* fcr,
      const FieldFaceb<Scal>* ffd, const FieldCell<Vect>* fcs,
      const FieldFaceb<Scal>* ffv, double t, double dt, Par par);
  ~ConvDiffVectGeneric();
  // ...
  void Assemble(
      const FieldCell<Vect>& fcw, const FieldFaceb<Scal>& ffv) override;
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
  MapEmbed<BCond<Scal>>& GetVelocityCond(size_t d);
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
