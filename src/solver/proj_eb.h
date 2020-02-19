// Created by Petr Karnakov on 18.01.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <memory>

#include "fluid.h"
#include "solver/proj.h"
#include "util/convdiff.h"

template <class EB_>
class ProjEmbed final : public FluidSolver<typename EB_::M> {
 public:
  using EB = EB_;
  using M = typename EB::M;
  using P = FluidSolver<M>; // parent
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  static constexpr size_t dim = M::dim;
  using Par = typename Proj<M>::Par;
  template <class T>
  using FieldFaceb = typename EmbedTraits<EB>::template FieldFaceb<T>;

  // Constructor.
  // fcw: initial velocity
  // mfc: face conditions
  // mcc: cell conditions
  // fcr: density
  // fcd: dynamic viscosity
  // fcf: force
  // ffbp: projections of balanced force
  // fcsv: volume source
  // fcsm: mass source
  // t: initial time
  // dt: time step
  // par: parameters
  ProjEmbed(
      M& m, const Embed<M>& eb, const FieldCell<Vect>& fcw,
      MapCondFaceFluid& mfc, const MapCell<std::shared_ptr<CondCellFluid>>& mcc,
      const FieldCell<Scal>* fcr, const FieldCell<Scal>* fcd,
      const FieldCell<Vect>* fcf, const FieldFace<Scal>* ffbp,
      const FieldCell<Scal>* fcsv, const FieldCell<Scal>* fcsm, double t,
      double dt, Par par);
  ~ProjEmbed();
  const Par& GetPar() const;
  void SetPar(Par);
  // ...
  void StartStep() override;
  // ...
  void MakeIteration() override;
  // ...
  void FinishStep() override;
  // ...
  const FieldCell<Vect>& GetVelocity(Step) const override;
  // ...
  using P::GetVelocity;
  // ...
  const FieldCell<Scal>& GetPressure(Step) const override;
  // ...
  using P::GetPressure;
  // ...
  const FieldFace<Scal>& GetVolumeFlux(Step) const override;
  // ...
  using P::GetVolumeFlux;
  // ...
  double GetAutoTimeStep() const override;
  // ...
  double GetError() const override;
  // ...
  const MapCondFace& GetVelocityCond() const override;

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};
