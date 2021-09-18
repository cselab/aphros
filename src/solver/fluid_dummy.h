// Created by Petr Karnakov on 30.09.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <memory>

#include "fluid.h"
#include "parse/vars.h"

template <class EB_>
class FluidDummy final : public FluidSolver<typename EB_::M> {
 public:
  using EB = EB_;
  using M = typename EB::M;
  using Base = FluidSolver<M>;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  static constexpr size_t dim = M::dim;
  template <class T>
  using FieldFaceb = typename EmbedTraits<EB>::template FieldFaceb<T>;

  // Constructor.
  // fcvel: initial velocity
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
  FluidDummy(
      M& m, const EB& eb, const FieldCell<Vect>& fcvel,
      const MapEmbed<BCondFluid<Vect>>& mebc, const FieldCell<Scal>* fcr,
      const FieldCell<Scal>* fcd, const FieldCell<Vect>* fcf,
      const FieldEmbed<Scal>* ffbp, const FieldCell<Scal>* fcsv,
      const FieldCell<Scal>* fcsm, double t, double dt, const Vars& var);
  ~FluidDummy();
  void MakeIteration() override;
  const FieldCell<Vect>& GetVelocity(Step) const override;
  using Base::GetVelocity;
  const FieldCell<Scal>& GetPressure(Step) const override;
  using Base::GetPressure;
  const FieldEmbed<Scal>& GetVolumeFlux(Step) const override;
  double GetAutoTimeStep() const override;
  using Base::GetVolumeFlux;
  const MapEmbed<BCond<Vect>>& GetVelocityCond() const override;

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};
