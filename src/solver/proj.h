// Created by Petr Karnakov on 30.09.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <memory>

#include "fluid.h"
#include "util/convdiff.h"

template <class EB_>
class Proj final : public FluidSolver<typename EB_::M> {
 public:
  using EB = EB_;
  using M = typename EB::M;
  using Base = FluidSolver<M>;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  static constexpr size_t dim = M::dim;
  template <class T>
  using FieldFaceb = typename EmbedTraits<EB>::template FieldFaceb<T>;

  struct Par {
    Scal vrelax = 1; // velocity relaxation factor [0,1]
    Scal prelax = 1.; // pressure relaxation factor [0,1]
    bool second = true; // second order in time
    Vect meshvel = Vect(0); // relative mesh velocity
    size_t inletflux_numid = 0; // reduction for id from 0 to numid-1
    ConvSc convsc = ConvSc::quick; // convection scheme
    Scal convdf = 1.; // deferred correction factor
    bool linreport = false; // report linear solvers
    Conv conv = Conv::imp; // convection-diffusion solver
  };

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
  Proj(
      M& m, const EB& eb, const FieldCell<Vect>& fcw, MapCondFaceFluid& mfc,
      const MapCell<std::shared_ptr<CondCellFluid>>& mcc, FieldCell<Scal>* fcr,
      FieldCell<Scal>* fcd, FieldCell<Vect>* fcf, FieldFace<Scal>* ffbp,
      FieldCell<Scal>* fcsv, FieldCell<Scal>* fcsm, double t, double dt,
      Par par);
  ~Proj();
  const Par& GetPar() const;
  void SetPar(Par);
  void StartStep() override;
  void MakeIteration() override;
  void FinishStep() override;
  const FieldCell<Vect>& GetVelocity(Step) const override;
  using Base::GetVelocity;
  const FieldCell<Scal>& GetPressure(Step) const override;
  using Base::GetPressure;
  const FieldEmbed<Scal>& GetVolumeFlux(Step) const override;
  using Base::GetVolumeFlux;
  double GetAutoTimeStep() const override;
  double GetError() const override;
  const MapCondFace& GetVelocityCond() const override;

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};
