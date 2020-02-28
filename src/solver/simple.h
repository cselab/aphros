// Created by Petr Karnakov on 27.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <memory>

#include "fluid.h"
#include "util/convdiff.h"

template <class M_>
class Simple final : public FluidSolver<M_> {
 public:
  using M = M_;
  using P = FluidSolver<M>; // parent
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  static constexpr size_t dim = M::dim;

  struct Par {
    Scal vrelax = 0.8; // velocity relaxation factor [0,1]
    Scal prelax = 1.; // pressure relaxation factor [0,1]
    Scal rhie = 1.; // Rhie-Chow factor [0,1] (0 disable, 1 full)
    bool second = true; // second order in time
    bool simpler = false; // Use SIMPLER
    Vect meshvel = Vect(0); // relative mesh velocity
    size_t inletflux_numid = 0; // reduction for id from 0 to numid-1
    ConvSc convsc = ConvSc::quick; // convection scheme
    Scal convdf = 1.; // deferred correction factor
    bool stokes = false;
    bool explconv = false; // explicit convection for Conv::imp
    bool convsymm = false; // symmetric solver for linear system in convdiff
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
  Simple(
      M& m, const FieldCell<Vect>& fcw, const MapEmbed<BCondFluid<Vect>>& mebc,
      const MapCell<std::shared_ptr<CondCellFluid>>& mcc, FieldCell<Scal>* fcr,
      FieldCell<Scal>* fcd, FieldCell<Vect>* fcf, FieldFace<Scal>* ffbp,
      FieldCell<Scal>* fcsv, FieldCell<Scal>* fcsm, double t, double dt,
      Par par);
  ~Simple();
  const Par& GetPar() const;
  void SetPar(Par);
  void StartStep() override;
  void MakeIteration() override;
  void FinishStep() override;
  const FieldCell<Vect>& GetVelocity(Step) const override;
  using P::GetVelocity;
  const FieldCell<Scal>& GetPressure(Step) const override;
  using P::GetPressure;
  const FieldEmbed<Scal>& GetVolumeFlux(Step) const override;
  using P::GetVolumeFlux;
  double GetAutoTimeStep() const override;
  double GetError() const override;
  // Returns velocity boundary conditions
  const MapEmbed<BCond<Vect>>& GetVelocityCond() const override;

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};
