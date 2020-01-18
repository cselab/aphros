// Created by Petr Karnakov on 18.01.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <memory>

#include "fluid.h"
#include "linear/linear.h"
#include "util/convdiff.h"

template <class M_>
class ProjEmbed final : public FluidSolver<M_> {
 public:
  using M = M_;
  using P = FluidSolver<M>; // parent
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  static constexpr size_t dim = M::dim;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;

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
  ProjEmbed(
      M& m, const FieldCell<Vect>& fcw, MapCondFaceFluid& mfc,
      const Embed<M>& eb, size_t bc, Vect bcvel,
      const MapCell<std::shared_ptr<CondCellFluid>>& mcc, FieldCell<Scal>* fcr,
      FieldCell<Scal>* fcd, FieldCell<Vect>* fcf, FieldFace<Scal>* ffbp,
      FieldCell<Scal>* fcsv, FieldCell<Scal>* fcsm, double t, double dt,
      Par par);
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
