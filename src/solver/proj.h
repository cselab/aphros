// Created by Petr Karnakov on 30.09.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <memory>

#include "fluid.h"
#include "linear/linear.h"
#include "util/convdiff.h"

template <class Vect>
struct ProjPar {
  using Scal = typename Vect::value_type;
  Scal vrelax = 1; // velocity relaxation factor [0,1]
  Scal prelax = 1.; // pressure relaxation factor [0,1]
  bool second = true; // second order in time
  Vect meshvel = Vect(0); // relative mesh velocity
  size_t inletflux_numid = 0; // reduction for id from 0 to numid-1
  ConvSc convsc = ConvSc::quick; // convection scheme
  Scal convdf = 1.; // deferred correction factor
  bool stokes = false; // Stokes flow
  bool convsymm = false; // symmetric solver for linear system in convdiff
  Conv conv = Conv::imp; // convection-diffusion solver
  bool bcg = true; // Bell-Colella-Glaz scheme
  Scal outlet_relax = 1;
  Scal inletpressure_factor =
      0; // correction factor on inlet with given pressure
  bool explviscous = false; // enable explicit viscous terms
  bool redistr_adv = false; // use RedistributeCutCellsAdvection()
                            // if true else RedistributeCutCells()
  size_t diffusion_iters = 1; // number of iterations in implicit diffusion
  bool diffusion_consistent_guess = false;
};

template <class M>
struct ProjArgs {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  const FieldCell<Vect>& fcvel;
  MapEmbed<BCondFluid<Vect>>& mebc;
  const MapCell<std::shared_ptr<CondCellFluid>>& mcc;
  const FieldCell<Scal>* fcr;
  const FieldCell<Scal>* fcd;
  const FieldCell<Vect>* fcf;
  const FieldEmbed<Scal>* ffbp;
  const FieldCell<Scal>* fcsv;
  const FieldCell<Scal>* fcsm;
  double t;
  double dt;
  std::shared_ptr<linear::Solver<M>> linsolver;
  ProjPar<Vect> par;
};

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
  using Par = ProjPar<Vect>;
  using Args = ProjArgs<M>;

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
  Proj(M& m, const EB& eb, const Args& args);
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
  const MapEmbed<BCond<Vect>>& GetVelocityCond() const override;

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};
