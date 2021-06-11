// Created by Petr Karnakov on 27.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <memory>

#include "fluid.h"
#include "util/convdiff.h"

template <class Vect>
struct SimplePar {
  using Scal = typename Vect::Scal;
  Scal vrelax = 0.8; // velocity relaxation factor [0,1]
  Scal prelax = 1.; // pressure relaxation factor [0,1]
  Scal rhie = 1.; // Rhie-Chow factor [0,1] (0 disable, 1 full)
  bool second = true; // second order in time
  Vect meshvel = Vect(0); // relative mesh velocity
  size_t inletflux_numid = 0; // reduction for id from 0 to numid-1
  ConvSc convsc = ConvSc::quick; // convection scheme
  Scal convdf = 1.; // deferred correction factor
  bool stokes = false;
  bool explconv = false; // explicit convection for Conv::imp
  bool convsymm = false; // symmetric solver for linear system in convdiff
  Conv conv = Conv::imp; // convection-diffusion solver
  Scal outlet_relax = 1;
  bool explviscous = true; // enable explicit viscous terms
};

template <class M>
struct SimpleArgs {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  const FieldCell<Vect>& fcvel;
  const MapEmbed<BCondFluid<Vect>>& mebc;
  const MapCell<std::shared_ptr<CondCellFluid>>& mcc;
  const FieldCell<Scal>* fcr;
  const FieldCell<Scal>* fcd;
  const FieldCell<Vect>* fcf;
  const FieldEmbed<Scal>* febp;
  const FieldCell<Scal>* fcsv;
  const FieldCell<Scal>* fcsm;
  double t;
  double dt;
  std::shared_ptr<linear::Solver<M>> linsolver_symm;
  std::shared_ptr<linear::Solver<M>> linsolver_gen;
  SimplePar<Vect> par;
};

template <class EB_>
class Simple final : public FluidSolver<typename EB_::M> {
 public:
  using EB = EB_;
  using M = typename EB::M;
  using Base = FluidSolver<M>;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  static constexpr size_t dim = M::dim;
  template <class T>
  using FieldFaceb = typename EmbedTraits<EB>::template FieldFaceb<T>;
  using Par = SimplePar<Vect>;
  using Args = SimpleArgs<M>;

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
  Simple(M& m, const EB& eb, const Args& args);
  ~Simple();
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
