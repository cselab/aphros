// Created by Petr Karnakov on 27.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include "cond.h"
#include "embed.h"
#include "geom/mesh.h"
#include "solver.h"

// Solver for convection-diffusion equation
// du/dt + div(v * u) = 1/r(d * grad^2 u + s)
// where:
// u: scalar field solved for
// v: velocity field
// d: diffusion rate
template <class M_>
class ConvDiffScal : public UnsteadyIterativeSolver {
 public:
  using M = M_;
  static constexpr size_t dim = M::dim;
  using Scal = typename M::Scal;

  struct Par {
    Scal relax = 1.; // relaxation factor [0,1] (1 -- no relaxation)
    bool second = true; // second order in time
    ConvSc sc = ConvSc::quick; // scheme for convective flux (see convdiffi.h)
    Scal df = 1.; // deferred correction factor
    Scal th = 1e-10; // threshold for flow direction
    bool linreport = false; // report linear solvers
  };

  // par: parameters
  // fcr: density
  // ffd: diffusion coefficient
  // fcs: source
  // ffv: volume flux
  ConvDiffScal(
      double t, double dt, M& m, Par par, const FieldCell<Scal>* fcr,
      const FieldFace<Scal>* ffd, const FieldCell<Scal>* fcs,
      const FieldFace<Scal>* ffv)
      : UnsteadyIterativeSolver(t, dt)
      , m(m)
      , par(par)
      , fcr_(fcr)
      , ffd_(ffd)
      , fcs_(fcs)
      , ffv_(ffv) {}
  virtual const FieldCell<Scal>& GetField(Step) const = 0;
  virtual const FieldCell<Scal>& GetField() const {
    return GetField(Step::time_curr);
  }
  // Corrects field and exchanges halos.
  // l: layer to correct
  // uc: correction [i]
  virtual void CorrectField(Step l, const FieldCell<Scal>& uc) = 0;
  // Assembles linear system
  // fcu: field from previous iteration [a]
  // ffv: volume flux
  // Output:
  // fcucs_: linear system, overwritten, returned by GetDiag and GetConst
  virtual void Assemble(
      const FieldCell<Scal>& fcu, const FieldFace<Scal>& ffv) = 0;
  // Returns the diagonal coefficient of the equation
  virtual FieldCell<Scal> GetDiag() const = 0;
  // Returns the constant term of the equation
  virtual FieldCell<Scal> GetConst() const = 0;
  virtual const Par& GetPar() const {
    return par;
  }
  virtual void SetPar(Par par0) {
    par = par0;
  }

 protected:
  M& m;
  Par par; // parameters
  const FieldCell<Scal>* fcr_; // density
  const FieldFace<Scal>* ffd_; // diffusion
  const FieldCell<Scal>* fcs_; // source
  const FieldFace<Scal>* ffv_; // volume flux
};
