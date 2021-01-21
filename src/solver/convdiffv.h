// Created by Petr Karnakov on 27.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include "convdiff.h"
#include "solver.h"

template <class EB_>
class ConvDiffVect : public UnsteadyIterativeSolver {
 public:
  using EB = EB_;
  using M = typename EB::M;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Par = typename ConvDiffScal<M>::Par;
  template <class T>
  using FieldFaceb = typename EmbedTraits<EB>::template FieldFaceb<T>;
  static constexpr size_t dim = M::dim;

  // par: parameters
  // fcr: density
  // ffd: diffusion coefficient
  // fcs: source
  // ffv: volume flux
  ConvDiffVect(
      double time, double dt, M& m_, const EB& eb_, Par par_,
      const FieldCell<Scal>* fcr, const FieldFaceb<Scal>* ffd,
      const FieldCell<Vect>* fcs, const FieldFaceb<Scal>* ffv)
      : UnsteadyIterativeSolver(time, dt)
      , m(m_)
      , eb(eb_)
      , par(par_)
      , fcr_(fcr)
      , ffd_(ffd)
      , fcs_(fcs)
      , ffv_(ffv) {}
  virtual const FieldCell<Vect>& GetVelocity(Step) const = 0;
  virtual const FieldCell<Vect>& GetVelocity() const {
    return GetVelocity(Step::time_curr);
  }
  virtual void CorrectVelocity(Step, const FieldCell<Vect>&) = 0;
  // Assembles linear system
  // fcw: current velocity
  // ffv: volume flux
  // Output:
  // linear system returned by GetDiag() and GetConst()
  virtual void Assemble(
      const FieldCell<Vect>& fcw, const FieldFaceb<Scal>& ffv) = 0;
  // Returns the diagonal coefficient of the equation in direction d
  virtual FieldCell<Scal> GetDiag(size_t d) const = 0;
  // Returns the constant term of the equation in direction d
  virtual FieldCell<Scal> GetConst(size_t d) const = 0;
  virtual const Par& GetPar() const {
    return par;
  }
  virtual void SetPar(Par par0) {
    par = par0;
  }

 protected:
  M& m;
  const EB& eb;
  Par par; // parameters
  const FieldCell<Scal>* fcr_; // density
  const FieldFaceb<Scal>* ffd_; // dynamic viscosity
  const FieldCell<Vect>* fcs_; // force
  const FieldFaceb<Scal>* ffv_; // volume flux
};
