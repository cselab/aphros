// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include "geom/mesh.h"
#include "solver/embed.h"
#include "solver/solver.h"

template <class Scal>
struct BCondAdvection {
  enum class Halo { fill, reflect };

  BCondAdvection() = default;
  BCondAdvection(size_t nci) : nci(nci) {}
  size_t GetNci() const {
    return nci;
  }

  size_t nci;
  Scal clear0 = 0; // snap to 0 if vf<clear0
  Scal clear1 = 1; // snap to 1 if vf>clear1
  Halo halo = Halo::reflect;
  Scal fill_vf; // volume fraction in halo cells if halo=fill
  Scal fill_cl; // color in halo cells if halo=fill
  Scal contang = -1; // contact angle [0..pi], negative to disable
};

template <class M_>
class AdvectionSolver : public UnsteadyIterativeSolver {
  using M = M_;
  using Scal = typename M::Scal;

 protected:
  M& m;
  const FieldEmbed<Scal>* fev_; // volume flux [velocity*area]
  const FieldCell<Scal>* fcs_; // source [value/time]

 public:
  // fev: volume flux
  // fcs: source
  AdvectionSolver(
      double t, double dt, M& m, const FieldEmbed<Scal>* fev,
      const FieldCell<Scal>* fcs)
      : UnsteadyIterativeSolver(t, dt), m(m), fev_(fev), fcs_(fcs) {}
  // Postprocessing after time step (dumps)
  virtual void PostStep() {}
  // Volume fraction
  virtual const FieldCell<Scal>& GetField(Step) const = 0;
  // Volume fraction at last time step
  virtual const FieldCell<Scal>& GetField() const {
    return GetField(Step::time_curr);
  }
  virtual void DumpInterface(std::string /*filename*/) const {};
  virtual void DumpInterfaceMarch(std::string /*filename*/) const {};
};
