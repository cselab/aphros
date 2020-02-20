// Created by Petr Karnakov on 27.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <memory>

#include "advection.h"
#include "solver/cond.h"

template <class M_>
class Tvd final : public AdvectionSolver<M_> {
 public:
  using M = M_;
  using P = AdvectionSolver<M>;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

 public:
  struct Par {
    Scal sharp = 0.;
    Scal sharpo = 0.;
    Scal sharp_max = 1.;
    bool split = false;
  };
  // Constructor
  Tvd(M& m, const FieldCell<Scal>& fcu, const MapCondFaceAdvection<Scal>& mfc,
      const FieldEmbed<Scal>* fev, const FieldCell<Scal>* fcs, double t,
      double dt, Par par);
  ~Tvd();
  const Par& GetPar() const;
  void SetPar(Par);
  void StartStep() override;
  void MakeIteration() override;
  void FinishStep() override;
  const FieldCell<Scal>& GetField(Step l) const override;
  using P::GetField;

 private:
  struct Imp; // implementation
  std::unique_ptr<Imp> imp;

  using P::m;
};
