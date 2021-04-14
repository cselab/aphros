// Created by Petr Karnakov on 14.04.2021
// Copyright 2021 ETH Zurich

#pragma once

#include <memory>

#include "linear.h"

namespace linear {

template <class M>
class SolverConjugateCL : public Solver<M> {
 public:
  using Base = Solver<M>;
  using Conf = typename Base::Conf;
  using Info = typename Base::Info;
  using Scal = typename M::Scal;
  using Expr = typename M::Expr;
  struct Extra {
    bool residual_max = false; // if true, use max-norm of residual, else L2
  };
  SolverConjugateCL(const Conf& conf, const Extra& extra, const M&);
  ~SolverConjugateCL();
  Info Solve(
      const FieldCell<Expr>& fc_system, const FieldCell<Scal>* fc_init,
      FieldCell<Scal>& fc_sol, M& m) override;

 private:
  struct Imp;
  const std::unique_ptr<Imp> imp;
};

} // namespace linear
