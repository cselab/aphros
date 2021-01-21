// Created by Petr Karnakov on 23.10.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <memory>
#include <string>

#include "linear.h"

namespace linear {

template <class M>
class SolverHypre : public Solver<M> {
 public:
  using Base = Solver<M>;
  using Conf = typename Base::Conf;
  using Info = typename Base::Info;
  using Scal = typename M::Scal;
  using Expr = typename M::Expr;
  struct Extra { // extra config
    std::string solver = "pcg"; // name of the solver to use
    int print = 0; // print level, 0 for none
  };
  SolverHypre(const Conf& conf, const Extra& extra, const M&);
  ~SolverHypre();
  Info Solve(
      const FieldCell<Expr>& fc_system, const FieldCell<Scal>* fc_init,
      FieldCell<Scal>& fc_sol, M& m) override;

 private:
  struct Imp;
  const std::unique_ptr<Imp> imp;
};

} // namespace linear
