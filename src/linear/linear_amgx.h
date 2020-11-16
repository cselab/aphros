// Created by Petr Karnakov on 15.11.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <memory>
#include <string>

#include "linear.h"

namespace linear {

template <class M>
class SolverAmgx : public Solver<M> {
 public:
  using Base = Solver<M>;
  using Conf = typename Base::Conf;
  using Info = typename Base::Info;
  using Scal = typename M::Scal;
  using Expr = typename M::Expr;
  struct Extra { // extra config
    std::string config_path = "amgx.json"; // path to json config file
    std::string config_extra = ""; // extra config as string
    std::string mode = "dDDI";
    std::string log_path = "amgx.log"; // path to log file
  };
  SolverAmgx(const Conf& conf, const Extra& extra, const M& m);
  ~SolverAmgx();
  Info Solve(
      const FieldCell<Expr>& fc_system, const FieldCell<Scal>* fc_init,
      FieldCell<Scal>& fc_sol, M& m) override;

 private:
  struct Imp;
  const std::unique_ptr<Imp> imp;
};

} // namespace linear
