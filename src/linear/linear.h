// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

#include "debug/isnan.h"
#include "debug/linear.h"
#include "geom/mesh.h"

namespace linear {

template <class M>
class Solver {
 public:
  using Scal = typename M::Scal;
  using Expr = typename M::Expr;

  struct Conf {
    Scal tol = 0;
    int maxiter = 100;
  };

  struct Info {
    Scal residual;
    int iter;
  };

  Solver(const Conf& conf_) : conf(conf_) {}
  virtual ~Solver() = default;
  // Solves linear system
  //   system(x) = 0
  //
  // Input:
  // fc_system: coefficients `e` of linear expressions
  //   e[0] * x[c] + sum_i(e[i + 1] * x[c(i)]) + e.back() = 0
  // fc_init: initial guess, may equal &fc_sol, nullptr to assume zero guess.
  //
  // Output:
  // fc_sol: solution
  // Info: final residual and number of iterations
  virtual Info Solve(
      const FieldCell<Expr>& fc_system, const FieldCell<Scal>* fc_init,
      FieldCell<Scal>& fc_sol, M& m) = 0;
  virtual void SetConf(const Conf& c) {
    conf = c;
  }
  virtual const Conf& GetConf() {
    return conf;
  }

 protected:
  Conf conf;
};

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
  SolverHypre(const Conf& conf, const Extra& extra);
  ~SolverHypre();
  Info Solve(
      const FieldCell<Expr>& fc_system, const FieldCell<Scal>* fc_init,
      FieldCell<Scal>& fc_sol, M& m) override;

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};

template <class M>
class SolverConjugate : public Solver<M> {
 public:
  using Base = Solver<M>;
  using Conf = typename Base::Conf;
  using Info = typename Base::Info;
  using Scal = typename M::Scal;
  using Expr = typename M::Expr;
  struct Extra {};
  SolverConjugate(const Conf& conf, const Extra& extra);
  ~SolverConjugate();
  Info Solve(
      const FieldCell<Expr>& fc_system, const FieldCell<Scal>* fc_init,
      FieldCell<Scal>& fc_sol, M& m) override;

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};

template <class M>
class SolverJacobi : public Solver<M> {
 public:
  using Base = Solver<M>;
  using Conf = typename Base::Conf;
  using Info = typename Base::Info;
  using Scal = typename M::Scal;
  using Expr = typename M::Expr;
  struct Extra {};
  SolverJacobi(const Conf& conf, const Extra& extra);
  ~SolverJacobi();
  Info Solve(
      const FieldCell<Expr>& fc_system, const FieldCell<Scal>* fc_init,
      FieldCell<Scal>& fc_sol, M& m) override;

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};

} // namespace linear

