// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <memory>

#include "debug/linear.h"
#include "geom/mesh.h"
#include "parse/vars.h"
#include "util/module.h"

namespace linear {

template <class M>
class Solver {
 public:
  using Scal = typename M::Scal;
  using Expr = typename M::Expr;

  struct Conf {
    Scal tol = 0;
    int miniter = 1;
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
class ModuleLinear : public Module<ModuleLinear<M>> {
 public:
  using Vect = typename M::Vect;
  using Module<ModuleLinear>::Module;
  virtual std::unique_ptr<Solver<M>> Make(
      const Vars&, std::string prefix, const M& m) = 0;
  static typename Solver<M>::Conf GetConf(const Vars& var, std::string prefix) {
    auto addprefix = [prefix](std::string name) {
      return "hypre_" + prefix + "_" + name;
    };
    typename linear::Solver<M>::Conf conf;
    conf.tol = var.Double[addprefix("tol")];
    conf.maxiter = var.Int[addprefix("maxiter")];
    conf.miniter = var.Int(addprefix("miniter"), 0);
    return conf;
  }
};

template <class M>
class SolverConjugate : public Solver<M> {
 public:
  using Base = Solver<M>;
  using Conf = typename Base::Conf;
  using Info = typename Base::Info;
  using Scal = typename M::Scal;
  using Expr = typename M::Expr;
  struct Extra {
    bool residual_max = false; // if true, use max-norm of residual, else L2
  };
  SolverConjugate(const Conf& conf, const Extra& extra, const M&);
  ~SolverConjugate();
  Info Solve(
      const FieldCell<Expr>& fc_system, const FieldCell<Scal>* fc_init,
      FieldCell<Scal>& fc_sol, M& m) override;

 private:
  struct Imp;
  const std::unique_ptr<Imp> imp;
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
  SolverJacobi(const Conf& conf, const Extra& extra, const M&);
  ~SolverJacobi();
  Info Solve(
      const FieldCell<Expr>& fc_system, const FieldCell<Scal>* fc_init,
      FieldCell<Scal>& fc_sol, M& m) override;

 private:
  struct Imp;
  const std::unique_ptr<Imp> imp;
};

} // namespace linear
