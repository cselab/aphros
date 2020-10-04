// Created by Petr Karnakov on 04.10.2020
// Copyright 2020 ETH Zurich

#pragma once

#include "linear.h"
#include "util/logger.h"

template <class M_>
auto ULinear<M_>::MakeLinearSolver(const Vars& var, std::string prefix)
    -> std::unique_ptr<linear::Solver<M>> {
  auto addprefix = [prefix](std::string name) {
    return "hypre_" + prefix + "_" + name;
  };

  const std::string name = var.String("linsolver_" + prefix, "hypre");
  typename linear::Solver<M>::Conf conf;
  conf.tol = var.Double[addprefix("tol")];
  conf.maxiter = var.Int[addprefix("maxiter")];
  if (name == "hypre") {
    typename linear::SolverHypre<M>::Extra extra;
    extra.solver = var.String[addprefix("solver")];
    extra.print = var.Int["hypre_print"];
    return std::make_unique<linear::SolverHypre<M>>(conf, extra);
  } else if (name == "conjugate") {
    return std::make_unique<linear::SolverConjugate<M>>(
        conf, typename linear::SolverConjugate<M>::Extra());
  } else if (name == "jacobi") {
    return std::make_unique<linear::SolverJacobi<M>>(
        conf, typename linear::SolverJacobi<M>::Extra());
  } else if (name == "default") {
    typename linear::SolverDefault<M>::Extra extra;
    using T = typename M::LS::T;
    extra.type = (prefix == "symm" ? T::symm : T::gen);
    extra.prefix = prefix;
    return std::make_unique<linear::SolverDefault<M>>(conf, extra);
  }
  fassert(false, "Unknown linsolver_" + prefix + "=" + name);
}
