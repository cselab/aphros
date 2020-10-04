// Created by Petr Karnakov on 04.10.2020
// Copyright 2020 ETH Zurich

#pragma once

#include "linear.h"
#include "util/logger.h"

template <class M_>
auto ULinear<M_>::GetLinearSolverHypre(const Vars& var, std::string prefix)
    -> std::unique_ptr<linear::Solver<M>> {
  fassert(
      prefix == "symm" || prefix == "gen" || prefix == "vort",
      "Unknown prefix=" + prefix);

  auto addprefix = [prefix](std::string name) {
    return "hypre_" + prefix + "_" + name;
  };

  typename linear::Solver<M>::Conf conf;
  conf.tol = var.Double[addprefix("tol")];
  conf.maxiter = var.Int[addprefix("maxiter")];
  typename linear::SolverHypre<M>::Extra extra;
  extra.solver = var.String[addprefix("solver")];
  extra.print = var.Int["hypre_print"];
  return std::make_unique<linear::SolverHypre<M>>(conf, extra);
}
