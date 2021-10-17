// Created by Petr Karnakov on 04.10.2020
// Copyright 2020 ETH Zurich

#pragma once

#include "linear.h"
#include "logger.h"

template <class M_>
auto ULinear<M_>::MakeLinearSolver(
    const Vars& var, std::string prefix, const M& m)
    -> std::unique_ptr<linear::Solver<M>> {
#if USEFLAG(HYPRE)
  FORCE_LINK(linear_hypre);
#endif
#if USEFLAG(AMGX)
  FORCE_LINK(linear_amgx);
#endif
  FORCE_LINK(linear_conjugate);
  FORCE_LINK(linear_jacobi);

  auto addprefix = [prefix](std::string name) {
    return "hypre_" + prefix + "_" + name;
  };

  const std::string name = var.String("linsolver_" + prefix, "hypre");
  if (auto* mod = linear::ModuleLinear<M>::GetInstance(name)) {
    return mod->Make(var, prefix, m);
  }
  fassert(false, "Unknown linsolver_" + prefix + "=" + name);
}
