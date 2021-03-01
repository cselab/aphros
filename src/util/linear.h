// Created by Petr Karnakov on 04.10.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <memory>

#include "linear/linear.h"
#include "parse/vars.h"

template <class M_>
class ULinear {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  // Constructs a linear solver from configuration (defaults defined in
  // base.conf) prefix: prefix for system type (symm, gen, vort)
  static std::unique_ptr<linear::Solver<M>> MakeLinearSolver(
      const Vars& var, std::string prefix, const M&);
};
