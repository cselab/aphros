// Created by Petr Karnakov on 20.03.2021
// Copyright 2021 ETH Zurich

#pragma once

#include "hydro.h"

template <class M>
struct HydroPost {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  static FieldCell<Scal> GetField(
      const Hydro<M>* hydro, std::string name, const M& m);

  struct Imp;
};
