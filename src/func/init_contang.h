// Created by Petr Karnakov on 04.05.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <cmath>
#include <functional>
#include <stdexcept>

#include "geom/field.h"
#include "parse/vars.h"
#include "util/module.h"

template <class M>
class ModuleInitContang : public Module<ModuleInitContang<M>> {
 public:
  using Scal = typename M::Scal;
  using Module<ModuleInitContang>::Module;
  virtual void operator()(
      FieldCell<Scal>& fc_contang, const Vars& var, const M& m) = 0;
};
