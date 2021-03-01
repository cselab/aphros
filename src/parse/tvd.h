// Created by Petr Karnakov on 04.04.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <string>

#include <solver/tvd.h>
#include "parse/solver.h"

template <class M>
struct ParsePar<Tvd<M>> {
  using Par = typename Tvd<M>::Par;
  Par operator()(const Vars& var) {
    Par p;
    p.sharp = var.Double["sharp"];
    p.sharpo = var.Double["sharpo"];
    p.split = var.Int["split"];
    return p;
  }
};
