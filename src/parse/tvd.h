#pragma once

#include <string>

#include "parse/solver.h"
#include "solver/tvd.h"

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
