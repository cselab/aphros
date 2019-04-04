#pragma once

#include <string>

#include "parse/vars.h"
#include "solver/tvd.h"

template <class M>
void Parse(typename solver::Tvd<M>::Par* p, const Vars& var) {
  p->sharp = var.Double["sharp"];
  p->sharpo = var.Double["sharpo"];
  p->split = var.Int["split"];
}
