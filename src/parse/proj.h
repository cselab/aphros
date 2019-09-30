#pragma once

#include <string>
#include <stdexcept>

#include "parse/vars.h"
#include "solver/proj.h"

template <class M>
void Parse(typename solver::Proj<M>::Par* p, const Vars& var) {
  p->prelax = var.Double["prelax"];
  p->vrelax = var.Double["vrelax"];
  p->second = var.Int["second_order"];
  p->guessextra = var.Double["guessextra"];
  // TODO: add check for inletflux_numid
  p->inletflux_numid = var.Int["inletflux_numid"];
  p->convsc = solver::GetConvSc(var.String["convsc"]);
  p->convdf = var.Double["convdf"];
  p->linreport = var.Int["linreport"];
  using Conv = typename solver::Proj<M>::Conv;
  std::string conv = var.String["conv"];
  if (conv == "imp") {
    p->conv = Conv::imp;
  } else if (conv == "exp") {
    p->conv = Conv::exp;
  } else {
    throw std::runtime_error("Parse: unknown conv=" + conv);
  }
}
