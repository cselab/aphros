#pragma once

#include <stdexcept>
#include <string>

#include "parse/solver.h"
#include "solver/simple.h"

template <class M>
struct ParsePar<Simple<M>> {
  using Par = typename Simple<M>::Par;
  Par operator()(const Vars& var) {
    Par p;
    p.prelax = var.Double["prelax"];
    p.vrelax = var.Double["vrelax"];
    p.rhie = var.Double["rhie"];
    p.second = var.Int["second_order"];
    p.simpler = var.Int["simpler"];
    // TODO: add check for inletflux_numid
    p.inletflux_numid = var.Int["inletflux_numid"];
    p.convsc = GetConvSc(var.String["convsc"]);
    p.convdf = var.Double["convdf"];
    p.linreport = var.Int["linreport"];
    std::string conv = var.String["conv"];
    if (conv == "imp") {
      p.conv = Conv::imp;
    } else if (conv == "exp") {
      p.conv = Conv::exp;
    } else {
      throw std::runtime_error("Parse: unknown conv=" + conv);
    }
    return p;
  }
};
