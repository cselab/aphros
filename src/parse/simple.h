// Created by Petr Karnakov on 04.04.2019
// Copyright 2019 ETH Zurich

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
    // TODO: add check for inletflux_numid
    p.inletflux_numid = var.Int["inletflux_numid"];
    p.convsc = GetConvSc(var.String["convsc"]);
    p.convdf = var.Double["convdf"];
    p.stokes = var.Int["stokes"];
    p.convsymm = var.Int["convsymm"];
    p.explconv = var.Int["explconv"];
    p.explviscous = var.Int["explviscous"];
    std::string conv = var.String["conv"];
    if (conv == "imp") {
      p.conv = Conv::imp;
    } else if (conv == "exp") {
      p.conv = Conv::exp;
    } else {
      fassert(false, "Parse: unknown conv=" + conv);
    }
    p.outlet_relax = var.Double["outlet_relax"];
    return p;
  }
};
