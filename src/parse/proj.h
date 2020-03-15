// Created by Petr Karnakov on 30.09.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <stdexcept>
#include <string>

#include "parse/vars.h"
#include "solver/proj.h"

template <class M>
struct ParsePar<Proj<M>> {
  using Par = typename Proj<M>::Par;
  Par operator()(const Vars& var) {
    Par p;
    p.prelax = var.Double["prelax"];
    p.vrelax = var.Double["vrelax"];
    p.second = var.Int["second_order"];
    // TODO: add check for inletflux_numid
    p.inletflux_numid = var.Int["inletflux_numid"];
    p.convsc = GetConvSc(var.String["convsc"]);
    p.convdf = var.Double["convdf"];
    p.stokes = var.Int["stokes"];
    p.convsymm = var.Int["convsymm"];
    p.explconv = var.Int["explconv"];
    p.bcg = var.Int["proj_bcg"];
    std::string conv = var.String["conv"];
    if (conv == "imp") {
      p.conv = Conv::imp;
    } else if (conv == "exp") {
      p.conv = Conv::exp;
    } else {
      throw std::runtime_error("Parse: unknown conv=" + conv);
    }
    p.outlet_relax = var.Double["outlet_relax"];
    return p;
  }
};
