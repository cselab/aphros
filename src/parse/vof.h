// Created by Petr Karnakov on 04.04.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <string>

#include "parse/solver.h"
#include "solver/vof.h"
#include "solver/vofm.h"

template <class M>
struct ParsePar<Vof<M>> {
  using Par = typename Vof<M>::Par;
  using Vect = typename M::Vect;
  Par operator()(const Vars& var) {
    Par p;
    p.verb = var.Int["vof_verb"];
    p.recolor_unionfind = var.Int["vof_recolor_unionfind"];
    p.recolor_reduce = var.Int["vof_recolor_reduce"];
    p.recolor_grid = var.Int["vof_recolor_grid"];
    p.vtkbin = var.Int["vtkbin"];
    p.vtkmerge = var.Int["vtkmerge"];
    p.vtkiso = var.Double["vtkiso"];
    p.clipth = var.Double["clipth"];
    p.filterth = var.Double["filterth"];
    p.dim = var.Int["dim"];
    p.bcc_reflectpoly = var.Int["bcc_reflectpoly"];
    p.dumppolymarch_fill = var.Double["dumppolymarch_fill"];
    p.sharpen = var.Int["sharpen"];
    p.sharpen_cfl = var.Double["sharpen_cfl"];
    p.avgnorm0 = var.Double["avgnorm0"];
    p.avgnorm1 = var.Double["avgnorm1"];
    p.clfixed = var.Double["clfixed"];
    p.clfixed_x = Vect(var.Vect["clfixed_x"]);
    p.cloverride = var.Int["cloverride"];
    p.layers = var.Int["vofm_layers"];
    p.coalth = var.Double["vofm_coalth"];
    p.extrapolate_boundaries = var.Int["vof_extrapolate_boundaries"];

    using Scheme = typename Par::Scheme;
    std::string s = var.String["vof_scheme"];
    if (s == "plain") {
      p.scheme = Scheme::plain;
    } else if (s == "aulisa") {
      p.scheme = Scheme::aulisa;
    } else if (s == "weymouth") {
      p.scheme = Scheme::weymouth;
    } else {
      fassert(false, "Update: unknown vof_scheme=" + s);
    }
    return p;
  }
};

template <class EB>
struct ParsePar<Vofm<EB>> {
  using M = typename EB::M;
  using Par = typename Vofm<M>::Par;
  Par operator()(const Vars& var) {
    return ParsePar<Vof<M>>()(var);
  }
};
