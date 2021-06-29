// Created by Petr Karnakov on 25.12.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <string>

#include "parse/solver.h"
#include "solver/curv.h"
#include "solver/partstr.h"
#include "solver/partstrmeshm.h"

template <class Scal>
struct ParsePar<PartStr<Scal>> {
  using Par = typename PartStr<Scal>::Par;
  Par operator()(Scal hc, const Vars& var) {
    Par p;
    p.leq = var.Double["part_h"];
    p.relax = var.Double["part_relax"];
    p.npmax = var.Int["part_np"];
    p.segcirc = var.Double["part_segcirc"];
    p.hc = hc;
    p.dn = var.Int["part_dn"];
    return p;
  }
};

template <class M>
struct ParsePar<PartStrMeshM<M>> {
  using Scal = typename M::Scal;
  using Par = typename PartStrMeshM<M>::Par;
  Par operator()(typename PartStr<Scal>::Par ps, const Vars& var) {
    Par p;
    p.ps = ps;
    p.ns = var.Int["part_ns"];
    p.tol = var.Double["part_tol"];
    p.itermax = var.Int["part_itermax"];
    p.verb = var.Int["part_verb"];
    p.dim = var.Int["dim"];
    p.dump_fr = var.Int["part_dump_fr"];
    p.vtkbin = var.Int["vtkbin"];
    p.vtkmerge = var.Int["vtkmerge"];
    return p;
  }
};
