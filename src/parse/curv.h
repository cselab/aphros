#pragma once

#include <string>

#include "parse/vars.h"
#include "solver/partstrmeshm.h"
#include "solver/partstr.h"

using namespace solver;

template <class Scal>
void Parse(typename PartStr<Scal>::Par& p, Scal hc, const Vars& var) {
  p.leq = var.Double["part_h"];
  p.relax = var.Double["part_relax"];
  p.npmax = var.Int["part_np"];
  p.segcirc = var.Double["part_segcirc"];
  p.hc = hc;
  p.dn = var.Int["part_dn"];
}

template <class M, class Scal = typename M::Scal>
void Parse(
    typename PartStrMeshM<M>::Par& p, const typename PartStr<Scal>::Par& ps,
    const Vars& var) {
  p.ps = std::make_shared<typename PartStr<Scal>::Par>(ps);
  p.dump_fr = var.Int["part_dump_fr"];
  p.ns = var.Int["part_ns"];
  p.tol = var.Double["part_tol"];
  p.itermax = var.Int["part_itermax"];
  p.verb = var.Int["part_verb"];
  p.dim = var.Int["dim"];
  p.vtkbin = var.Int["vtkbin"];
  p.vtkmerge = var.Int["vtkmerge"];
}
