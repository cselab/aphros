#pragma once

#include <string>

#include "parse/vars.h"
#include "solver/vof.h"

template <class M>
void Parse(typename solver::Vof<M>::Par* p, const Vars& var) {
  using Par = typename solver::Vof<M>::Par;
  p->curvgrad = var.Int["curvgrad"];
  p->part = var.Int["part"];
  p->part_verb = var.Int["part_verb"];
  p->part_relax = var.Double["part_relax"];
  p->part_h = var.Double["part_h"];
  p->part_kstr = var.Double["part_kstr"];
  p->part_kattr = var.Double["part_kattr"];
  p->part_kbend = var.Double["part_kbend"];
  p->part_bendmean = var.Int["part_bendmean"];
  p->part_dump_fr = var.Int["part_dump_fr"];
  p->part_report_fr = var.Int["part_report_fr"];
  p->part_k = var.Int["part_k"];
  p->part_intth = var.Double["part_intth"];
  p->poly_intth = var.Double["poly_intth"];
  p->clipth = var.Double["clipth"];
  p->dim = var.Int["dim"];
  p->dumppoly = var.Int["dumppoly"];
  p->dumppart = var.Int["dumppart"];
  p->dumppartinter = var.Int["dumppartinter"];
  p->bcc_reflect = var.Int["bcc_reflect"];
  p->part_constr = var.Int["part_constr"];
  p->part_segcirc = var.Double["part_segcirc"];
  p->part_itermax = var.Int["part_itermax"];
  p->part_tol = var.Double["part_tol"];
  p->part_np = var.Int["part_np"];
  p->part_ns = var.Int["part_ns"];
  p->part_tmax = var.Double["part_tmax"];
  p->part_dtmax = var.Double["part_dtmax"];
  p->part_anglim = var.Double["part_anglim"];
  p->part_dn = var.Int["part_dn"];

  {
    std::string s = var.String["part_attrforce"];
    if (s == "line") {
      p->part_attrforce = Par::AF::line;
    } else if (s == "center") {
      p->part_attrforce = Par::AF::center;
    } else if (s == "volume") {
      p->part_attrforce = Par::AF::volume;
    } else {
      throw std::runtime_error("Update: unknown part_attrforce=" + s);
    }
  }
  {
    std::string s = var.String["part_attrreconst"];
    if (s == "line") {
      p->part_attrreconst = Par::AR::line;
    } else if (s == "volume") {
      p->part_attrreconst = Par::AR::volume;
    } else {
      throw std::runtime_error("Update: unknown part_attrreconst=" + s);
    }
  }
}
