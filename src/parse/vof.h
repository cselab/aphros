#pragma once

#include <string>

#include "parse/vars.h"

template <class M, class Vof>
void Parse(typename Vof::Par* p, const Vars& var) {
  using Vect = typename M::Vect;
  p->verb = var.Int["vof_verb"];
  p->recolor_unionfind = var.Int["vof_recolor_unionfind"];
  p->recolor_reduce = var.Int["vof_recolor_reduce"];
  p->recolor_grid = var.Int["vof_recolor_grid"];
  p->vtkbin = var.Int["vtkbin"];
  p->vtkmerge = var.Int["vtkmerge"];
  p->vtkiso = var.Double["vtkiso"];
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
  p->dumppolymarch = var.Int["dumppolymarch"];
  p->dumppart = var.Int["dumppart"];
  p->dumppartinter = var.Int["dumppartinter"];
  p->bcc_reflect = var.Int["bcc_reflect"];
  p->bcc_reflectpoly = var.Int["bcc_reflectpoly"];
  p->dumppolymarch_fill = var.Double["dumppolymarch_fill"];
  p->bcc_fill = var.Double["bcc_fill"];
  p->bcc_clear = var.Int["bcc_clear"];
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
  p->part_maxr = var.Double["part_maxr"];
  p->sharpen = var.Int["sharpen"];
  p->sharpen_cfl = var.Double["sharpen_cfl"];
  p->avgnorm0 = var.Double["avgnorm0"];
  p->avgnorm1 = var.Double["avgnorm1"];
  p->clfixed = var.Double["clfixed"];
  p->clfixed_x = Vect(var.Vect["clfixed_x"]);
  p->layers = var.Int["vofm_layers"];
  p->coalth = var.Double["vofm_coalth"];

  using Par = typename Vof::Par;
  {
    using AF = typename Par::AF;
    std::string s = var.String["part_attrforce"];
    if (s == "line") {
      p->part_attrforce = AF::line;
    } else if (s == "center") {
      p->part_attrforce = AF::center;
    } else if (s == "volume") {
      p->part_attrforce = AF::volume;
    } else {
      throw std::runtime_error("Update: unknown part_attrforce=" + s);
    }
  }
  {
    using AR = typename Par::AR;
    std::string s = var.String["part_attrreconst"];
    if (s == "line") {
      p->part_attrreconst = AR::line;
    } else if (s == "volume") {
      p->part_attrreconst = AR::volume;
    } else {
      throw std::runtime_error("Update: unknown part_attrreconst=" + s);
    }
  }
  {
    using Scheme = typename Par::Scheme;
    std::string s = var.String["vof_scheme"];
    if (s == "plain") {
      p->scheme = Scheme::plain;
    } else if (s == "aulisa") {
      p->scheme = Scheme::aulisa;
    } else if (s == "weymouth") {
      p->scheme = Scheme::weymouth;
    } else {
      throw std::runtime_error("Update: unknown vof_scheme=" + s);
    }
  }
}
