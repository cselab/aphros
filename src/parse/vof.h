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
  p->clipth = var.Double["clipth"];
  p->dim = var.Int["dim"];
  p->bcc_reflectpoly = var.Int["bcc_reflectpoly"];
  p->dumppolymarch_fill = var.Double["dumppolymarch_fill"];
  p->sharpen = var.Int["sharpen"];
  p->sharpen_cfl = var.Double["sharpen_cfl"];
  p->avgnorm0 = var.Double["avgnorm0"];
  p->avgnorm1 = var.Double["avgnorm1"];
  p->clfixed = var.Double["clfixed"];
  p->clfixed_x = Vect(var.Vect["clfixed_x"]);
  p->cloverride = var.Int["cloverride"];
  p->layers = var.Int["vofm_layers"];
  p->coalth = var.Double["vofm_coalth"];

  using Par = typename Vof::Par;
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
