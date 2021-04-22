// Created by Petr Karnakov on 10.04.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <memory>

#include "geom/map.h"
#include "solver/cond.h"
#include "solver/convdiff.h"
#include "solver/convdiffvg.h"

enum class Conv { exp, imp };

template <class MEB>
struct GetConvDiff {
  using Scal = typename MEB::Scal;
  using Vect = typename MEB::Vect;
  using Par = typename ConvDiffScal<MEB>::Par;

  std::unique_ptr<ConvDiffVect<MEB>> operator()(
      Conv conv, typename MEB::M& m, const MEB& eb,
      const ConvDiffVectArgs<MEB>& args);
};

template <class ConvDiffPar, class FluidPar>
void SetConvDiffPar(ConvDiffPar& d, const FluidPar& p) {
  d.relax = p.vrelax;
  d.second = p.second;
  d.sc = p.convsc;
  d.df = p.convdf;
  d.stokes = p.stokes;
  d.explconv = p.explconv;
  d.symm = p.convsymm;
}

template <class ConvDiffPar, class FluidPar>
ConvDiffPar UpdateConvDiffPar(ConvDiffPar d, const FluidPar& p) {
  SetConvDiffPar(d, p);
  return d;
}
