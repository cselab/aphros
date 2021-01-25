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

// Converts vector conditions to scalar.
// mfv: vector velocity conditions
// d: direction, 0..2
template <class M>
MapEmbed<BCond<typename M::Scal>> GetScalarCond(
    const MapEmbed<BCond<typename M::Vect>>& mev, size_t d, const M& m) {
  using Scal = typename M::Scal;
  MapEmbed<BCond<Scal>> mes;

  for (auto& p : mev.GetMapFace()) {
    const IdxFace f = p.first;
    const auto& bcv = p.second;
    auto& bcs = mes[f];
    bcs.type = bcv.type;
    bcs.nci = bcv.nci;
    switch (bcv.type) {
      case BCondType::dirichlet:
      case BCondType::neumann:
        bcs.val = bcv.val[d];
        break;
      case BCondType::mixed:
      case BCondType::reflect:
        if (size_t(m.GetDir(f)) == d) {
          bcs.type = BCondType::dirichlet;
        } else {
          bcs.type = BCondType::neumann;
        }
        bcs.val = bcv.val[d];
        break;
      case BCondType::extrap:
        // nop
        break;
    }
  }
  for (auto& p : mev.GetMapCell()) {
    const IdxCell c = p.first;
    const auto& bcv = p.second;
    auto& bcs = mes[c];
    bcs.type = bcv.type;
    bcs.nci = bcv.nci;
    switch (bcv.type) {
      case BCondType::dirichlet:
      case BCondType::neumann:
        bcs.val = bcv.val[d];
        break;
      case BCondType::mixed:
      case BCondType::reflect:
        bcs.type = BCondType::neumann;
        // TODO revise to have zero normal component
      case BCondType::extrap:
        // nop
        break;
    }
  }
  return mes;
}
