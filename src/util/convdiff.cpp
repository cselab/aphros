// Created by Petr Karnakov on 20.11.2019
// Copyright 2019 ETH Zurich

#include "convdiff.h"

#include "solver/convdiffe.h"
#include "solver/convdiffi.h"
#include "solver/convdiffvg.h"
#include "solver/embed.h"

template <class MEB>
std::unique_ptr<ConvDiffVect<MEB>> GetConvDiff<MEB>::operator()(
    Conv conv, typename MEB::M& m, const MEB& eb,
    const ConvDiffVectArgs<MEB>& args) {
  using CDI = ConvDiffVectGeneric<MEB, ConvDiffScalImp<MEB>>; // implicit
  using CDE = ConvDiffVectGeneric<MEB, ConvDiffScalExp<MEB>>; // explicit

  switch (conv) {
    case Conv::imp:
      return std::make_unique<CDI>(m, eb, args);
    case Conv::exp:
      return std::make_unique<CDE>(m, eb, args);
  }
  fassert(false);
}

using M = MeshStructured<double, 3>;
using EB = Embed<M>;

template struct GetConvDiff<M>;
template struct GetConvDiff<EB>;
