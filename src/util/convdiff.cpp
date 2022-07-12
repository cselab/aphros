// Created by Petr Karnakov on 20.11.2019
// Copyright 2019 ETH Zurich

#include "convdiff.h"

#include "solver/convdiffe.h"
#include "solver/convdiffi.h"
#include "solver/convdiffvg.h"
#include "solver/embed.h"
#include "util/make_unique.h"

template <class MEB>
std::unique_ptr<ConvDiffVect<MEB>> GetConvDiff<MEB>::operator()(
    Conv conv, typename MEB::M& m, const MEB& eb,
    const ConvDiffVectArgs<MEB>& args) {
  using CDI = ConvDiffVectGeneric<MEB, ConvDiffScalImp<MEB>>; // implicit
  using CDE = ConvDiffVectGeneric<MEB, ConvDiffScalExp<MEB>>; // explicit

  switch (conv) {
    case Conv::imp:
      return MakeUnique<CDI>(m, eb, args);
    case Conv::exp:
      return MakeUnique<CDE>(m, eb, args);
  }
  fassert(false);
}

#define X(dim) template struct GetConvDiff<MeshCartesian<double, dim>>;
MULTIDIMX
#undef X

#define X(dim) template struct GetConvDiff<Embed<MeshCartesian<double, dim>>>;
MULTIDIMX
#undef X
