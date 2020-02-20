// Created by Petr Karnakov on 20.11.2019
// Copyright 2019 ETH Zurich

#include "convdiff.h"

#include "solver/convdiffe.h"
#include "solver/convdiffi.h"
#include "solver/convdiffvg.h"
#include "solver/embed.h"

template <class M>
std::unique_ptr<ConvDiffVect<M>> GetConvDiff<M>::operator()(
    Conv conv, M& m, const M& eb, const FieldCell<Vect>& fcw,
    const MapCondFace& mfc, const FieldCell<Scal>* fcr,
    const FieldFace<Scal>* ffd, const FieldCell<Vect>* fcs,
    const FieldEmbed<Scal>* fev, double t, double dt, Par par) {
  using CDI = ConvDiffVectGeneric<M, ConvDiffScalImp<M>>; // implicit
  using CDE = ConvDiffVectGeneric<M, ConvDiffScalExp<M>>; // explicit

  switch (conv) {
    case Conv::imp:
      return std::unique_ptr<CDI>(new CDI(
          m, eb, fcw, mfc, fcr, ffd, fcs, &fev->GetFieldFace(), t, dt, par));
    case Conv::exp:
      return std::unique_ptr<CDE>(new CDE(
          m, eb, fcw, mfc, fcr, ffd, fcs, &fev->GetFieldFace(), t, dt, par));
  }
  return nullptr;
}

template <class M>
std::unique_ptr<ConvDiffVect<Embed<M>>> GetConvDiff<Embed<M>>::operator()(
    Conv conv, M& m, const Embed<M>& eb, const FieldCell<Vect>& fcw,
    const MapCondFace& mfc, const FieldCell<Scal>* fcr,
    const FieldEmbed<Scal>* fed, const FieldCell<Vect>* fcs,
    const FieldEmbed<Scal>* fev, double t, double dt, Par par) {
  using EB = Embed<M>;
  using CDE = ConvDiffVectGeneric<EB, ConvDiffScalExp<EB>>; // explicit

  switch (conv) {
    case Conv::exp:
      return std::unique_ptr<CDE>(
          new CDE(m, eb, fcw, mfc, fcr, fed, fcs, fev, t, dt, par));
  }
  return nullptr;
}

using M = MeshStructured<double, 3>;
using EB = Embed<M>;

template struct GetConvDiff<M>;
template struct GetConvDiff<EB>;
