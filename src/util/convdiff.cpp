// Created by Petr Karnakov on 20.11.2019
// Copyright 2019 ETH Zurich

#include "convdiff.h"

#include "solver/convdiffe.h"
#include "solver/convdiffi.h"
#include "solver/convdiffvg.h"

template <class M>
std::unique_ptr<ConvDiffVect<M>> GetConvDiff<M>::operator()(
    Conv conv, M& m, const FieldCell<Vect>& fcw, const MapCondFace& mfc,
    const FieldCell<Scal>* fcr, const FieldFace<Scal>* ffd,
    const FieldCell<Vect>* fcs, const FieldFace<Scal>* ffv, double t, double dt,
    Par par) {
  using CDI = ConvDiffVectGeneric<M, ConvDiffScalImp<M>>; // implicit
  using CDE = ConvDiffVectGeneric<M, ConvDiffScalExp<M>>; // explicit

  switch (conv) {
    case Conv::imp:
      return std::unique_ptr<CDI>(
          new CDI(m, m, fcw, mfc, fcr, ffd, fcs, ffv, t, dt, par));
    case Conv::exp:
      return std::unique_ptr<CDE>(
          new CDE(m, m, fcw, mfc, fcr, ffd, fcs, ffv, t, dt, par));
  }
  return nullptr;
}

using M = MeshStructured<double, 3>;

template std::unique_ptr<ConvDiffVect<M>> GetConvDiff<M>::operator()(
    Conv conv, M& m, const FieldCell<Vect>& fcw, const MapCondFace& mfc,
    const FieldCell<Scal>* fcr, const FieldFace<Scal>* ffd,
    const FieldCell<Vect>* fcs, const FieldFace<Scal>* ffv, double t, double dt,
    Par par);
