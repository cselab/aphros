// Created by Petr Karnakov on 10.04.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <memory>

#include "geom/map.h"
#include "solver/cond.h"
#include "solver/convdiff.h"
#include "solver/convdiffv.h"

template <class T>
MapEmbed<BCond<T>> GetBCond(const MapCondFace& mff) {
  MapEmbed<BCond<T>> mebc;
  for (auto& p : mff) {
    const IdxFace f = p.first;
    const auto& cb = p.second;

    auto& bc = mebc[f];
    bc.nci = cb->GetNci();

    if (auto cd = cb.template Get<CondFaceValFixed<T>>()) {
      bc.type = BCondType::dirichlet;
      bc.val = cd->second();
    } else if (auto cd = cb.template Get<CondFaceGradFixed<T>>()) {
      bc.type = BCondType::neumann;
      bc.val = cd->GetGrad();
    } else if (auto cd = cb.template Get<CondFaceReflect>()) {
      bc.type = BCondType::reflect;
    } else if (auto cd = cb.template Get<CondFaceExtrap>()) {
      bc.type = BCondType::extrap;
    } else {
      throw std::runtime_error("GetBCond: unknown condition");
    }
  }
  return mebc;
}

template <class T>
MapCondFace GetCond(const MapEmbed<BCond<T>>& mebc) {
  MapCondFace mf;
  for (auto& p : mebc.GetMapFace()) {
    const IdxFace f = p.first;
    auto& bc = p.second;
    auto nci = bc.nci;
    auto& cond = mf[f];
    switch (bc.type) {
      case BCondType::dirichlet:
        cond.Set<CondFaceValFixed<T>>(bc.val, nci);
        break;
      case BCondType::neumann:
        cond.Set<CondFaceGradFixed<T>>(bc.val, nci);
        break;
      case BCondType::reflect:
        cond.Set<CondFaceReflect>(nci);
        break;
      case BCondType::extrap:
        cond.Set<CondFaceExtrap>(nci);
        break;
      default:
        throw std::runtime_error("GetCond: unknown condition");
    }
  }
  return mf;
}


// Converts vector conditions to scalar.
// mfv: vector velocity conditions
// d: direction, 0..2
template <class M>
MapCondFace GetScalarCond(const MapCondFace& mfv, size_t d, const M& m) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  MapCondFace mfs;
  // Face conditions for each velocity component
  for (auto& it : mfv) {
    IdxFace f = it.first;
    auto& cb = it.second;
    if (cb.template Get<CondFaceVal<Vect>>()) {
      mfs[f] = EvalComp<Vect>(cb, d);
    } else if (cb.template Get<CondFaceGrad<Vect>>()) {
      mfs[f] = EvalComp<Vect>(cb, d);
    } else if (cb.template Get<CondFaceReflect>()) {
      auto nci = cb->GetNci();
      // XXX: adhoc for cartesian grid
      if (d == m.GetNormal(f).abs().argmax()) {
        // normal, zero value
        mfs[f] = UniquePtr<CondFaceValFixed<Scal>>(0., nci);
      } else {
        // tangential, zero gradient
        mfs[f] = UniquePtr<CondFaceGradFixed<Scal>>(0., nci);
      }
    } else {
      throw std::runtime_error("GetScalarCond: unknown face condition");
    }
  }
  return mfs;
}

enum class Conv { exp, imp };

template <class M>
struct GetConvDiff {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Par = typename ConvDiffScal<M>::Par;

  std::unique_ptr<ConvDiffVect<M>> operator()(
      Conv conv, M& m, const M& eb, const FieldCell<Vect>& fcvel,
      const MapEmbed<BCond<Vect>>& mebc, const FieldCell<Scal>* fcr,
      const FieldFace<Scal>* ffd, const FieldCell<Vect>* fcs,
      const FieldEmbed<Scal>* fev, double t, double dt, Par par);
};

template <class M>
struct GetConvDiff<Embed<M>> {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Par = typename ConvDiffScal<M>::Par;

  std::unique_ptr<ConvDiffVect<Embed<M>>> operator()(
      Conv conv, M& m, const Embed<M>& eb, const FieldCell<Vect>& fcvel,
      const MapEmbed<BCond<Vect>>& mebc, const FieldCell<Scal>* fcr,
      const FieldEmbed<Scal>* fed, const FieldCell<Vect>* fcs,
      const FieldEmbed<Scal>* fev, double t, double dt, Par par);
};

template <class ConvDiffPar, class FluidPar>
void SetConvDiffPar(ConvDiffPar& d, const FluidPar& p) {
  d.relax = p.vrelax;
  d.second = p.second;
  d.sc = p.convsc;
  d.df = p.convdf;
  d.linreport = p.linreport;
  d.stokes = p.stokes;
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
        throw std::runtime_error("GetCond: not implemented");
      case BCondType::extrap:
        // nop
        break;
    }
  }
  return mes;
}
