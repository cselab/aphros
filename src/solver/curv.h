// Created by Petr Karnakov on 25.12.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <memory>

#include "advection.h"
#include "geom/mesh.h"
#include "partstrmeshm.h"

template <class M_>
struct UCurv {
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Plic = generic::Plic<Scal>;

  // Computes curvature with height functions.
  // fcu: volume fraction [a]
  // fcn: normal
  // edim: effective dimension (2 or 3)
  // Output:
  // fck: curvature
  static void CalcCurvHeight(
      const FieldCell<Scal>& fcu, const FieldCell<Vect>& fcn, size_t edim,
      FieldCell<Scal>& fck, M& m);

  template <class EB>
  static std::unique_ptr<PartStrMeshM<M>> CalcCurvPart(
      const Plic& plic, const typename PartStrMeshM<M>::Par& par,
      const FieldCell<Scal>* fc_contang, const Multi<FieldCell<Scal>*>& fck,
      M& m, const EB& eb);

  static std::unique_ptr<PartStrMeshM<M>> CalcCurvPart(
      const AdvectionSolver<M>* asbase,
      const typename PartStrMeshM<M>::Par& par,
      const FieldCell<Scal>* fc_contang, const Multi<FieldCell<Scal>*>& fck,
      M& m);

 private:
  struct Imp;
};
