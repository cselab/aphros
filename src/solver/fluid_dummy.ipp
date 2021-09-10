// Created by Petr Karnakov on 30.09.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <cmath>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>

#include "approx.h"
#include "approx_eb.h"
#include "debug/isnan.h"
#include "debug/linear.h"
#include "fluid.h"
#include "fluid_dummy.h"
#include "linear/linear.h"
#include "util/posthook.h"

template <class EB_>
struct FluidDummy<EB_>::Imp {
  using Owner = FluidDummy<EB>;
  using ExprFace = typename M::ExprFace;
  using Expr = typename M::Expr;
  using UEB = UEmbed<M>;

  Imp(Owner* owner, const EB& eb0, const FieldCell<Vect>& fcvel,
      const MapEmbed<BCondFluid<Vect>>& mebc, const Vars& var0)
      : owner_(owner), m(owner_->m), eb(eb0), var(var0) {
    fcvel_ = fcvel;
    fcp_.Reinit(m, 0);
    const auto ffwe = UEB::Interpolate(fcvel_, {}, eb);
    fev_.Reinit(m, 0);
    eb.LoopFaces([&](auto cf) { //
      fev_[cf] = ffwe[cf].dot(eb.GetSurface(cf));
    });
    mebc_vel_ = ConvertBCondFluidToVelocity<M>(mebc);
  }
  double GetAutoTimeStep() const {
    Scal dt = std::numeric_limits<Scal>::max();
    for (auto f : eb.Faces()) {
      const Scal vel = fev_[f] / m.GetArea(f);
      if (vel != 0.) {
        dt = std::min<Scal>(dt, std::abs(m.GetCellSize()[0] / vel));
      }
    }
    return dt;
  }
  void MakeIteration() {
    auto sem = m.GetSem("fluid-iter");
    const Scal t = owner_->GetTime();
    const Scal dt = owner_->GetTimeStep();
    auto& ffv = fev_.template Get<FieldFace<Scal>>();
    if (sem("local")) {
      FluidDummyHook(fcvel_, ffv, t, dt, var, m);
      owner_->IncIter();
    }
  }

  Owner* owner_;
  M& m;
  const EB& eb;
  const Vars& var;

  FieldEmbed<Scal> fev_; // volume flux
  FieldCell<Scal> fcp_; // pressure
  FieldCell<Vect> fcvel_; // velocity
  MapEmbed<BCond<Vect>> mebc_vel_;
};

template <class EB_>
FluidDummy<EB_>::FluidDummy(
    M& m_, const EB& eb, const FieldCell<Vect>& fcvel,
    const MapEmbed<BCondFluid<Vect>>& mebc, const FieldCell<Scal>* fcr,
    const FieldCell<Scal>* fcd, const FieldCell<Vect>* fcf,
    const FieldEmbed<Scal>* febp, const FieldCell<Scal>* fcsv,
    const FieldCell<Scal>* fcsm, double t, double dt, const Vars& var)
    : Base(t, dt, m_, fcr, fcd, fcf, febp, fcsv, fcsm)
    , imp(new Imp(this, eb, fcvel, mebc, var)) {}

template <class EB_>
FluidDummy<EB_>::~FluidDummy() = default;

template <class EB_>
void FluidDummy<EB_>::MakeIteration() {
  return imp->MakeIteration();
}

template <class EB_>
auto FluidDummy<EB_>::GetVelocity(Step) const -> const FieldCell<Vect>& {
  return imp->fcvel_;
}

template <class EB_>
auto FluidDummy<EB_>::GetPressure(Step) const -> const FieldCell<Scal>& {
  return imp->fcp_;
}

template <class EB_>
auto FluidDummy<EB_>::GetVolumeFlux(Step) const -> const FieldEmbed<Scal>& {
  return imp->fev_;
}

template <class EB_>
double FluidDummy<EB_>::GetAutoTimeStep() const {
  return imp->GetAutoTimeStep();
}

template <class EB_>
auto FluidDummy<EB_>::GetVelocityCond() const -> const MapEmbed<BCond<Vect>>& {
  return imp->mebc_vel_;
}
