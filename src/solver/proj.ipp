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
#include "linear/linear.h"
#include "proj.h"
#include "util/convdiff.h"
#include "util/fluid.h"
#include "util/metrics.h"

// ranges (cells/faces)
// [i]: inner
// [s]: support
// [a]: all

// fields
// p: pressure
// gp: pressure gradient
// w: velocity
// v: volume flux
// we: predicted velocity (after solving velocity equations)
// ve: predicted volume flux

template <class EB_>
struct Proj<EB_>::Imp {
  using Owner = Proj<EB>;
  using ExprFace = typename M::ExprFace;
  using Expr = typename M::Expr;
  using UEB = UEmbed<M>;

  Imp(Owner* owner, const EB& eb0, const FieldCell<Vect>& fcvel,
      const MapEmbed<BCondFluid<Vect>>& mebc,
      const MapCell<std::shared_ptr<CondCellFluid>>& mcc, Par par)
      : owner_(owner)
      , par(par)
      , m(owner_->m)
      , eb(eb0)
      , edim_range_(m.GetEdim())
      , mebc_(mebc)
      , mcc_(mcc) {
    using namespace fluid_condition;

    UpdateDerivedConditions();

    fc_accel_.Reinit(m, Vect(0));
    fcvel_.time_curr = fcvel;

    fcp_.time_curr.Reinit(m, 0.);

    // Calc initial volume fluxes
    const auto ffwe = UEB::Interpolate(fcvel_.time_curr, me_vel_, eb);
    fev_.time_curr.Reinit(m, 0.);
    eb.LoopFaces([&](auto cf) { //
      fev_.time_curr[cf] = ffwe[cf].dot(eb.GetSurface(cf));
    });
    // Apply meshvel
    const Vect& meshvel = par.meshvel;
    eb.LoopFaces([&](auto cf) { //
      fev_.time_curr[cf] -= meshvel.dot(eb.GetSurface(cf));
    });
  }

  void UpdateDerivedConditions() {
    using namespace fluid_condition;

    me_vel_ = GetVelCond<M>(mebc_);
    mebc_.LoopPairs([&](auto p) {
      const auto cf = p.first;
      const auto& bc = p.second;
      const auto nci = bc.nci;
      me_pressure_[cf] = BCond<Scal>(BCondType::neumann, nci);
      me_visc_[cf] = BCond<Scal>(BCondType::neumann, nci);
      if (bc.type == BCondFluidType::slipwall ||
          bc.type == BCondFluidType::symm) {
        me_force_[cf] = BCond<Vect>(BCondType::mixed, nci);
      } else {
        me_force_[cf] = BCond<Vect>(BCondType::neumann, nci);
      }
    });

    mccp_.clear();
    for (auto& it : mcc_) {
      const IdxCell c = it.first;
      const CondCellFluid* cb = it.second.get();

      if (auto cd = dynamic_cast<const GivenPressure<M>*>(cb)) {
        mccp_[c] = std::make_shared<CondCellValFixed<Scal>>(cd->GetPressure());
      } else {
        throw std::runtime_error(FILELINE + ": unknown cell condition");
      }
    }

    is_boundary_.Reinit(m, false);
    mebc_.LoopPairs([&](auto p) { //
      is_boundary_[p.first] = true;
    });

    ffvisc_ = UEB::Interpolate(*owner_->fcd_, me_visc_, eb);
    ffdens_ = UEB::Interpolate(*owner_->fcr_, me_visc_, eb);
  }
  void Advection(
      FieldCell<Vect>& fcvel, const FieldCell<Vect>& fcvel_time_prev,
      const FieldFaceb<Scal>& ffv, const FieldCell<Vect>* fc_accel,
      const Scal dt) {
    auto sem = m.GetSem("advection");
    for (size_t d = 0; d < dim; ++d) {
      if (sem("local" + std::to_string(d))) {
        const auto mebc = GetScalarCond(me_vel_, d, m);
        const auto fcu = GetComponent(fcvel, d);
        const FieldCell<Vect> fcg =
            UEB::AverageGradient(UEB::Gradient(fcu, mebc, eb), eb);
        const FieldFaceb<Scal> ffu =
            UEB::InterpolateUpwind(fcu, mebc, par.convsc, fcg, ffv, eb);
        FieldCell<Scal> fc_sum(eb, 0);
        for (auto c : eb.Cells()) {
          Scal sum = 0;
          eb.LoopNci(c, [&](auto q) {
            const auto cf = eb.GetFace(c, q);
            const auto flux = -ffu[cf] * ffv[cf];
            sum += flux * eb.GetOutwardFactor(c, q);
          });
          fc_sum[c] = sum;
        }
        fc_sum = UEB::RedistributeCutCells(fc_sum, eb);
        for (auto c : eb.Cells()) {
          fcvel[c][d] =
              fcvel_time_prev[c][d] + fc_sum[c] * dt / eb.GetVolume(c);
        }
      }
    }
    if (sem("comm")) {
      if (fc_accel) {
        for (auto c : eb.Cells()) {
          fcvel[c] += (*fc_accel)[c] * dt;
        }
      }
      m.Comm(&fcvel);
    }
  }
  void DiffusionExplicit(
      FieldCell<Vect>& fcvel, const FieldCell<Vect>& fcvel_time_prev,
      const FieldCell<Scal>& fc_dens, const Scal dt) {
    auto sem = m.GetSem("diffusion");
    for (size_t d = 0; d < dim; ++d) {
      if (sem("local" + std::to_string(d))) {
        const auto mebc = GetScalarCond(me_vel_, d, m);
        const auto fcu = GetComponent(fcvel, d);
        const FieldFaceb<Scal> ffg = UEB::Gradient(fcu, mebc, eb);
        FieldCell<Scal> fc_sum(eb, 0);
        for (auto c : eb.Cells()) {
          Scal sum = 0;
          eb.LoopNci(c, [&](auto q) {
            const auto cf = eb.GetFace(c, q);
            const auto flux = ffg[cf] * ffvisc_[cf] * eb.GetArea(cf);
            sum += flux * eb.GetOutwardFactor(c, q);
          });
          fc_sum[c] = sum;
        }
        fc_sum = UEB::RedistributeCutCells(fc_sum, eb);
        for (auto c : eb.Cells()) {
          fcvel[c][d] = fcvel_time_prev[c][d] +
                        fc_sum[c] * dt / (fc_dens[c] * eb.GetVolume(c));
        }
      }
    }
    if (sem("comm")) {
      m.Comm(&fcvel);
    }
  }
  void DiffusionImplicit(
      FieldCell<Vect>& fcvel, const FieldCell<Vect>& fcvel_time_prev,
      const FieldCell<Scal>& fc_dens, const Scal dt) {
    auto sem = m.GetSem("diffusion");
    struct {
      FieldCell<Scal> fcu;
      FieldCell<Expr> fcl;
    } * ctx(sem);
    auto& fcu = ctx->fcu;
    auto& fcl = ctx->fcl;
    for (size_t d = 0; d < dim; ++d) {
      if (sem("local" + std::to_string(d))) {
        const auto mebc = GetScalarCond(me_vel_, d, m);
        fcu = GetComponent(fcvel, d);
        const FieldFaceb<ExprFace> ffg = UEB::GradientImplicit(fcu, mebc, eb);
        fcl.Reinit(eb, Expr::GetUnit(0));
        for (auto c : eb.Cells()) {
          Expr sum(0);
          eb.LoopNci(c, [&](auto q) {
            const auto cf = eb.GetFace(c, q);
            const ExprFace flux = ffg[cf] * ffvisc_[cf] * eb.GetArea(cf) *
                                  eb.GetOutwardFactor(c, q);
            eb.AppendExpr(sum, flux, q);
          });
          Expr td(0); // time derivative
          td[0] = 1 / dt;
          td[Expr::dim - 1] = -fcvel_time_prev[c][d] / dt;
          fcl[c] = td * eb.GetVolume(c) * fc_dens[c] - sum;
        }
        fcl.SetName("velocity" + std::to_string(d));
      }
      if (sem.Nested("solve" + std::to_string(d))) {
        Solve(fcl, &fcu, fcu, M::LS::T::symm, m);
      }
      if (sem("copy")) {
        SetComponent(fcvel, d, fcu);
      }
    }
    if (sem("comm")) {
      m.Comm(&fcvel);
    }
  }

  void StartStep() {
    auto sem = m.GetSem("fluid-start");
    if (sem("convdiff-init")) {
      owner_->ClearIter();
      CHECKNAN(fcp_.time_curr, m.CN())
    }

    if (sem("convdiff-start")) {
      // rotate layers
      fcp_.iter_curr = fcp_.time_curr;
      fev_.iter_curr = fev_.time_curr;
      fcvel_.iter_curr = fcvel_.time_curr;
    }
  }
  // Apply cell conditions for pressure.
  // fcpb: base pressure [i]
  // fcs: linear system for pressure [i]
  void ApplyCellCond(const FieldCell<Scal>& fcpb, FieldCell<Expr>& fcs) {
    (void)fcpb;
    for (auto& it : mccp_) {
      const IdxCell c = it.first; // target cell
      const CondCell* cb = it.second.get(); // cond base
      if (auto cd = dynamic_cast<const CondCellVal<Scal>*>(cb)) {
        auto& e = fcs[c];
        const Scal pc = cd->second(); // new value for p[c]
        e[0] += 1;
        e[Expr::dim - 1] -= pc;
      }
    }
  }
  // Adds pressure gradient term to flux.
  // fcp: pressure
  // fck: diagonal coefficient
  // ffv: flux
  FieldFaceb<ExprFace> GetFlux(
      const FieldCell<Scal>& fcp, const FieldFaceb<Scal>& ffv, Scal dt) {
    FieldFaceb<ExprFace> ff_flux = UEB::GradientImplicit(fcp, me_pressure_, eb);
    eb.LoopFaces([&](auto cf) { //
      ff_flux[cf] *= -eb.GetArea(cf) / ffdens_[cf] * dt;
      ff_flux[cf][2] += ffv[cf];
    });
    return ff_flux;
  }
  // Expressions for sum of fluxes and source:
  //   sum(v) - source * volume
  // ffv: fluxes [i]
  // fcsv: volume source [i]
  // Output:
  // fce: result [i]
  FieldCell<Expr> GetFluxSum(
      const FieldFaceb<ExprFace>& ffv, const FieldCell<Scal>& fcsv) {
    // initialize as diagonal system
    FieldCell<Expr> fce(m, Expr::GetUnit(0));
    for (auto c : eb.Cells()) {
      Expr sum(0);
      eb.LoopNci(c, [&](auto q) {
        const auto cf = eb.GetFace(c, q);
        const ExprFace v = ffv[cf] * eb.GetOutwardFactor(c, q);
        eb.AppendExpr(sum, v, q);
      });
      sum[Expr::dim - 1] -= fcsv[c] * eb.GetVolume(c);
      fce[c] = sum;
    }
    return fce;
  }
  // Append explicit part of viscous force.
  // fcvel: velocity [a]
  // Output:
  // fcf += viscous term [i]
  void AppendExplViscous(
      const FieldCell<Vect>& fcvel, FieldCell<Vect>& fcf, const M& m) {
    const auto wf = UEB::Interpolate(fcvel, me_vel_, m);
    for (auto d : edim_range_) {
      const auto wfo = GetComponent(wf, d);
      const auto gc = ::Gradient(wfo, m);
      const auto gf = UEB::Interpolate(gc, me_force_, m);
      for (auto c : eb.Cells()) {
        Vect s(0);
        for (auto q : eb.Nci(c)) {
          const IdxFace f = m.GetFace(c, q);
          s += gf[f] * (ffvisc_[f] * m.GetOutwardSurface(c, q)[d]);
        }
        fcf[c] += s / m.GetVolume(c);
      }
    }
  }
  void Project(FieldCell<Scal>& fcp, FieldFaceb<Scal>& ffv, Scal dt) {
    auto sem = m.GetSem("project");
    struct {
      FieldFaceb<ExprFace> ffvc; // expression for corrected volume flux [i]
      FieldCell<Expr> fcpcs; // pressure correction linear system [i]
    } * ctx(sem);
    if (sem("local")) {
      ctx->ffvc = GetFlux(fcp, ffv, dt);
      ctx->fcpcs = GetFluxSum(ctx->ffvc, *owner_->fcsv_);
      ApplyCellCond(fcp, ctx->fcpcs);
      ctx->fcpcs.SetName("pressure");
    }
    if (sem.Nested("solve")) {
      Solve(ctx->fcpcs, &fcp, fcp, M::LS::T::symm, m);
    }
    if (sem("fluxes")) {
      eb.LoopFaces([&](auto cf) { //
        ffv[cf] = UEB::Eval(ctx->ffvc[cf], cf, fcp, eb);
      });
    }
  }
  FieldCell<Vect> GetAcceleration(const FieldCell<Scal>& fcp) {
    const auto ffgp = UEB::Gradient(fcp, me_pressure_, eb);
    FieldFaceb<Scal> ff_forcep(eb, 0);
    eb.LoopFaces([&](auto cf) { //
      ff_forcep[cf] = (*owner_->febp_)[cf] - ffgp[cf];
    });
    auto fca = UEB::AverageGradient(ff_forcep, eb);
    // cell force
    for (auto c : eb.Cells()) {
      fca[c] += (*owner_->fcf_)[c];
    }
    // AppendExplViscous(fcvel_.iter_curr, fc_accel_, eb); XXX
    // divide by density
    for (auto c : eb.Cells()) {
      fca[c] /= (*owner_->fcr_)[c];
    }
    return fca;
  }
  void AppendExplViscous(
      const FieldCell<Vect>& fcvel, FieldCell<Vect>& fcf, const Embed<M>& eb) {
    // FIXME: not implemented
    (void)fcvel;
    (void)fcf;
    (void)eb;
    return;
  }
  void UpdateBc(typename M::Sem& sem) {
    if (sem("bc-derived")) {
      UpdateDerivedConditions();
    }
    if (sem.Nested("bc-inletflux")) {
      UFluid<M>::UpdateInletFlux(
          m, eb, fcvel_.iter_curr, mebc_, par.inletflux_numid, me_vel_);
    }
    if (sem.Nested("bc-outlet")) {
      UFluid<M>::template UpdateOutletVelocity(
          m, eb, fcvel_.iter_curr, mebc_, *owner_->fcsv_, me_vel_);
    }
  }
  void MakeIteration() {
    auto sem = m.GetSem("fluid-iter");
    auto& fcp_prev = fcp_.iter_prev;
    auto& fcp_curr = fcp_.iter_curr;
    const Scal dt = owner_->GetTimeStep();
    auto& ffv = fev_.iter_curr.template Get<FieldFaceb<Scal>>();
    auto& ffvt = fev_.time_curr.template Get<FieldFaceb<Scal>>();
    auto& fcvel = fcvel_.iter_curr;
    if (sem("init")) {
      // rotate layers
      fcp_prev = fcp_curr;
      fev_.iter_prev = fev_.iter_curr;
      fcvel_.iter_prev = fcvel_.iter_curr;
    }

    UpdateBc(sem);

    if (sem("forceinit")) {
      fc_accel_ = GetAcceleration(fcp_curr);
    }
    if (sem.Nested("proj")) {
      // ff_accel_: acceleration from force terms and pressure gradient
      // fcvelmid_: initial guess for velocity at t+dt/2 to estimate fluxes
      // ffv: divergence-free volume flux at t+dt/2 from previous iteration
      Advection(fcvel, fcvel_.time_curr, ffvt, &fc_accel_, dt * 0.5);
      // fcvelmid_: predicted velocity at t+dt/2
    }
    if (sem()) {
      const FieldFaceb<Vect> ffvel = UEB::Interpolate(fcvel, me_vel_, eb);
      eb.LoopFaces([&](auto cf) { //
        ffv[cf] = ffvel[cf].dot(eb.GetSurface(cf));
      });
      // ffv: predicted volume flux at t+dt/2
    }
    if (sem.Nested("project")) {
      Project(fcp_curr, ffv, dt * 0.5);
      // ffv: divergence-free predicted volume flux at t+dt/2
    }
    if (sem()) {
      fcvel = fcvel_.time_curr;
    }
    if (sem.Nested("proj")) {
      // ff_accel_: acceleration from force terms and pressure gradient
      // ffv: divergence-free volume flux at t+dt/2 from previous iteration
      Advection(fcvel, fcvel_.time_curr, ffv, &fc_accel_, dt);
      // fcvel: velocity at t+dt
    }
    if (sem.Nested("diffusion")) {
      switch (par.conv) {
        case Conv::exp:
          DiffusionExplicit(fcvel, fcvel, *owner_->fcr_, dt);
          break;
        case Conv::imp:
          DiffusionImplicit(fcvel, fcvel, *owner_->fcr_, dt);
          break;
      }
    }
    if (sem("subtract")) {
      // subtract old acceleration from cell velocity
      for (auto c : eb.Cells()) {
        fcvel[c] -= fc_accel_[c] * dt;
      }
      m.Comm(&fcvel);
    }
    if (sem("face-acceleration")) {
      // compute volume flux from center velocity
      // and add acceleration from body force
      const FieldFaceb<Vect> ffvel = UEB::Interpolate(fcvel, me_vel_, eb);
      eb.LoopFaces([&](auto cf) { //
        auto& v = ffv[cf];
        v = ffvel[cf].dot(eb.GetSurface(cf));
        if (!is_boundary_[cf]) {
          v += (*owner_->febp_)[cf] / ffdens_[cf] * dt * eb.GetArea(cf);
        }
      });
    }
    if (sem.Nested("project")) {
      Project(fcp_curr, ffv, dt);
      // ffv: divergence-free volume flux at t+dt
      // fcp_curr: final pressure
    }
    if (sem("pcorr-apply")) {
      fc_accel_ = GetAcceleration(fcp_curr);
      // Add new acceleration to cell velocity
      for (auto c : eb.Cells()) {
        fcvel[c] += fc_accel_[c] * dt;
      }
      m.Comm(&fcvel);
    }
    if (sem("iter_diff")) {
      iter_diff_ = 0;
      auto& fc = fcvel_.iter_curr;
      auto& fcm = fcvel_.iter_prev;
      for (auto c : eb.Cells()) {
        iter_diff_ = std::max(iter_diff_, (fc[c] - fcm[c]).norminf());
      }
      m.Reduce(&iter_diff_, "max");
      owner_->IncIter();
    }
  }
  void FinishStep() {
    auto sem = m.GetSem("fluid-finish");
    if (sem("inctime")) {
      fcp_.time_curr = fcp_.iter_curr;
      fev_.time_curr = fev_.iter_curr;
      fcvel_.time_curr = fcvel_.iter_curr;
      CHECKNAN(fcp_.time_curr, m.CN())
      owner_->IncTime();
    }
  }
  double GetAutoTimeStep() const {
    Scal dt = std::numeric_limits<Scal>::max();
    for (auto f : eb.Faces()) {
      const Scal vel = fev_.time_curr[f] / m.GetArea(f);
      if (vel != 0.) {
        dt = std::min<Scal>(dt, std::abs(m.GetCellSize()[0] / vel));
      }
    }
    return dt;
  }
  const FieldCell<Vect>& GetVelocity(Step l) const {
    return fcvel_.Get(l);
  }

  Owner* owner_;
  Par par;
  M& m;
  const EB& eb;
  const GRange<size_t> edim_range_; // effective dimension range

  // Face conditions
  const MapEmbed<BCondFluid<Vect>>& mebc_;
  MapEmbed<BCond<Vect>> me_vel_;
  MapEmbed<BCond<Scal>> me_pressure_;
  MapEmbed<BCond<Vect>> me_force_;
  MapEmbed<BCond<Scal>> me_visc_;

  // Cell conditions
  MapCell<std::shared_ptr<CondCellFluid>> mcc_; // fluid cell cond
  MapCell<std::shared_ptr<CondCell>> mccp_; // pressure cell cond

  StepData<FieldEmbed<Scal>> fev_; // volume flux
  StepData<FieldCell<Scal>> fcp_; // pressure
  StepData<FieldCell<Vect>> fcvel_; // velocity

  FieldEmbed<bool> is_boundary_; // true on faces with boundary conditions

  // Cell fields:
  FieldCell<Vect> fc_accel_; // acceleration

  // Face fields:
  FieldFaceb<Scal> ffvisc_; // dynamic viscosity
  FieldFaceb<Scal> ffdens_; // density

  Scal iter_diff_;
};

template <class EB_>
Proj<EB_>::Proj(
    M& m, const EB& eb, const FieldCell<Vect>& fcvel,
    const MapEmbed<BCondFluid<Vect>>& mebc,
    const MapCell<std::shared_ptr<CondCellFluid>>& mcc,
    const FieldCell<Scal>* fcr, const FieldCell<Scal>* fcd,
    const FieldCell<Vect>* fcf, const FieldEmbed<Scal>* febp,
    const FieldCell<Scal>* fcsv, const FieldCell<Scal>* fcsm, double t,
    double dt, Par par)
    : Base(t, dt, m, fcr, fcd, fcf, febp, fcsv, fcsm)
    , imp(new Imp(this, eb, fcvel, mebc, mcc, par)) {}

template <class EB_>
Proj<EB_>::~Proj() = default;

template <class EB_>
auto Proj<EB_>::GetPar() const -> const Par& {
  return imp->par;
}

template <class EB_>
void Proj<EB_>::SetPar(Par par) {
  imp->par = par;
}

template <class EB_>
void Proj<EB_>::StartStep() {
  return imp->StartStep();
}

template <class EB_>
void Proj<EB_>::MakeIteration() {
  return imp->MakeIteration();
}

template <class EB_>
void Proj<EB_>::FinishStep() {
  return imp->FinishStep();
}

template <class EB_>
auto Proj<EB_>::GetVelocity(Step l) const -> const FieldCell<Vect>& {
  return imp->GetVelocity(l);
}

template <class EB_>
auto Proj<EB_>::GetPressure(Step l) const -> const FieldCell<Scal>& {
  return imp->fcp_.Get(l);
}

template <class EB_>
auto Proj<EB_>::GetVolumeFlux(Step l) const -> const FieldEmbed<Scal>& {
  return imp->fev_.Get(l);
}

template <class EB_>
double Proj<EB_>::GetAutoTimeStep() const {
  return imp->GetAutoTimeStep();
}

template <class EB_>
double Proj<EB_>::GetError() const {
  return imp->iter_diff_;
}

template <class EB_>
auto Proj<EB_>::GetVelocityCond() const -> const MapEmbed<BCond<Vect>>& {
  return imp->me_vel_;
}
