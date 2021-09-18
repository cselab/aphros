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

  Imp(Owner* owner, const EB& eb_, const Args& args)
      : owner_(owner)
      , par(args.par)
      , m(owner_->m)
      , eb(eb_)
      , edim_range_(m.GetEdim())
      , linsolver_(args.linsolver)
      , mebc_(args.mebc)
      , mcc_(args.mcc) {
    UpdateDerivedConditions();

    fc_accel_.Reinit(m, Vect(0));
    fcvel_.time_curr = args.fcvel;

    fcvel_diffusion_guess_ = fcvel_.time_curr;

    fcp_.time_curr.Reinit(m, 0.);
    fcp_predict_ = fcp_.time_curr;

    // Calc initial volume fluxes
    const auto ffwe = UEB::Interpolate(fcvel_.time_curr, me_vel_, eb);
    fev_.time_curr.Reinit(m, 0.);
    eb.LoopFaces([&](auto cf) { //
      fev_.time_curr[cf] = (ffwe[cf] - par.meshvel).dot(eb.GetSurface(cf));
    });
  }

  void UpdateDerivedConditions() {
    me_vel_ = ConvertBCondFluidToVelocity<M>(mebc_);
    mebc_.LoopPairs([&](auto p) {
      const auto cf = p.first;
      const auto& bc = p.second;
      const auto nci = bc.nci;
      me_visc_[cf] = BCond<Scal>(BCondType::neumann, nci);
      me_pressure_flux_[cf] = BCond<Scal>(BCondType::neumann, nci);
      if (bc.type == BCondFluidType::slipwall ||
          bc.type == BCondFluidType::symm) {
        me_force_[cf] = BCond<Vect>(BCondType::mixed, nci);
        me_pressure_[cf] = BCond<Scal>(BCondType::neumann, nci);
      } else {
        me_force_[cf] = BCond<Vect>(BCondType::neumann, nci);
        me_pressure_[cf] = BCond<Scal>(BCondType::extrap, nci);
      }
    });

    mccp_.clear();
    for (auto& it : mcc_) {
      const IdxCell c = it.first;
      const CondCellFluid* cb = it.second.get();

      using namespace fluid_condition;
      if (auto cd = dynamic_cast<const GivenPressure<M>*>(cb)) {
        mccp_[c] = std::make_shared<CondCellValFixed<Scal>>(cd->GetPressure());
      } else {
        fassert(false, "unknown cell condition");
      }
    }

    is_boundary_.Reinit(m, false);
    mebc_.LoopPairs([&](auto p) { //
      is_boundary_[p.first] = true;
    });

    ffvisc_ = UEB::Interpolate(*owner_->fcd_, me_visc_, eb);
    ffdens_ = UEB::InterpolateHarmonic(*owner_->fcr_, me_visc_, eb);
  }
  void CalcInitialPressure(
      FieldCell<Scal>& fcp, FieldFaceb<Scal>& ffv,
      const FieldCell<Vect>& fcvel) {
    auto sem = m.GetSem("initial-pressure");
    const Scal dt = owner_->GetTimeStep();
    if (sem("face-acceleration")) {
      const FieldFaceb<Vect> ffvel = UEB::Interpolate(fcvel, me_vel_, eb);
      eb.LoopFaces([&](auto cf) { //
        auto& v = ffv[cf];
        v = (ffvel[cf] - par.meshvel).dot(eb.GetSurface(cf));
        if (!is_boundary_[cf]) {
          v += (*owner_->febp_)[cf] / ffdens_[cf] * dt * eb.GetArea(cf);
        }
      });
    }
    if (sem.Nested("project")) {
      Project(fcp, ffv, dt);
    }
  }
  void Advection(
      FieldCell<Vect>& fcvel, const FieldCell<Vect>& fcvel_time_prev,
      const FieldFaceb<Scal>& ffv, const FieldCell<Vect>& fc_accel,
      const Scal dt) {
    auto sem = m.GetSem("advection");
    struct {
      FieldCell<Scal> fc_sum;
    } * ctx(sem);
    auto& fc_sum = ctx->fc_sum;
    for (size_t d = 0; d < dim; ++d) {
      if (sem("local" + std::to_string(d))) {
        const auto mebc = GetScalarCond(me_vel_, d, m);
        fcvel.CheckHalo(1);
        const auto fcu = GetComponent(fcvel, d);
        FieldFaceb<Scal> ffu;
        if (par.bcg) {
          const auto fc_src = GetComponent(fc_accel, d);
          ffu = UEB::InterpolateBcg(fcu, mebc, ffv, fc_src, dt, eb);
        } else {
          const FieldCell<Vect> fcg =
              UEB::AverageGradient(UEB::Gradient(fcu, mebc, eb), eb);
          ffu = UEB::InterpolateUpwind(fcu, mebc, par.convsc, fcg, ffv, eb);
        }
        ffu.SetName(FILELINE + ":ffu");
        ffu.CheckHalo(0);
        fc_sum.Reinit(eb, 0);
        for (auto c : eb.Cells()) {
          Scal sum = 0;
          eb.LoopNci(c, [&](auto q) {
            const auto cf = eb.GetFace(c, q);
            const auto flux = -ffu[cf] * ffv[cf];
            sum += flux * eb.GetOutwardFactor(c, q);
          });
          fc_sum[c] = sum;
        }
        m.Comm(&fc_sum);
      }
      if (sem()) {
        // FIXME: empty stage to finish communication in halo cells
        // fc_sum gets different buffer in next assignment
      }
      if (sem("redist" + std::to_string(d))) {
        if (par.redistr_adv) {
          fc_sum = UEB::RedistributeCutCellsAdvection(fc_sum, ffv, 1, dt, eb);
        } else {
          fc_sum = UEB::RedistributeCutCells(fc_sum, eb);
        }
        for (auto c : eb.Cells()) {
          fcvel[c][d] =
              fcvel_time_prev[c][d] + fc_sum[c] * dt / eb.GetVolume(c);
        }
      }
    }
    if (sem("comm")) {
      for (auto c : eb.Cells()) {
        fcvel[c] += fc_accel[c] * dt;
      }
      m.Comm(&fcvel);
      fcvel.SetHalo(2);
    }
    if (sem()) {
      // FIXME: empty stage to finish communication in halo cells
    }
  }
  void DiffusionExplicit(
      FieldCell<Vect>& fcvel, const FieldCell<Vect>& fcvel_time_prev,
      const FieldCell<Scal>& fc_dens, const Scal dt) {
    auto sem = m.GetSem("diffusion");
    struct {
      FieldCell<Scal> fc_sum;
    } * ctx(sem);
    auto& fc_sum = ctx->fc_sum;
    for (size_t d = 0; d < dim; ++d) {
      if (sem("local" + std::to_string(d))) {
        const auto mebc = GetScalarCond(me_vel_, d, m);
        const auto fcu = GetComponent(fcvel, d);
        const FieldFaceb<Scal> ffg = UEB::Gradient(fcu, mebc, eb);
        ffg.CheckHalo(0);
        fc_sum.Reinit(eb, 0);
        for (auto c : eb.Cells()) {
          Scal sum = 0;
          eb.LoopNci(c, [&](auto q) {
            const auto cf = eb.GetFace(c, q);
            const auto flux = ffg[cf] * ffvisc_[cf] * eb.GetArea(cf);
            sum += flux * eb.GetOutwardFactor(c, q);
          });
          fc_sum[c] = sum;
        }
        m.Comm(&fc_sum);
      }
      if (sem()) {
        // FIXME: empty stage to finish communication in halo cells
        // fc_sum gets different buffer in the next stage
      }
      if (sem("redist" + std::to_string(d))) {
        fc_sum = UEB::RedistributeCutCells(fc_sum, eb);
        for (auto c : eb.Cells()) {
          fcvel[c][d] = fcvel_time_prev[c][d] +
                        fc_sum[c] * dt / (fc_dens[c] * eb.GetVolume(c));
        }
      }
    }
    if (sem("comm")) {
      m.Comm(&fcvel);
      fcvel.SetHalo(2);
    }
    if (sem()) {
      // FIXME: empty stage to finish communication in halo cells
    }
  }
  void DiffusionImplicit(
      FieldCell<Vect>& fcvel, const FieldCell<Vect>& fcvel_time_prev,
      const FieldCell<Vect>& fcvel_guess, const FieldCell<Scal>& fc_dens,
      const Scal dt) {
    auto sem = m.GetSem("diffusion");
    struct {
      FieldCell<Scal> fcu;
      FieldCell<Scal> fcu_guess;
      FieldCell<Expr> fcl;
    } * ctx(sem);
    auto& fcu = ctx->fcu;
    auto& fcu_guess = ctx->fcu_guess;
    auto& fcl = ctx->fcl;
    for (size_t d = 0; d < dim; ++d) {
      if (sem("local" + std::to_string(d))) {
        fcu = GetComponent(fcvel, d);
        fcu_guess = GetComponent(fcvel_guess, d);
      }
      for (size_t iter = 0; iter < par.diffusion_iters; ++iter) {
        if (sem("local" + std::to_string(d))) {
          const auto mebc = GetScalarCond(me_vel_, d, m);
          const FieldFaceb<ExprFace> ffg =
              UEB::GradientImplicit(fcu_guess, mebc, eb);
          fcl.Reinit(eb, Expr::GetUnit(0));
          ffg.CheckHalo(0);
          for (auto c : eb.Cells()) {
            Expr sum(0);
            eb.LoopNci(c, [&](auto q) {
              const auto cf = eb.GetFace(c, q);
              const ExprFace flux = ffg[cf] * ffvisc_[cf] * eb.GetArea(cf);
              eb.AppendExpr(sum, flux * eb.GetOutwardFactor(c, q), q);
            });
            fcl[c] = -sum;
          }
          fcl.SetName("velocity" + std::to_string(d));
        }
        if (sem("time")) {
          for (auto c : eb.Cells()) {
            Expr td(0); // time derivative
            td[0] = 1 / dt;
            td.back() = -fcvel_time_prev[c][d] / dt;
            fcl[c] += td * eb.GetVolume(c) * fc_dens[c];
          }
        }
        if (sem.Nested("solve")) {
          linsolver_->Solve(fcl, &fcu_guess, fcu, m);
        }
        if (par.diffusion_iters > 1) {
          if (sem("updguess")) {
            fcu_guess = fcu;
          }
        }
      }
      if (sem("copy")) {
        SetComponent(fcvel, d, fcu);
      }
    }
    if (sem("comm")) {
      m.Comm(&fcvel);
      fcvel.SetHalo(2);
    }
    if (sem()) {
      // FIXME: empty stage to finish communication in halo cells
    }
  }

  void StartStep() {
    auto sem = m.GetSem("fluid-start");
    if (sem("convdiff-init")) {
      owner_->ClearIter();
      CHECKNAN(fcp_.time_curr, m.flags.check_nan)
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
        e.back() -= pc;
      }
    }
  }
  // Adds pressure gradient term to flux.
  // fcp: pressure
  // fck: diagonal coefficient
  // ffv: flux
  FieldFaceb<ExprFace> GetFlux(const FieldFaceb<Scal>& ffv, Scal dt) {
    // implicit pressure gradient
    FieldFaceb<ExprFace> ffe = UEB::GradientImplicit({}, eb);
    mebc_.LoopBCond(eb, [&](auto cf, IdxCell, auto& bc) {
      if ((bc.type == BCondFluidType::inletpressure ||
           bc.type == BCondFluidType::outletpressure)) {
        ffe[cf].back() += ffe[cf][1 - bc.nci] * bc.pressure;
        ffe[cf][1 - bc.nci] = 0;
      } else {
        ffe[cf] = ExprFace(0);
      }
    });
    eb.LoopFaces([&](auto cf) { //
      ffe[cf] *= -eb.GetArea(cf) / ffdens_[cf] * dt;
      ffe[cf][2] += ffv[cf];
    });
    return ffe;
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
      sum.back() -= fcsv[c] * eb.GetVolume(c);
      fce[c] = sum;
    }
    return fce;
  }
  void Project(FieldCell<Scal>& fcp, FieldFaceb<Scal>& ffv, Scal dt) {
    auto sem = m.GetSem("project");
    struct {
      FieldFaceb<ExprFace> ffvc; // expression for corrected volume flux [i]
      FieldCell<Expr> fcpcs; // linear system for pressure [i]
    } * ctx(sem);
    if (sem("local")) {
      ctx->ffvc = GetFlux(ffv, dt);
      ctx->fcpcs = GetFluxSum(ctx->ffvc, *owner_->fcsv_);
      ApplyCellCond(fcp, ctx->fcpcs);
      ctx->fcpcs.SetName("pressure");
    }
    if (sem.Nested("solve")) {
      linsolver_->Solve(ctx->fcpcs, &fcp, fcp, m);
    }
    if (sem("fluxes")) {
      eb.LoopFaces([&](auto cf) { //
        ffv[cf] = UEB::Eval(ctx->ffvc[cf], cf, fcp, eb);
      });
    }
    if (sem()) {
      // FIXME: empty stage to finish communication in halo cells
    }
  }
  FieldCell<Vect> GetAcceleration(const FieldCell<Scal>& fcp) {
    const auto ffgp = UEB::Gradient(fcp, me_pressure_, eb);
    FieldFaceb<Scal> ffa(eb, 0);
    eb.LoopFaces([&](auto cf) { //
      ffa[cf] = ((*owner_->febp_)[cf] - ffgp[cf]) / ffdens_[cf];
    });
    ffa.SetHalo(0);
    auto fca = UEB::AverageGradient(ffa, eb);
    fca.SetName(FILELINE + ":fca");
    // cell force
    for (auto c : eb.Cells()) {
      fca[c] += (*owner_->fcf_)[c] / (*owner_->fcr_)[c];
    }
    if (par.explviscous) {
      FieldCell<Vect> fc_force(m, Vect(0));
      UFluid<M>::AppendExplViscous(
          fc_force, fcvel_.iter_curr, me_vel_, ffvisc_, eb);
      for (auto c : eb.Cells()) {
        fca[c] += fc_force[c] / (*owner_->fcr_)[c];
      }
    }
    fca.SetHalo(0);
    return fca;
  }
  void UpdateBc(typename M::Sem& sem) {
    if (sem("bc-derived")) {
      UpdateDerivedConditions();
    }
    if (sem.Nested("bc-inletflux")) {
      UFluid<M>::UpdateVelocityOnPressureBoundaries(
          me_vel_, m, eb, fev_.iter_curr, mebc_, par.outlet_relax);
    }
    if (sem.Nested("bc-inletflux")) {
      UFluid<M>::UpdateInletFlux(
          m, eb, fcvel_.iter_curr, mebc_, par.inletflux_numid, me_vel_);
    }
    if (sem.Nested("bc-outlet")) {
      UFluid<M>::template UpdateOutletVelocity(
          m, eb, fcvel_.iter_curr, mebc_, *owner_->fcsv_, par.outlet_relax,
          me_vel_);
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

    if (sem.Nested()) {
      if (owner_->GetTime() == 0) {
        CalcInitialPressure(fcp_curr, ffv, fcvel);
      }
    }

    if (sem("forceinit")) {
      fc_accel_ = GetAcceleration(fcp_curr);
      m.Comm(&fc_accel_);
      fc_accel_.SetHalo(2);
    }
    if (!par.stokes) {
      if (sem.Nested()) {
        CommFieldFace(ffvt, m);
      }
      if (sem("predicted-flux")) {
        FieldFaceb<Vect> ffvel(eb, Vect(0));
        for (size_t d = 0; d < dim; ++d) {
          const auto mebc = GetScalarCond(me_vel_, d, m);
          const auto fcu = GetComponent(fcvel_.time_curr, d);
          const auto fc_src = GetComponent(fc_accel_, d);
          const auto ffu = UEB::InterpolateBcg(fcu, mebc, ffvt, fc_src, dt, eb);
          // ffu: velocity advected by dt/2
          eb.LoopFaces([&](auto cf) { //
            ffvel[cf][d] = ffu[cf];
          });
        }
        eb.LoopFaces([&](auto cf) { //
          ffv[cf] = (ffvel[cf] - par.meshvel).dot(eb.GetSurface(cf));
        });
      }
      if (sem.Nested("project")) {
        Project(fcp_predict_, ffv, dt * 0.5);
      }
      if (sem()) {
        fcvel = fcvel_.time_curr;
      }
      if (sem.Nested()) {
        CommFieldFace(ffv, m);
      }
      if (sem.Nested("proj")) {
        Advection(fcvel, fcvel_.time_curr, ffv, fc_accel_, dt);
      }
    } else {
      if (sem("acceleration")) {
        for (auto c : eb.SuCells()) {
          fcvel[c] += fc_accel_[c] * dt;
        }
        m.Comm(&fcvel);
        fcvel.SetHalo(2);
      }
    }
    if (!par.diffusion_consistent_guess) {
      if (sem("diffusion-update-guess")) {
        fcvel_diffusion_guess_ = fcvel;
      }
    }
    if (sem.Nested("diffusion")) {
      switch (par.conv) {
        case Conv::exp:
          DiffusionExplicit(fcvel, fcvel, *owner_->fcr_, dt);
          break;
        case Conv::imp:
          DiffusionImplicit(
              fcvel, fcvel, fcvel_diffusion_guess_, *owner_->fcr_, dt);
          break;
      }
    }
    if (par.diffusion_consistent_guess) {
      if (sem("diffusion-update-guess")) {
        fcvel_diffusion_guess_ = fcvel;
      }
    }
    if (sem("face-acceleration")) {
      if (par.diffusion_consistent_guess) {
        fcvel_diffusion_guess_ = fcvel;
      }
      // subtract old acceleration from cell velocity
      for (auto c : eb.AllCells()) {
        fcvel[c] -= fc_accel_[c] * dt;
      }
      // compute volume flux from cell velocity
      // and add acceleration from body force
      // TODO: include cell body force
      const FieldFaceb<Vect> ffvel = UEB::Interpolate(fcvel, me_vel_, eb);
      eb.LoopFaces([&](auto cf) { //
        auto& v = ffv[cf];
        v = (ffvel[cf] - par.meshvel).dot(eb.GetSurface(cf));
        if (!is_boundary_[cf]) {
          v += (*owner_->febp_)[cf] / ffdens_[cf] * dt * eb.GetArea(cf);
        }
      });
    }
    if (sem.Nested("project")) {
      Project(fcp_curr, ffv, dt);
    }
    if (sem("cell-acceleration")) {
      // add new acceleration to cell velocity
      fc_accel_ = GetAcceleration(fcp_curr);
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
      CHECKNAN(fcp_.time_curr, m.flags.check_nan)
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
  std::shared_ptr<linear::Solver<M>> linsolver_;

  // Face conditions
  MapEmbed<BCondFluid<Vect>>& mebc_;
  MapEmbed<BCond<Vect>> me_vel_;
  MapEmbed<BCond<Scal>> me_pressure_;
  MapEmbed<BCond<Scal>> me_pressure_flux_; // conditions pressure in GetFlux
  MapEmbed<BCond<Vect>> me_force_;
  MapEmbed<BCond<Scal>> me_visc_;

  // Cell conditions
  MapCell<std::shared_ptr<CondCellFluid>> mcc_; // fluid cell cond
  MapCell<std::shared_ptr<CondCell>> mccp_; // pressure cell cond

  StepData<FieldEmbed<Scal>> fev_; // volume flux
  StepData<FieldCell<Scal>> fcp_; // pressure
  FieldCell<Scal> fcp_predict_; // pressure at prediction step,
                                // used as initial guess
  StepData<FieldCell<Vect>> fcvel_; // velocity

  FieldCell<Vect> fcvel_diffusion_guess_; // velocity after diffusion step,
                                          // used as initial guess

  FieldEmbed<bool> is_boundary_; // true on faces with boundary conditions

  // Cell fields:
  FieldCell<Vect> fc_accel_; // acceleration

  // Face fields:
  FieldFaceb<Scal> ffvisc_; // dynamic viscosity
  FieldFaceb<Scal> ffdens_; // density

  Scal iter_diff_;
};

template <class EB_>
Proj<EB_>::Proj(M& m_, const EB& eb, const Args& args)
    : Base(
          args.t, args.dt, m_, args.fcr, args.fcd, args.fcf, args.ffbp,
          args.fcsv, args.fcsm)
    , imp(new Imp(this, eb, args)) {}

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
