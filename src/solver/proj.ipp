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
#include "convdiffv.h"
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
  using CD = ConvDiffVect<EB>; // convdiff solver
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

    fcfcd_.Reinit(m, Vect(0));
    typename CD::Par p;
    SetConvDiffPar(p, par);
    cd_ = GetConvDiff<EB>()(
        par.conv, m, eb, fcvel, me_vel_, owner_->fcr_, &ffd_, &fcfcd_,
        &fev_.iter_prev, owner_->GetTime(), owner_->GetTimeStep(), p);

    fcp_.time_curr.Reinit(m, 0.);
    fcp_.time_prev = fcp_.time_curr;

    // Calc initial volume fluxes
    const auto ffwe = UEB::Interpolate(cd_->GetVelocity(), me_vel_, eb);
    fev_.time_curr.Reinit(m, 0.);
    eb.LoopFaces([&](auto cf) { //
      fev_.time_curr[cf] = ffwe[cf].dot(eb.GetSurface(cf));
    });
    // Apply meshvel
    const Vect& meshvel = par.meshvel;
    eb.LoopFaces([&](auto cf) { //
      fev_.time_curr[cf] -= meshvel.dot(eb.GetSurface(cf));
    });

    fev_.time_prev = fev_.time_curr;
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
  }
  void StartStep() {
    auto sem = m.GetSem("fluid-start");
    if (sem("convdiff-init")) {
      owner_->ClearIter();
      CHECKNAN(fcp_.time_curr, m.CN())
      cd_->SetTimeStep(owner_->GetTimeStep());
    }

    if (sem.Nested("convdiff-start")) {
      cd_->StartStep();
    }

    if (sem("convdiff-start")) {
      // rotate layers
      fcp_.iter_curr = fcp_.time_curr;
      fev_.iter_curr = fev_.time_curr;
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
      const FieldCell<Scal>& fcp, const FieldEmbed<Scal>& fek,
      const FieldEmbed<Scal>& ffv) {
    FieldFaceb<ExprFace> ff_flux = UEB::GradientImplicit(fcp, me_pressure_, eb);
    eb.LoopFaces([&](auto cf) { //
      ff_flux[cf] *= -eb.GetArea(cf) / fek[cf];
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
  // Get diagcoeff from current convdiff equations
  void GetDiagCoeff(FieldEmbed<Scal>& fek) {
    auto sem = m.GetSem("diag");
    struct {
      FieldCell<Scal> fck;
    } * ctx(sem);
    auto& fck = ctx->fck;
    if (sem("local")) {
      fck.Reinit(m, 0);
      for (auto d : edim_range_) {
        const auto& fct = cd_->GetDiag(d);
        for (auto c : eb.Cells()) {
          fck[c] += fct[c];
        }
      }
      for (auto c : eb.Cells()) {
        fck[c] /= edim_range_.size();
      }
      CHECKNAN(fck, m.CN());
      m.Comm(&fck);
    }
    if (sem("interp")) {
      fek = UEB::Interpolate(fck, me_visc_, eb);
    }
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
          s += gf[f] * (ffd_[f] * m.GetOutwardSurface(c, q)[d]);
        }
        fcf[c] += s / m.GetVolume(c);
      }
    }
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
          m, eb, cd_->GetVelocity(Step::iter_curr), mebc_, par.inletflux_numid,
          me_vel_);
    }
    if (sem.Nested("bc-outlet")) {
      UFluid<M>::template UpdateOutletVelocity(
          m, eb, cd_->GetVelocity(Step::iter_curr), mebc_, *owner_->fcsv_,
          me_vel_);
    }
  }
  void MakeIteration() {
    auto sem = m.GetSem("fluid-iter");
    struct {
      FieldCell<Vect> fcvel_corr;
      FieldFaceb<ExprFace> ffvc; // expression for corrected volume flux [i]
      FieldCell<Expr> fcpcs; // pressure correction linear system [i]
      FieldEmbed<Scal> fek; // diag coeff of velocity equation
    } * ctx(sem);
    auto& fcp_prev = fcp_.iter_prev;
    auto& fcp_curr = fcp_.iter_curr;
    if (sem("init")) {
      cd_->SetPar(UpdateConvDiffPar(cd_->GetPar(), par));

      // interpolate visosity
      ffd_ = UEB::Interpolate(*owner_->fcd_, me_visc_, eb);

      // rotate layers
      fcp_prev = fcp_curr;
      fev_.iter_prev = fev_.iter_curr;
    }

    UpdateBc(sem);

    if (sem("forceinit")) {
      fcfcd_.Reinit(m, Vect(0));
      for (auto c : eb.Cells()) {
        fcfcd_[c] += (*owner_->fcf_)[c];
      }
      AppendExplViscous(cd_->GetVelocity(Step::iter_curr), fcfcd_, eb);
    }

    if (sem.Nested("convdiff-iter")) {
      cd_->MakeIteration();
    }

    if (sem.Nested("diag")) {
      GetDiagCoeff(ctx->fek);
    }

    if (sem("pcorr-assemble")) {
      const auto& febp = *owner_->febp_;
      // Acceleration
      const FieldFaceb<Vect> ffvel =
          UEB::Interpolate(cd_->GetVelocity(Step::iter_curr), me_vel_, eb);
      eb.LoopFaces([&](auto cf) { //
        auto& v = fev_.iter_curr[cf];
        v = ffvel[cf].dot(eb.GetSurface(cf));
        if (!is_boundary_[cf]) { // inner
          v += febp[cf] * eb.GetArea(cf) / ctx->fek[cf];
        } else { // boundary
          // nop, keep the mean flux
        }
      });

      // Projection
      ctx->ffvc = GetFlux(fcp_curr, ctx->fek, fev_.iter_curr);
      ctx->fcpcs = GetFluxSum(ctx->ffvc, *owner_->fcsv_);
      ApplyCellCond(fcp_curr, ctx->fcpcs);
      ctx->fcpcs.SetName("pressure");
    }
    if (sem.Nested("pcorr-solve")) {
      Solve(ctx->fcpcs, &fcp_curr, fcp_curr, M::LS::T::symm, m);
    }
    if (sem("pcorr-apply")) {
      eb.LoopFaces([&](auto cf) {
        fev_.iter_curr[cf] = UEB::Eval(ctx->ffvc[cf], cf, fcp_curr, eb);
      });

      // Acceleration and correction of velocity
      const auto& febp = *owner_->febp_;
      auto ffgp = UEB::Gradient(fcp_curr, me_pressure_, eb);
      ctx->fcvel_corr.Reinit(m, Vect(0));
      for (auto c : eb.Cells()) {
        Vect sum(0);
        eb.LoopNci(c, [&](auto q) {
          const auto cf = eb.GetFace(c, q);
          if (!is_boundary_[cf]) {
            const Scal a = (febp[cf] - ffgp[cf]) / std::max(1e-8, ctx->fek[cf]);
            sum += eb.GetNormal(cf) * (a * 0.5);
          }
        });
        ctx->fcvel_corr[c] = sum;
      }
    }

    if (sem.Nested("convdiff-corr")) {
      cd_->CorrectVelocity(Step::iter_curr, ctx->fcvel_corr);
    }

    if (sem("inc-iter")) {
      owner_->IncIter();
    }
  }
  void FinishStep() {
    auto sem = m.GetSem("fluid-finish");
    if (sem("inctime")) {
      fcp_.time_prev = fcp_.time_curr;
      fev_.time_prev = fev_.time_curr;
      fcp_.time_curr = fcp_.iter_curr;
      fev_.time_curr = fev_.iter_curr;
      CHECKNAN(fcp_.time_curr, m.CN())
      owner_->IncTime();
    }
    if (sem.Nested("convdiff-finish")) {
      cd_->FinishStep();
    }
  }
  double GetAutoTimeStep() {
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
    return cd_->GetVelocity(l);
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

  std::shared_ptr<CD> cd_;

  FieldEmbed<bool> is_boundary_; // true on faces with boundary conditions

  // Cell fields:
  FieldCell<Vect> fcfcd_; // force for convdiff [i]

  // Face fields:
  FieldFaceb<Scal> ffd_; // dynamic viscosity
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
  return imp->cd_->GetError();
}

template <class EB_>
auto Proj<EB_>::GetVelocityCond() const -> const MapEmbed<BCond<Vect>>& {
  return imp->me_vel_;
}
