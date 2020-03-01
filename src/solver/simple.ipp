// Created by Petr Karnakov on 31.07.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <cmath>
#include <memory>
#include <sstream>
#include <stdexcept>

#include "approx.h"
#include "approx_eb.h"
#include "convdiffv.h"
#include "debug/isnan.h"
#include "fluid.h"
#include "linear/linear.h"
#include "simple.h"
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
struct Simple<EB_>::Imp {
  using Owner = Simple<EB_>;
  using CD = ConvDiffVect<EB>;
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
      , edim_range_(0, m.GetEdim())
      , drr_(m.GetEdim(), dim)
      , mebc_(mebc)
      , mcc_(mcc) {
    using namespace fluid_condition;

    UpdateDerivedConditions();

    fcfcd_.Reinit(m, Vect(0));
    typename CD::Par p;
    SetConvDiffPar(p, par);
    cd_ = GetConvDiff<EB>()(
        par.conv, m, eb, fcvel, me_vel_, owner_->fcr_, &ffvisc_, &fcfcd_,
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

    const auto ffp = UEB::Interpolate(fcp_.time_curr, me_pressure_, eb);
  }

  void UpdateDerivedConditions() {
    using namespace fluid_condition;

    me_vel_ = GetVelCond<M>(mebc_);
    mebc_.LoopPairs([&](auto p) {
      const auto cf = p.first;
      const auto& bc = p.second;
      const auto nci = bc.nci;
      me_visc_[cf] = BCond<Scal>(BCondType::neumann, nci);

      if (bc.type == BCondFluidType::slipwall ||
          bc.type == BCondFluidType::symm) {
        me_force_[cf] = BCond<Vect>(BCondType::mixed, nci);
        me_pressure_[cf] = BCond<Scal>(BCondType::neumann, nci);
        me_pcorr_[cf] = BCond<Scal>(BCondType::neumann, nci);
      } else {
        me_force_[cf] = BCond<Vect>(BCondType::neumann, nci);
        me_pressure_[cf] = BCond<Scal>(BCondType::extrap, nci);
        me_pcorr_[cf] = BCond<Scal>(BCondType::extrap, nci);
      }
    });

    mc_pressure_.clear();
    for (auto& it : mcc_) {
      const IdxCell c = it.first;
      const CondCellFluid* cb = it.second.get();

      if (auto cd = dynamic_cast<const GivenPressure<M>*>(cb)) {
        mc_pressure_[c] =
            std::make_shared<CondCellValFixed<Scal>>(cd->GetPressure());
      } else {
        throw std::runtime_error(FILELINE + ": unknown cell condition");
      }
    }

    is_boundary_.Reinit(m, false);
    mebc_.LoopPairs([&](auto p) { //
      is_boundary_[p.first] = true;
    });

    ffvisc_ = UEB::Interpolate(*owner_->fcd_, me_visc_, eb);
  }
  // Restore force from projections.
  // febp: force projections, bp=b.dot(n)
  // Output:
  // fcb: restored force
  void CalcExtForce(const FieldEmbed<Scal>& febp, FieldCell<Vect>& fcb) {
    auto sem = m.GetSem("extforce");
    if (sem("loc")) {
      fcb.Reinit(m);
      for (auto c : eb.Cells()) {
        Vect sum(0);
        eb.LoopNci(c, [&](auto q) {
          const auto cf = eb.GetFace(c, q);
          sum += eb.GetNormal(cf) * (febp[cf] * 0.5);
        });
        fcb[c] = sum;
      }
      m.Comm(&fcb);
    }
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
  // Rhie-Chow interpolation of predicted volume flux
  // including balanced force (hydrostatics and surface tension)
  // fcvel: predicted velocity field [s]
  // fcp: pressure field [s]
  // fcgp: gradient of pressure field [s]
  // fck, ffk: diag coeff [s]
  // Output:
  // fev: result [i]
  FieldFaceb<Scal> RhieChow(
      const FieldCell<Vect>& fcvel, const FieldCell<Scal>& fcp,
      const FieldCell<Vect>& fcgp, const FieldCell<Scal>& fck,
      const FieldFaceb<Scal>& fek) {
    auto fftv = UEB::Interpolate(fcvel, me_vel_, eb);
    FieldFaceb<Scal> fev(m);

    const Scal rh = par.rhie; // rhie factor
    auto& febp = *owner_->febp_;
    auto& fcbp = fcb_;

    eb.LoopFaces([&](auto cf) { //
      // mean flux
      auto& v = fev[cf];
      v = fftv[cf].dot(eb.GetSurface(cf));
      if (!is_boundary_[cf]) {
        const IdxCell cm = eb.GetCell(cf, 0);
        const IdxCell cp = eb.GetCell(cf, 1);

        // compact pressure gradient
        const Scal gp = (fcp[cp] - fcp[cm]) / eb.GetCellSize()[0];

        // compact
        const Scal o = (febp[cf] - gp) * eb.GetArea(cf) / fek[cf];

        // wide
        const Vect wm = (fcbp[cm] - fcgp[cm]) / fck[cm];
        const Vect wp = (fcbp[cp] - fcgp[cp]) / fck[cp];
        const Scal w = (wm + wp).dot(eb.GetSurface(cf)) * 0.5;

        // apply
        v += rh * (o - w);
      } else { // boundary
        // nop, keep mean flux
      }
    });

    // Apply meshvel
    eb.LoopFaces([&](auto cf) { //
      fev[cf] -= par.meshvel.dot(eb.GetSurface(cf));
    });
    return fev;
  }
  // Apply cell conditions for pressure.
  // fcpb: base pressure [i]
  // fcs: linear system in terms of correction of base pressure [i]
  void ApplyCellCond(const FieldCell<Scal>& fcpb, FieldCell<Expr>& fcs) {
    for (auto& it : mc_pressure_) {
      const IdxCell c = it.first; // target cell
      const CondCell* cb = it.second.get(); // cond base
      if (auto cd = dynamic_cast<const CondCellVal<Scal>*>(cb)) {
        auto& e = fcs[c];
        const Scal pc = cd->second() - fcpb[c]; // new value for p[c]
        e[0] += 1;
        e[Expr::dim - 1] -= pc;
      }
    }
  }

  // Flux expressions in terms of pressure:
  //   /  grad(p) * area / k + v, inner
  //   \  a, boundary
  // ffk: diag coeff [i]
  // fev: addition to flux [i]
  FieldFaceb<ExprFace> GetFlux(
      const FieldFaceb<Scal>& ffk, const FieldFaceb<Scal>& ffv) {
    FieldFaceb<ExprFace> ffe =
        UEB::GradientImplicit(FieldCell<Scal>(eb, 0), {}, eb);
    eb.LoopFaces([&](auto cf) { //
      if (!is_boundary_[cf]) {
        ffe[cf] *= -eb.GetArea(cf) / ffk[cf];
      } else {
        ffe[cf] = ExprFace(0);
      }
      ffe[cf][2] += ffv[cf];
    });
    return ffe;
  }
  // Expressions for sum of fluxes and source:
  //   sum(v) - sv * vol
  // fev: fluxes [i]
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
  void GetDiagCoeff(FieldCell<Scal>& fck, FieldFaceb<Scal>& ffk) {
    auto sem = m.GetSem("diag");
    if (sem("local")) {
      fck.Reinit(m, 0);
      for (auto d : edim_range_) {
        auto fct = cd_->GetDiag(d);
        for (auto c : eb.Cells()) {
          fck[c] += fct[c];
        }
      }
      for (auto c : eb.Cells()) {
        fck[c] /= edim_range_.size();
      }

      CHECKNAN(fck, m.CN())

      m.Comm(&fck);
    }
    if (sem("interp")) {
      ffk = UEB::Interpolate(fck, me_visc_, eb);
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
      const auto gc = Gradient(wfo, m);
      const auto gf = UEB::Interpolate(gc, me_force_, m);
      for (auto c : m.Cells()) {
        Vect s(0);
        for (auto q : m.Nci(c)) {
          IdxFace f = m.GetFace(c, q);
          s += gf[f] * (ffvisc_[f] * m.GetOutwardSurface(c, q)[d]);
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
      UFluid<M>::UpdateOutletVelocity(
          m, eb, cd_->GetVelocity(Step::iter_curr), mebc_, *owner_->fcsv_,
          me_vel_);
    }
  }
  void CalcForce(typename M::Sem& sem) {
    if (sem.Nested("forceinit")) {
      CalcExtForce(*owner_->febp_, fcb_);
    }

    if (sem("forceinit")) {
      // initialize force for convdiff
      fcfcd_.Reinit(m, Vect(0));
      // append force and balanced force
      auto& fcf = *owner_->fcf_;
      for (auto c : eb.Cells()) {
        fcfcd_[c] += fcf[c] + fcb_[c];
      }
    }

    if (sem("explvisc")) {
      AppendExplViscous(cd_->GetVelocity(Step::iter_curr), fcfcd_, eb);
    }
  }
  // TODO: rewrite norm() using dist() where needed
  void MakeIteration() {
    auto sem = m.GetSem("fluid-iter");
    struct {
      FieldCell<Vect> fcgp; // pressure gradient
      FieldCell<Expr> fcpcs; // pressure correction linear system
      FieldFaceb<ExprFace> ffvc; // expression for corrected volume flux [i]
      FieldCell<Vect> fcvel_corr;
      FieldCell<Scal> fck; // diag coeff of velocity equation
      FieldCell<Scal> fcpc; // pressurue correction
      FieldFaceb<Scal> ffk; // diag coeff of velocity equation
    } * ctx(sem);
    auto& fcgp = ctx->fcgp;
    auto& fcvel_corr = ctx->fcvel_corr;
    auto& fck = ctx->fck;
    auto& ffk = ctx->ffk;
    auto& fcpcs = ctx->fcpcs;
    auto& ffvc = ctx->ffvc;
    auto& fcpc = ctx->fcpc;
    auto& fcp_prev = fcp_.iter_prev;
    auto& fcp_curr = fcp_.iter_curr;
    if (sem("init")) {
      cd_->SetPar(UpdateConvDiffPar(cd_->GetPar(), par));

      // rotate layers
      fcp_prev = fcp_curr;
      fev_.iter_prev = fev_.iter_curr;
    }

    UpdateBc(sem);

    CalcForce(sem);

    if (sem("pgrad")) {
      fcgp =
          UEB::AverageGradient(UEB::Gradient(fcp_curr, me_pressure_, eb), eb);

      // append pressure gradient to force
      for (auto c : eb.Cells()) {
        fcfcd_[c] += fcgp[c] * (-1.);
      }
    }

    if (sem.Nested("convdiff-iter")) {
      cd_->MakeIteration();
    }

    if (sem.Nested("diag")) {
      GetDiagCoeff(fck, ffk);
    }

    if (sem("pcorr-assemble")) {
      const auto fev_rhie =
          RhieChow(cd_->GetVelocity(Step::iter_curr), fcp_curr, fcgp, fck, ffk);
      CHECKNAN(fev_rhie, m.CN())

      ffvc = GetFlux(ffk, fev_rhie);
      CHECKNAN(ffvc, m.CN())

      fcpcs = GetFluxSum(ffvc, *owner_->fcsv_);
      CHECKNAN(fcpcs, m.CN())

      ApplyCellCond(fcp_curr, fcpcs);
      fcpcs.SetName("pressure");
    }

    if (sem.Nested("pcorr-solve")) {
      Solve(fcpcs, nullptr, fcpc, M::LS::T::symm, m);
    }
    if (sem("pcorr-apply")) {
      CHECKNAN(fcpc, m.CN())

      // Correct pressure
      Scal pr = par.prelax; // pressure relaxation
      for (auto c : eb.AllCells()) {
        fcp_curr[c] += pr * fcpc[c];
      }

      // Calc divergence-free volume fluxes
      eb.LoopFaces([&](auto cf) {
        fev_.iter_curr[cf] = UEB::Eval(ffvc[cf], cf, fcpc, eb);
      });
      CHECKNAN(fev_.iter_curr, m.CN())

      const auto fcgpc =
          UEB::AverageGradient(UEB::Gradient(fcpc, me_pcorr_, eb), eb);

      // Calc velocity correction
      fcvel_corr.Reinit(m);
      for (auto c : eb.Cells()) {
        fcvel_corr[c] = -fcgpc[c] / fck[c];
      }
      CHECKNAN(fcvel_corr, m.CN())
    }
    if (sem.Nested("convdiff-corr")) {
      cd_->CorrectVelocity(Step::iter_curr, fcvel_corr);
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
    return cd_->GetVelocity(l);
  }

  Owner* owner_;
  Par par;
  M& m;
  const EB& eb;
  const GRange<size_t> edim_range_; // effective dimension range
  const GRange<size_t> drr_; // remaining dimensions

  // Boundary conditions
  const MapEmbed<BCondFluid<Vect>>& mebc_;
  MapEmbed<BCond<Vect>> me_vel_;
  MapEmbed<BCond<Scal>> me_pressure_;
  MapEmbed<BCond<Scal>> me_pcorr_;
  MapEmbed<BCond<Vect>> me_force_;
  MapEmbed<BCond<Scal>> me_visc_;

  // Cell conditions
  MapCell<std::shared_ptr<CondCellFluid>> mcc_; // fluid cell cond
  MapCell<std::shared_ptr<CondCell>> mc_pressure_; // pressure cell cond

  StepData<FieldEmbed<Scal>> fev_; // volume flux
  StepData<FieldCell<Scal>> fcp_; // pressure

  std::unique_ptr<CD> cd_;

  FieldEmbed<bool> is_boundary_; // true on faces with boundary conditions

  // Cell fields:
  FieldCell<Vect> fcb_; // restored balanced force [s]
  FieldCell<Vect> fcfcd_; // force for convdiff [i]

  // Face fields:
  FieldFaceb<Scal> ffvisc_; // dynamic viscosity
};

template <class EB_>
Simple<EB_>::Simple(
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
Simple<EB_>::~Simple() = default;

template <class EB_>
auto Simple<EB_>::GetPar() const -> const Par& {
  return imp->par;
}

template <class EB_>
void Simple<EB_>::SetPar(Par par) {
  imp->par = par;
}

template <class EB_>
void Simple<EB_>::StartStep() {
  return imp->StartStep();
}

template <class EB_>
void Simple<EB_>::MakeIteration() {
  return imp->MakeIteration();
}

template <class EB_>
void Simple<EB_>::FinishStep() {
  return imp->FinishStep();
}

template <class EB_>
auto Simple<EB_>::GetVelocity(Step l) const -> const FieldCell<Vect>& {
  return imp->GetVelocity(l);
}

template <class EB_>
auto Simple<EB_>::GetPressure(Step l) const -> const FieldCell<Scal>& {
  return imp->fcp_.Get(l);
}

template <class EB_>
auto Simple<EB_>::GetVolumeFlux(Step l) const -> const FieldEmbed<Scal>& {
  return imp->fev_.Get(l);
}

template <class EB_>
double Simple<EB_>::GetAutoTimeStep() const {
  return imp->GetAutoTimeStep();
}

template <class EB_>
double Simple<EB_>::GetError() const {
  return imp->cd_->GetError();
}

template <class EB_>
auto Simple<EB_>::GetVelocityCond() const -> const MapEmbed<BCond<Vect>>& {
  return imp->me_vel_;
}
