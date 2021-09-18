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

  Imp(Owner* owner, const EB& eb_, const Args& args)
      : owner_(owner)
      , par(args.par)
      , m(owner_->m)
      , eb(eb_)
      , edim_range_(0, m.GetEdim())
      , drr_(m.GetEdim(), dim)
      , linsolver_symm_(args.linsolver_symm)
      , linsolver_gen_(args.linsolver_gen)
      , mebc_(args.mebc)
      , mcc_(args.mcc) {
    UpdateDerivedConditions();

    fcfcd_.Reinit(m, Vect(0));
    typename CD::Par p;
    SetConvDiffPar(p, par);

    fassert(linsolver_symm_);
    fassert(linsolver_gen_);

    const auto& ffv = fev_.iter_prev.template Get<FieldFaceb<Scal>>();
    const ConvDiffVectArgs<EB> cdargs{
        args.fcvel,     me_vel_, owner_->fcr_,      &ffvisc_,
        &fcfcd_,        &ffv,    owner_->GetTime(), owner_->GetTimeStep(),
        linsolver_gen_, p};

    cd_ = GetConvDiff<EB>()(par.conv, m, eb, cdargs);

    fcp_.time_curr.Reinit(m, 0.);
    fcp_.time_prev = fcp_.time_curr;

    // Calc initial volume fluxes
    const auto ffwe = UEB::Interpolate(cd_->GetVelocity(), me_vel_, eb);
    fev_.time_curr.Reinit(m, 0.);
    eb.LoopFaces([&](auto cf) { //
      fev_.time_curr[cf] = (ffwe[cf] - par.meshvel).dot(eb.GetSurface(cf));
    });
    fev_.time_prev = fev_.time_curr;
  }

  void UpdateDerivedConditions() {
    me_vel_ = ConvertBCondFluidToVelocity<M>(mebc_);
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

      using namespace fluid_condition;
      if (auto cd = dynamic_cast<const GivenPressure<M>*>(cb)) {
        mc_pressure_[c] =
            std::make_shared<CondCellValFixed<Scal>>(cd->GetPressure());
      } else {
        fassert(false, "unknown cell condition");
      }
    }

    is_boundary_.Reinit(m, false);
    mebc_.LoopPairs([&](auto p) { //
      is_boundary_[p.first] = true;
    });

    ffvisc_ = UEB::Interpolate(*owner_->fcd_, me_visc_, eb);
  }
  void StartStep() {
    auto sem = m.GetSem("fluid-start");
    if (sem("convdiff-init")) {
      owner_->ClearIter();
      CHECKNAN(fcp_.time_curr, m.flags.check_nan)
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
    auto& fcbp = fc_bforce_;

    for (auto f : eb.Faces()) {
      // mean flux
      auto& v = fev[f];
      v = fftv[f].dot(eb.GetSurface(f));
      if (!is_boundary_[f]) {
        const IdxCell cm = eb.GetCell(f, 0);
        const IdxCell cp = eb.GetCell(f, 1);

        // compact pressure gradient
        const Scal gp = (fcp[cp] - fcp[cm]) / eb.GetCellSize()[0];

        // compact
        const Scal o = (febp[f] - gp) * eb.GetArea(f) / fek[f];

        // wide
        const Vect wm = (fcbp[cm] - fcgp[cm]) / fck[cm];
        const Vect wp = (fcbp[cp] - fcgp[cp]) / fck[cp];
        const Scal w = (wm + wp).dot(eb.GetSurface(f)) * 0.5;

        // apply
        v += rh * (o - w);
      } else { // boundary
        // nop, keep mean flux
      }
    }

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
        e.back() -= pc;
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
    FieldFaceb<ExprFace> ffe = UEB::GradientImplicit({}, eb);
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
      sum.back() -= fcsv[c] * eb.GetVolume(c);
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

      CHECKNAN(fck, m.flags.check_nan)

      m.Comm(&fck);
    }
    if (sem("interp")) {
      ffk = UEB::InterpolateHarmonic(fck, me_visc_, eb);
    }
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
          par.outlet_relax, me_vel_);
    }
  }
  void CalcForce(typename M::Sem& sem) {
    if (sem("forceinit")) {
      fc_bforce_ = UEB::AverageGradient(
          owner_->febp_->template Get<FieldFaceb<Scal>>(), eb);
      m.Comm(&fc_bforce_);
    }

    if (sem("forceinit")) {
      // initialize force for convdiff
      fcfcd_.Reinit(m, Vect(0));
      // append force and balanced force
      auto& fcf = *owner_->fcf_;
      for (auto c : eb.Cells()) {
        fcfcd_[c] += fcf[c] + fc_bforce_[c];
      }
    }
    if (par.explviscous && sem("explvisc")) {
      UFluid<M>::AppendExplViscous(
          fcfcd_, cd_->GetVelocity(Step::iter_curr), me_vel_, ffvisc_, eb);
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
      CHECKNAN(fev_rhie, m.flags.check_nan)

      ffvc = GetFlux(ffk, fev_rhie);
      CHECKNAN(ffvc, m.flags.check_nan)

      fcpcs = GetFluxSum(ffvc, *owner_->fcsv_);
      CHECKNAN(fcpcs, m.flags.check_nan)

      ApplyCellCond(fcp_curr, fcpcs);
      fcpcs.SetName("pressure");
    }

    if (sem.Nested("pcorr-solve")) {
      linsolver_symm_->Solve(fcpcs, nullptr, fcpc, m);
    }
    if (sem("pcorr-apply")) {
      CHECKNAN(fcpc, m.flags.check_nan)

      // Correct pressure
      Scal pr = par.prelax; // pressure relaxation
      for (auto c : eb.AllCells()) {
        fcp_curr[c] += pr * fcpc[c];
      }

      // Calc divergence-free volume fluxes
      eb.LoopFaces([&](auto cf) {
        fev_.iter_curr[cf] = UEB::Eval(ffvc[cf], cf, fcpc, eb);
      });
      CHECKNAN(fev_.iter_curr, m.flags.check_nan)

      const auto fcgpc =
          UEB::AverageGradient(UEB::Gradient(fcpc, me_pcorr_, eb), eb);

      // Calc velocity correction
      fcvel_corr.Reinit(m);
      for (auto c : eb.Cells()) {
        fcvel_corr[c] = -fcgpc[c] / fck[c];
      }
      CHECKNAN(fcvel_corr, m.flags.check_nan)
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
      CHECKNAN(fcp_.time_curr, m.flags.check_nan)
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
  std::shared_ptr<linear::Solver<M>> linsolver_symm_;
  std::shared_ptr<linear::Solver<M>> linsolver_gen_;

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
  FieldCell<Vect> fc_bforce_; // restored balanced force [s]
  FieldCell<Vect> fcfcd_; // force for convdiff [i]

  // Face fields:
  FieldFaceb<Scal> ffvisc_; // dynamic viscosity
};

template <class EB_>
Simple<EB_>::Simple(M& m_, const EB& eb, const Args& args)
    : Base(
          args.t, args.dt, m_, args.fcr, args.fcd, args.fcf, args.febp,
          args.fcsv, args.fcsm)
    , imp(new Imp(this, eb, args)) {}

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
