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
  // Expression on face: v[0] * cm + v[1] * cp + v[2]
  using ExprFace = generic::Vect<Scal, 3>;
  // Expression on cell: v[0] * c + v[1] * cxm + ... + v[6] * czp + v[7]
  using Expr = generic::Vect<Scal, M::dim * 2 + 2>;
  using UEB = UEmbed<M>;

  Imp(Owner* owner, const EB& eb0, const FieldCell<Vect>& fcw,
      const MapEmbed<BCondFluid<Vect>>& mebc,
      const MapCell<std::shared_ptr<CondCellFluid>>& mcc, Par par,
      const FieldFace<Scal>* ffbp)
      : owner_(owner)
      , par(par)
      , m(owner_->m)
      , eb(eb0)
      , dr_(0, m.GetEdim())
      , drr_(m.GetEdim(), dim)
      , mebc_(mebc)
      , mcc_(mcc)
      , ffbp_(ffbp)
      , fcpcs_(m)
      , ffvc_(m) {
    using namespace fluid_condition;

    UpdateDerivedConditions();

    fcfcd_.Reinit(m, Vect(0));
    typename CD::Par p;
    SetConvDiffPar(p, par);
    cd_ = GetConvDiff<EB>()(
        par.conv, m, eb, fcw, mebc_vel_, owner_->fcr_, &ffd_, &fcfcd_,
        &fev_.iter_prev, owner_->GetTime(), owner_->GetTimeStep(), p);

    fcp_.time_curr.Reinit(m, 0.);
    fcp_.time_prev = fcp_.time_curr;

    // Calc initial volume fluxes
    auto ffwe = UEB::Interpolate(cd_->GetVelocity(), mebc_vel_, eb);
    fev_.time_curr.Reinit(m, 0.);
    for (auto f : eb.Faces()) {
      fev_.time_curr[f] = ffwe[f].dot(eb.GetSurface(f));
    }
    // Apply meshvel
    const Vect& meshvel = par.meshvel;
    for (auto f : eb.Faces()) {
      fev_.time_curr[f] -= meshvel.dot(eb.GetSurface(f));
    }

    fev_.time_prev = fev_.time_curr;
  }

  void UpdateDerivedConditions() {
    using namespace fluid_condition;

    mebc_vel_ = GetVelCond<M>(mebc_);
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
      const CondCellFluid* cb = it.second.get(); // cond base

      if (auto cd = dynamic_cast<const GivenPressure<M>*>(cb)) {
        mccp_[c] = std::make_shared<CondCellValFixed<Scal>>(cd->GetPressure());
      } else {
        throw std::runtime_error("proj: unknown cell condition");
      }
    }

    ffbd_.Reinit(m, false);
    for (auto p : mebc_.GetMapFace()) {
      ffbd_[p.first] = true;
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
        Scal pc = cd->second(); // new value for p[c]
        e = Expr(0);
        // override target cell
        e[0] = 1.;
        e[Expr::dim - 1] = pc;
        // override neighbours
        for (size_t q : eb.Nci(c)) {
          IdxCell cn = m.GetCell(c, q);
          auto& en = fcs[cn];
          size_t qn = (q % 2 == 0 ? q + 1 : q - 1); // id of c from cn
          en[Expr::dim - 1] += en[1 + qn] * pc;
          en[1 + qn] = 0;
        }
      }
    }
  }
  // Flux expressions in terms of pressure correction pc:
  //   /  grad(pc) * area / k + v, inner
  //   \  a, boundary
  // ffk: diag coeff [i]
  // ffv: constant term [i]
  // Returns:
  // ffe [i], v=ffe[c] stores expression: v[0]*cm + v[1]*cp + v[2]
  template <int dummy>
  FieldFace<ExprFace> CalcFlux(
      const FieldFace<Scal>& ffk, const FieldEmbed<Scal>& ffv) {
    FieldFace<ExprFace> ffe(m, ExprFace(0));
    for (auto f : eb.Faces()) {
      auto& e = ffe[f];
      const Scal hinv = 1. / m.GetCellSize()[0];
      if (!ffbd_[f]) { // inner
        const Scal a = -m.GetArea(f) * hinv / ffk[f];
        e[0] = -a;
        e[1] = a;
      } else { // boundary
        e[0] = 0;
        e[1] = 0;
      }
      e[2] = ffv[f];
    }
    return ffe;
  }
  template <int dummy>
  FieldEmbed<ExprFace> CalcFlux(
      const FieldEmbed<Scal>& fek, const FieldEmbed<Scal>& fev) {
    FieldEmbed<ExprFace> ffe(m, ExprFace(0));
    for (auto f : eb.Faces()) {
      auto& e = ffe[f];
      if (!ffbd_[f]) { // inner
        const IdxCell cm = m.GetCell(f, 0);
        const IdxCell cp = m.GetCell(f, 1);
        const Scal dn = eb.ClipGradDenom(
            eb.GetNormal(f).dot(eb.GetCellCenter(cp) - eb.GetCellCenter(cm)));
        const Scal a = -eb.GetArea(f) / (fek[f] * dn);
        e[0] = -a;
        e[1] = a;
      } else { // boundary
        e[0] = 0;
        e[1] = 0;
      }
      e[2] = fev[f];
    }
    for (auto c : eb.CFaces()) {
      auto& e = ffe[c];
      e[2] = fev[c];
    }
    return ffe;
  }
  // Expressions for sum of fluxes and source:
  //   sum(v) - sv * vol
  // ffv: fluxes [i]
  // fcsv: volume source [i]
  // Output:
  // fce: result [i]
  FieldCell<Expr> CalcFluxSum(
      const FieldFaceb<ExprFace>& ffv, const FieldCell<Scal>& fcsv) {
    // initialize as diagonal system
    FieldCell<Expr> fce(m, Expr::GetUnit(0));
    for (auto c : eb.Cells()) {
      Expr e(0);
      for (auto q : eb.Nci(c)) {
        const IdxFace f = m.GetFace(c, q);
        const ExprFace v = ffv[f] * m.GetOutwardFactor(c, q);
        e[0] += v[1 - q % 2];
        e[1 + q] += v[q % 2];
        e[Expr::dim - 1] += v[2];
      }
      e[Expr::dim - 1] -= fcsv[c] * m.GetVolume(c);
      fce[c] = e;
    }
    return fce;
  }
  // Get diagcoeff from current convdiff equations
  void GetDiagCoeff(FieldCell<Scal>& fck, FieldFaceb<Scal>& ffk) {
    auto sem = m.GetSem("diag");
    if (sem("local")) {
      fck.Reinit(m, 0);
      for (auto d : dr_) {
        auto fct = cd_->GetDiag(d);
        for (auto c : eb.Cells()) {
          fck[c] += fct[c];
        }
      }
      for (auto c : eb.Cells()) {
        fck[c] /= dr_.size();
      }

      CHECKNAN(fck, m.CN());

      m.Comm(&fck);
    }
    if (sem("interp")) {
      ffk = UEB::Interpolate(fck, {}, eb);
    }
  }
  // Append explicit part of viscous force.
  // fcw: velocity [a]
  // Output:
  // fcf += viscous term [i]
  void AppendExplViscous(
      const FieldCell<Vect>& fcw, FieldCell<Vect>& fcf, const M& m) {
    const auto wf = UEB::Interpolate(fcw, mebc_vel_, m);
    for (auto d : dr_) {
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
      const FieldCell<Vect>& fcw, FieldCell<Vect>& fcf, const Embed<M>& eb) {
    // FIXME: not implemented
    (void)fcw;
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
          m, m, GetVelocity(Step::iter_curr), mebc_, par.inletflux_numid,
          mebc_vel_);
    }
    if (sem.Nested("bc-outlet")) {
      UFluid<M>::template UpdateOutletVelocity<M>(
          m, m, GetVelocity(Step::iter_curr), mebc_, *owner_->fcsv_, mebc_vel_);
    }
  }
  // TODO: rewrite norm() using dist() where needed
  void MakeIteration() {
    auto sem = m.GetSem("fluid-iter");
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
      AppendExplViscous(cd_->GetVelocity(Step::iter_curr), fcfcd_, eb);
    }

    if (sem.Nested("convdiff-iter")) {
      // Convection-diffusion
      cd_->MakeIteration();
    }

    if (sem.Nested("diag")) {
      GetDiagCoeff(fck_, ffk_);
    }

    if (sem("pcorr-assemble")) {
      // Acceleration
      const auto ffvel =
          UEB::Interpolate(cd_->GetVelocity(Step::iter_curr), mebc_vel_, eb);
      auto& ffbp = *ffbp_;
      for (auto f : eb.Faces()) {
        Scal v = ffvel[f].dot(eb.GetSurface(f));
        if (!ffbd_[f]) { // inner
          v += ffbp[f] * eb.GetArea(f) / ffk_[f];
        } else { // boundary
          // nop, keep the mean flux
        }
        fev_.iter_curr[f] = v;
      }

      // Projection
      ffvc_ = CalcFlux<0>(ffk_, fev_.iter_curr);
      fcpcs_ = CalcFluxSum(ffvc_, *owner_->fcsv_);
      ApplyCellCond(fcp_curr, fcpcs_);
    }
    if (sem.Nested("pcorr-solve")) {
      Solve(fcpcs_, &fcp_curr, fcpc_, M::LS::T::symm, m);
    }
    if (sem("pcorr-apply")) {
      if (par.linreport && m.IsRoot()) {
        std::cout << "pcorr:"
                  << " res=" << m.GetResidual() << " iter=" << m.GetIter()
                  << std::endl;
      }
      // set pressure
      fcp_curr = fcpc_;

      for (auto f : eb.Faces()) {
        IdxCell cm = m.GetCell(f, 0);
        IdxCell cp = m.GetCell(f, 1);
        auto& e = ffvc_[f];
        fev_.iter_curr[f] = e[0] * fcpc_[cm] + e[1] * fcpc_[cp] + e[2];
      }

      // Acceleration and correction to center velocity
      // XXX adhoc , using me_visc_ but should be zero-derivative
      const auto ffgp = UEB::Gradient(fcp_curr, me_pressure_, eb);
      fcwc_.Reinit(m);
      auto& ffbp = *ffbp_;
      for (auto c : m.Cells()) {
        Vect s(0);
        for (auto q : eb.Nci(c)) {
          const IdxFace f = m.GetFace(c, q);
          if (!ffbd_[f]) { // inner
            const Scal a = (ffbp[f] - ffgp[f]) / ffk_[f];
            s += eb.GetNormal(f) * (a * 0.5);
          } else {
            // nop, no acceleration
          }
        }
        fcwc_[c] = s;
      }
    }

    if (sem.Nested("convdiff-corr")) {
      cd_->CorrectVelocity(Step::iter_curr, fcwc_);
    }

    if (sem("inc-iter")) {
      owner_->IncIter();
      fcwc_.Free();
      fck_.Free();
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
  GRange<size_t> dr_; // effective dimension range
  GRange<size_t> drr_; // remaining dimensions

  // Face conditions
  const MapEmbed<BCondFluid<Vect>>& mebc_;
  MapEmbed<BCond<Vect>> mebc_vel_;
  MapEmbed<BCond<Scal>> me_pressure_;
  MapEmbed<BCond<Vect>> me_force_;
  MapEmbed<BCond<Scal>> me_visc_;

  // Cell conditions
  MapCell<std::shared_ptr<CondCellFluid>> mcc_; // fluid cell cond
  MapCell<std::shared_ptr<CondCell>> mccp_; // pressure cell cond

  const FieldFace<Scal>* ffbp_;
  StepData<FieldEmbed<Scal>> fev_; // volume flux
  StepData<FieldCell<Scal>> fcp_; // pressure

  std::shared_ptr<CD> cd_;

  // TODO: Const specifier for CondFace*

  FieldFace<bool> ffbd_; // is boundary

  // Cell fields:
  FieldCell<Scal> fck_; // diag coeff of velocity equation
  FieldCell<Expr> fcpcs_; // pressure correction linear system [i]
  FieldCell<Scal> fcpc_; // pressure correction
  FieldCell<Vect> fcwc_; // velocity correction
  FieldCell<Vect> fcfcd_; // force for convdiff [i]

  // tmp
  FieldCell<Scal> fct_;
  FieldCell<Vect> fctv_;

  // Face fields:
  FieldFaceb<Scal> ffd_; // dynamic viscosity
  FieldFaceb<Scal> ffk_; // diag coeff of velocity equation
  FieldFaceb<ExprFace> ffvc_; // expression for corrected volume flux [i]
};

template <class EB_>
Proj<EB_>::Proj(
    M& m, const EB& eb, const FieldCell<Vect>& fcw,
    const MapEmbed<BCondFluid<Vect>>& mebc,
    const MapCell<std::shared_ptr<CondCellFluid>>& mcc, FieldCell<Scal>* fcr,
    FieldCell<Scal>* fcd, FieldCell<Vect>* fcf, FieldFace<Scal>* ffbp,
    FieldCell<Scal>* fcsv, FieldCell<Scal>* fcsm, double t, double dt, Par par)
    : FluidSolver<M>(t, dt, m, fcr, fcd, fcf, nullptr, fcsv, fcsm)
    , imp(new Imp(this, eb, fcw, mebc, mcc, par, ffbp)) {}

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
  return imp->mebc_vel_;
}
