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
      , mcc_(mcc)
      , fcpcs_(m)
      , ffvc_(m) {
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

    const auto ffp = UEB::Interpolate(fcp_.time_curr, me_pressure_, eb);
    fcgp_ = UEB::Gradient(ffp, eb);
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
  // Restore force from projections.
  // febp: force projections, bp=b.dot(n)
  // Output:
  // fcb: restored force
  void CalcExtForce(const FieldEmbed<Scal>& febp, FieldCell<Vect>& fcb) {
    auto sem = m.GetSem("extforce");

    if (sem("loc")) {
      // XXX specific for Cartesian mesh
      // TODO consider just a weighted average of fn * n
      //      weight should be proportional to accuracy gradient approx
      //      which is better if surface area is larger
      fcb.Reinit(m);
      for (auto c : m.Cells()) {
        Vect s(0);
        for (auto q : m.Nci(c)) {
          // TODO: revise for non-rectangular cell
          IdxFace f = m.GetFace(c, q);
          s +=
              m.GetSurface(f) * (febp[f] * m.GetVolume(c) / m.GetArea(f) * 0.5);
        }
        fcb[c] = s / m.GetVolume(c);
      }

      // TODO: comm febp_ instead
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
  void RhieChow(
      const FieldCell<Vect>& fcvel, const FieldCell<Scal>& fcp,
      const FieldCell<Vect>& fcgp, const FieldCell<Scal>& fck,
      const FieldFaceb<Scal>& fek, FieldFaceb<Scal>& fev) {
    auto fftv = UEB::Interpolate(fcvel, me_vel_, eb);

    const Scal rh = par.rhie; // rhie factor
    fev.Reinit(m);
    auto& febp = *owner_->febp_;
    eb.LoopFaces([&](auto cf) { //
      // mean flux
      auto& v = fev[cf];
      v = fftv[cf].dot(eb.GetSurface(cf));
      if (!is_boundary_[cf]) { // if not boundary
        const IdxCell cm = eb.GetCell(cf, 0);
        const IdxCell cp = eb.GetCell(cf, 1);

        // compact pressure gradient
        Scal gp = (fcp[cp] - fcp[cm]) / eb.GetCellSize()[0];

        // compact
        Scal o = (febp[cf] - gp) * eb.GetArea(cf) / fek[cf];

        // wide
        Vect wm = (fcb_[cm] - fcgp[cm]) / fck[cm];
        Vect wp = (fcb_[cp] - fcgp[cp]) / fck[cp];
        Scal w = (wm + wp).dot(eb.GetSurface(cf)) * 0.5;

        // apply
        v += rh * (o - w);
      } else { // if boundary
        // nop, keep mean flux
      }
    });

    // Apply meshvel
    eb.LoopFaces([&](auto cf) { //
      fev[cf] -= par.meshvel.dot(eb.GetSurface(cf));
    });
  }
  // Apply cell conditions for pressure.
  // fcpb: base pressure [i]
  // fcs: linear system in terms of correction of base pressure [i]
  void ApplyCellCond(const FieldCell<Scal>& fcpb, FieldCell<Expr>& fcs) {
    (void)fcpb;
    for (auto& it : mccp_) {
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
  // Flux expressions in terms of pressure correction pc:
  //   /  grad(pb+pc) * area / k + v, inner
  //   \  a, boundary
  // fcpb: base pressure [s]
  // ffk: diag coeff [i]
  // fev: addition to flux [i]
  // Output:
  // ffe: result [i], Vect v stores expression: v[0]*cm + v[1]*cp + v[2]
  FieldFaceb<ExprFace> GetFlux(
      const FieldCell<Scal>& fcpb, const FieldFaceb<Scal>& ffk,
      const FieldFaceb<Scal>& fev) {
    FieldFaceb<ExprFace> ffe(m);
    for (auto f : m.Faces()) {
      auto& e = ffe[f];
      IdxCell cm = m.GetCell(f, 0);
      IdxCell cp = m.GetCell(f, 1);
      Scal hr = m.GetArea(f) / m.GetVolume(cp);
      if (!is_boundary_[f]) { // inner
        Scal a = -m.GetArea(f) * hr / ffk[f];
        e[0] = -a;
        e[1] = a;
        e[2] = (fcpb[cp] - fcpb[cm]) * a + fev[f];
      } else { // boundary
        e[0] = 0;
        e[1] = 0;
        e[2] = fev[f];
      }
    }
    return ffe;
  }
  // Flux expressions in terms of pressure correction pc:
  //   /  grad(pc) * area / k + v, inner
  //   \  a, boundary
  // ffk: diag coeff [i]
  // fev: addition to flux [i]
  // Output:
  // ffe: result [i], Vect v stores expression: v[0]*cm + v[1]*cp + v[2]
  FieldFaceb<ExprFace> GetFlux(
      const FieldFaceb<Scal>& ffk, const FieldFaceb<Scal>& fev) {
    FieldFaceb<ExprFace> ffe(m);
    for (auto f : m.Faces()) {
      auto& e = ffe[f];
      IdxCell cp = m.GetCell(f, 1);
      Scal hr = m.GetArea(f) / m.GetVolume(cp);
      if (!is_boundary_[f]) { // inner
        Scal a = -m.GetArea(f) * hr / ffk[f];
        e[0] = -a;
        e[1] = a;
      } else { // boundary
        e[0] = 0;
        e[1] = 0;
      }
      e[2] = fev[f];
    }
    return ffe;
  }
  // Expressions for sum of fluxes and source:
  //   sum(v) - sv * vol
  // fev: fluxes [i]
  // fcsv: volume source [i]
  // Output:
  // fce: result [i]
  FieldCell<Expr> GetFluxSum(
      const FieldFaceb<ExprFace>& fev, const FieldCell<Scal>& fcsv) {
    FieldCell<Expr> fce(m, Expr(0));
    for (auto c : m.Cells()) {
      auto& e = fce[c];
      for (auto q : m.Nci(c)) {
        const IdxFace f = m.GetFace(c, q);
        const ExprFace v = fev[f] * m.GetOutwardFactor(c, q);
        e[0] += v[1 - q % 2];
        e[1 + q] += v[q % 2];
        e[Expr::dim - 1] += v[2];
      }
      e[Expr::dim - 1] -= fcsv[c] * m.GetVolume(c);
    }
    return fce;
  }
  // Expressions for sum of fluxes.
  //   sum(v)
  // fev: fluxes [i]
  // Output:
  // fce: result [i]
  FieldCell<Expr> GetFluxSum(const FieldFaceb<ExprFace>& fev) {
    FieldCell<Expr> fce(m, Expr(0));
    for (auto c : m.Cells()) {
      auto& e = fce[c];
      for (auto q : m.Nci(c)) {
        IdxFace f = m.GetFace(c, q);
        ExprFace v = fev[f] * m.GetOutwardFactor(c, q);
        e[0] += v[1 - q % 2];
        e[1 + q] += v[q % 2];
        e[Expr::dim - 1] += v[2];
      }
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
        for (auto c : m.Cells()) {
          fck[c] += fct[c];
        }
      }
      for (auto c : m.Cells()) {
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
  void AppendExplViscous(const FieldCell<Vect>& fcvel, FieldCell<Vect>& fcf) {
    const auto wf = UEB::Interpolate(fcvel, me_vel_, m);
    for (auto d : edim_range_) {
      const auto wfo = GetComponent(wf, d);
      const auto gc = Gradient(wfo, m);
      const auto gf = UEB::Interpolate(gc, me_force_, m);
      for (auto c : m.Cells()) {
        Vect s(0);
        for (auto q : m.Nci(c)) {
          IdxFace f = m.GetFace(c, q);
          s += gf[f] * (ffd_[f] * m.GetOutwardSurface(c, q)[d]);
        }
        fcf[c] += s / m.GetVolume(c);
      }
    }
  }
  // Restore pressure given velocity and volume flux
  // Assume MakeIteration() was called for convdiff solver
  // fcvel: given velocity
  // fev: given volume flux
  // fcpp: previous pressure
  // Output:
  // fcp: output pressure (may be aliased with fcpp)
  // fctv_: modified tmp
  void CalcPressure(
      const FieldCell<Vect>& fcvel, const FieldEmbed<Scal>& fev,
      const FieldCell<Scal>& fcpp, FieldCell<Scal>& fcp) {
    auto sem = m.GetSem("calcpressure");
    auto& fcl = fctv_; // evaluation of velocity equations
    if (sem.Nested("cd-asm")) {
      cd_->Assemble(fcvel, fev.template Get<FieldFaceb<Scal>>());
    }
    if (sem.Nested("diag")) {
      GetDiagCoeff(fck_, ffk_);
    }
    if (sem("eval")) {
      fcl.Reinit(m);
      for (auto d : edim_range_) {
        SetComponent(fcl, d, cd_->GetConst(d));
      }
      for (auto d : drr_) {
        SetComponent(fcl, d, 0);
      }
      m.Comm(&fcl);
    }
    if (sem("assemble")) {
      const Scal rh = par.rhie; // rhie factor
      FieldFaceb<Scal> ffa(m); // addition to flux TODO revise comment
      auto& febp = *owner_->febp_;
      for (auto f : m.Faces()) {
        IdxCell cm = m.GetCell(f, 0);
        IdxCell cp = m.GetCell(f, 1);

        if (!is_boundary_[f]) { // if not boundary
          auto s = m.GetSurface(f);
          auto sa = m.GetArea(f);
          auto kf = rh * sa / ffk_[f];
          Vect bm =
              fcvel[cm] - (fcl[cm] - fcgp_[cm] + fcb_[cm]) / fck_[cm] * rh;
          Vect bp =
              fcvel[cp] - (fcl[cp] - fcgp_[cp] + fcb_[cp]) / fck_[cp] * rh;
          ffa[f] = (bm + bp).dot(s) * 0.5 + febp[f] * kf - fev[f];
        } else { // if boundary
          ffa[f] = 0.;
        }
      }
      fcl.Free();

      ffvc_ = GetFlux(fcpp, ffk_, ffa);
      fcpcs_ = GetFluxSum(ffvc_);
      ApplyCellCond(fcpp, fcpcs_);
    }

    if (sem.Nested("pcorr-solve")) {
      Solve(fcpcs_, nullptr, fcpc_, M::LS::T::symm, m);
    }

    if (sem("apply")) {
      fcgpc_ = Gradient(UEB::Interpolate(fcpc_, me_pcorr_, m), m);

      // Correct pressure
      Scal pr = par.prelax; // pressure relaxation
      for (auto c : m.Cells()) {
        fcp[c] = fcpp[c] + fcpc_[c] * pr;
      }
      m.Comm(&fcp);
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
      for (auto c : m.Cells()) {
        fcfcd_[c] += fcf[c] + fcb_[c];
      }
    }

    if (sem("explvisc")) {
      AppendExplViscous(cd_->GetVelocity(Step::iter_curr), fcfcd_);
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
      ffd_ = UEB::Interpolate(*owner_->fcd_, me_visc_, m);

      // rotate layers
      fcp_prev = fcp_curr;
      fev_.iter_prev = fev_.iter_curr;
    }

    UpdateBc(sem);

    CalcForce(sem);

    if (sem("pgrad")) {
      auto ffp = UEB::Interpolate(fcp_curr, me_pressure_, m);
      fcgp_ = Gradient(ffp, m);

      // append pressure gradient to force
      for (auto c : m.Cells()) {
        fcfcd_[c] += fcgp_[c] * (-1.);
      }
    }

    if (par.simpler) {
      if (sem.Nested("simpler")) {
        CalcPressure(
            cd_->GetVelocity(Step::iter_curr), fev_.iter_curr, fcp_curr,
            fcp_curr);
      }

      if (sem("pgrad")) {
        auto ffp = UEB::Interpolate(fcp_curr, me_pressure_, m);
        fcgp_ = Gradient(ffp, m);

        // append pressure correction gradient to force
        for (auto c : m.Cells()) {
          fcfcd_[c] += fcgpc_[c] * (-1.);
        }
      }
    }

    if (sem.Nested("convdiff-iter")) {
      // Solve for predictor velocity
      cd_->MakeIteration();
    }

    if (sem.Nested("diag")) {
      GetDiagCoeff(fck_, ffk_);
    }

    if (sem("pcorr-assemble")) {
      RhieChow(
          cd_->GetVelocity(Step::iter_curr), fcp_curr, fcgp_, fck_, ffk_,
          ffve_);
      CHECKNAN(ffve_, m.CN())

      ffvc_ = GetFlux(ffk_, ffve_);
      CHECKNAN(ffvc_, m.CN())

      fcpcs_ = GetFluxSum(ffvc_, *owner_->fcsv_);
      CHECKNAN(fcpcs_, m.CN())

      ApplyCellCond(fcp_curr, fcpcs_);
    }

    if (sem.Nested("pcorr-solve")) {
      Solve(fcpcs_, nullptr, fcpc_, M::LS::T::symm, m);
    }

    if (sem("pcorr-apply")) {
      CHECKNAN(fcpc_, m.CN())

      if (!par.simpler) {
        // Correct pressure
        Scal pr = par.prelax; // pressure relaxation
        for (auto c : m.AllCells()) {
          fcp_curr[c] += pr * fcpc_[c];
        }
      }

      // Calc divergence-free volume fluxes
      for (auto f : m.Faces()) {
        IdxCell cm = m.GetCell(f, 0);
        IdxCell cp = m.GetCell(f, 1);
        auto& e = ffvc_[f];
        fev_.iter_curr[f] = e[0] * fcpc_[cm] + e[1] * fcpc_[cp] + e[2];
      }
      CHECKNAN(fev_.iter_curr, m.CN())

      auto ffpc = UEB::Interpolate(fcpc_, me_pcorr_, m);
      CHECKNAN(ffpc, m.CN())
      fcgpc_ = Gradient(ffpc, m);
      CHECKNAN(fcgpc_, m.CN())

      // Calc velocity correction
      fcvelcorr_.Reinit(m);
      for (auto c : m.Cells()) {
        fcvelcorr_[c] = fcgpc_[c] / (-fck_[c]);
      }
      CHECKNAN(fcvelcorr_, m.CN())
    }

    if (sem.Nested("convdiff-corr")) {
      // Correct velocity and comm
      cd_->CorrectVelocity(Step::iter_curr, fcvelcorr_);
    }

    if (sem("inc-iter")) {
      owner_->IncIter();
      fcvelcorr_.Free();
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
    auto& flux = fev_.time_curr;
    for (auto c : m.Cells()) {
      for (size_t i = 0; i < m.GetNumFaces(c); ++i) {
        IdxFace f = m.GetFace(c, i);
        if (flux[f] != 0.) {
          dt = std::min<Scal>(dt, std::abs(m.GetVolume(c) / flux[f]));
        }
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
  MapCell<std::shared_ptr<CondCell>> mccp_; // pressure cell cond

  StepData<FieldEmbed<Scal>> fev_; // volume flux
  StepData<FieldCell<Scal>> fcp_; // pressure

  std::unique_ptr<CD> cd_;

  FieldEmbed<bool> is_boundary_; // true on faces with boundary conditions

  // Cell fields:
  FieldCell<Vect> fcgp_; // gradient of pressure
  FieldCell<Scal> fck_; // diag coeff of velocity equation
  FieldCell<Expr> fcpcs_; // pressure correction linear system [i]
  FieldCell<Scal> fcpc_; // pressure correction
  FieldCell<Vect> fcgpc_; // gradient of pressure correction
  FieldCell<Vect> fcvelcorr_; // velocity correction
  FieldCell<Vect> fcb_; // restored balanced force [s]
  FieldCell<Vect> fcfcd_; // force for convdiff [i]

  // tmp
  FieldCell<Scal> fct_;
  FieldCell<Vect> fctv_;

  // Face fields:
  FieldFaceb<Scal> ffd_; // dynamic viscosity
  FieldFaceb<Scal> ffve_; // predicted volume flux [i]
  FieldFaceb<Scal> ffk_; // diag coeff of velocity equation
  FieldFaceb<ExprFace> ffvc_; // expression for corrected volume flux [i]
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
