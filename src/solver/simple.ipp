// Created by Petr Karnakov on 31.07.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <cmath>
#include <memory>
#include <sstream>
#include <stdexcept>

#include "approx.h"
#include "convdiffv.h"
#include "debug/isnan.h"
#include "fluid.h"
#include "simple.h"
#include "util/convdiff.h"
#include "util/fluid.h"
#include "util/metrics.h"

// Rules:
// - Each function assumes that all fields on layers
//     iter_prev, time_prev, time_curr
//   and force, source and viscosity
//   are known and doesn't modify them.
// - No function except for MakeIteration refers to iter_curr.
//
// domain (cells/faces)
// [i]: inner
// [s]: support
// [a]: all
//
// notation:
// p: pressure
// gp: pressure gradient
// w: velocity
// v: volume flux
// we: predicted velocity (after solving velocity equations)
// ve: predicted volume flux

template <class M_>
struct Simple<M_>::Imp {
  using Owner = Simple<M_>;
  using CD = ConvDiffVect<M>;
  // Expression on face: v[0] * cm + v[1] * cp + v[2]
  using ExprFace = generic::Vect<Scal, 3>;
  // Expression on cell: v[0] * c + v[1] * cxm + ... + v[6] * czp + v[7]
  using Expr = generic::Vect<Scal, M::dim * 2 + 2>;

  Imp(Owner* owner, const FieldCell<Vect>& fcw,
      const MapEmbed<BCondFluid<Vect>>& mebc, MapCondFaceFluid& mfc,
      const MapCell<std::shared_ptr<CondCellFluid>>& mcc, Par par)
      : owner_(owner)
      , par(par)
      , m(owner_->m)
      , dr_(0, m.GetEdim())
      , drr_(m.GetEdim(), dim)
      , mebc_(mebc)
      , mfc_(mfc)
      , mcc_(mcc)
      , fcpcs_(m)
      , ffvc_(m) {
    using namespace fluid_condition;

    mfc_ = GetCondFluid<M>(mebc_);

    ffbd_.Reinit(m, false);

    UpdateDerivedConditions();

    fcfcd_.Reinit(m, Vect(0));
    typename CD::Par p;
    SetConvDiffPar(p, par);
    cd_ = GetConvDiff<M>()(
        par.conv, m, m, fcw, mfcw_, owner_->fcr_, &ffd_, &fcfcd_,
        &ffv_.iter_prev, owner_->GetTime(), owner_->GetTimeStep(), p);

    fcp_.time_curr.Reinit(m, 0.);
    fcp_.time_prev = fcp_.time_curr;

    // Calc initial volume fluxes
    auto ffwe = Interpolate(cd_->GetVelocity(), mfcw_, m);
    ffv_.time_curr.Reinit(m, 0.);
    for (auto f : m.Faces()) {
      ffv_.time_curr[f] = ffwe[f].dot(m.GetSurface(f));
    }
    // Apply meshvel
    const Vect& meshvel = par.meshvel;
    for (auto f : m.Faces()) {
      ffv_.time_curr[f] -= meshvel.dot(m.GetSurface(f));
    }

    ffv_.time_prev = ffv_.time_curr;

    auto ffp = Interpolate(fcp_.time_curr, mfcp_, m);
    fcgp_ = Gradient(ffp, m);
  }

  void UpdateDerivedConditions() {
    using namespace fluid_condition;

    mfcw_ = GetVelCond(m, mfc_);
    mfcf_.clear();
    mfcp_.clear();
    mfcpc_.clear();
    mfcd_.clear();
    for (auto& it : mfc_) {
      const IdxFace f = it.first;
      ffbd_[f] = true;
      auto& cb = it.second;
      size_t nci = cb->GetNci();

      mfcf_[f].template Set<CondFaceGradFixed<Vect>>(Vect(0), nci);
      mfcp_[f].template Set<CondFaceExtrap>(nci);
      mfcpc_[f].template Set<CondFaceExtrap>(nci);
      mfcd_[f].template Set<CondFaceGradFixed<Scal>>(0., nci);

      if (cb.template Get<NoSlipWall<M>>()) {
        // nop
      } else if (cb.template Get<Inlet<M>>()) {
        // nop
      } else if (cb.template Get<Outlet<M>>()) {
        // nop
      } else if (cb.template Get<SlipWall<M>>() || cb.template Get<Symm<M>>()) {
        mfcf_[f].template Set<CondFaceReflect>(nci);
        mfcp_[f].template Set<CondFaceGradFixed<Scal>>(0., nci);
        mfcpc_[f].template Set<CondFaceGradFixed<Scal>>(0, nci);
      } else {
        throw std::runtime_error("proj: unknown condition");
      }
    }

    mccp_.clear();
    mccp_.clear();
    for (auto& it : mcc_) {
      const IdxCell c = it.first;
      CondCellFluid* cb = it.second.get(); // cond base

      if (auto cd = dynamic_cast<GivenPressure<M>*>(cb)) {
        mccp_[c] = std::make_shared<CondCellValFixed<Scal>>(cd->GetPressure());
      } else {
        throw std::runtime_error("proj: unknown cell condition");
      }
    }
  }
  // Restore force from projections.
  // ffbp: force projections, bp=b.dot(n)
  // Output:
  // fcb: restored force
  void CalcExtForce(const FieldFace<Scal>& ffbp, FieldCell<Vect>& fcb) {
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
              m.GetSurface(f) * (ffbp[f] * m.GetVolume(c) / m.GetArea(f) * 0.5);
        }
        fcb[c] = s / m.GetVolume(c);
      }

      // TODO: comm ffbp_ instead
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
      ffv_.iter_curr = ffv_.time_curr;
    }
  }
  // Rhie-Chow interpolation of predicted volume flux
  // including balanced force (hydrostatics and surface tension)
  // fcw: predicted velocity field [s]
  // fcp: pressure field [s]
  // fcgp: gradient of pressure field [s]
  // fck, ffk: diag coeff [s]
  // Output:
  // ffv: result [i]
  void RhieChow(
      const FieldCell<Vect>& fcw, const FieldCell<Scal>& fcp,
      const FieldCell<Vect>& fcgp, const FieldCell<Scal>& fck,
      const FieldFace<Scal>& ffk, FieldFace<Scal>& ffv) {
    auto fftv = Interpolate(fcw, mfcw_, m);

    const Scal rh = par.rhie; // rhie factor
    ffv.Reinit(m);
    auto& ffbp = *owner_->ffbp_;
    for (auto f : m.Faces()) {
      // Init with mean flux
      ffv[f] = fftv[f].dot(m.GetSurface(f));
      if (!ffbd_[f]) { // if not boundary
        IdxCell cm = m.GetCell(f, 0);
        IdxCell cp = m.GetCell(f, 1);

        // compact pressure gradient
        Scal hr = m.GetArea(f) / m.GetVolume(cp);
        Scal gp = (fcp[cp] - fcp[cm]) * hr;

        // compact
        Scal o = (ffbp[f] - gp) * m.GetArea(f) / ffk[f];

        // wide
        Vect wm = (fcb_[cm] - fcgp[cm]) / fck[cm];
        Vect wp = (fcb_[cp] - fcgp[cp]) / fck[cp];
        Scal w = (wm + wp).dot(m.GetSurface(f)) * 0.5;

        // apply
        ffv[f] += rh * (o - w);
      } else { // if boundary
        // nop, keep mean flux
      }
    }

    // Apply meshvel
    for (auto f : m.Faces()) {
      ffv[f] -= par.meshvel.dot(m.GetSurface(f));
    }
  }
  // Apply cell conditions for pressure.
  // fcpb: base pressure [i]
  // fcs: linear system in terms of correction of base pressure [i]
  void ApplyPcCond(const FieldCell<Scal>& fcpb, FieldCell<Expr>& fcs) {
    for (auto& it : mccp_) {
      const IdxCell c = it.first; // target cell
      CondCell* cb = it.second.get(); // cond base
      if (auto cd = dynamic_cast<CondCellVal<Scal>*>(cb)) {
        auto& e = fcs[c];
        Scal pc = cd->second() - fcpb[c]; // new value for p[c]
        e = Expr(0);
        // override target cell
        e[0] = 1.;
        e[Expr::dim - 1] = pc;
        // override neighbours
        for (size_t q : m.Nci(c)) {
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
  //   /  grad(pb+pc) * area / k + v, inner
  //   \  a, boundary
  // fcpb: base pressure [s]
  // ffk: diag coeff [i]
  // ffv: addition to flux [i]
  // Output:
  // ffe: result [i], Vect v stores expression: v[0]*cm + v[1]*cp + v[2]
  void GetFlux(
      const FieldCell<Scal>& fcpb, const FieldFace<Scal>& ffk,
      const FieldFace<Scal>& ffv, FieldFace<ExprFace>& ffe) {
    ffe.Reinit(m);
    for (auto f : m.Faces()) {
      auto& e = ffe[f];
      IdxCell cm = m.GetCell(f, 0);
      IdxCell cp = m.GetCell(f, 1);
      Scal hr = m.GetArea(f) / m.GetVolume(cp);
      if (!ffbd_[f]) { // inner
        Scal a = -m.GetArea(f) * hr / ffk[f];
        e[0] = -a;
        e[1] = a;
        e[2] = (fcpb[cp] - fcpb[cm]) * a + ffv[f];
      } else { // boundary
        e[0] = 0;
        e[1] = 0;
        e[2] = ffv[f];
      }
    }
  }
  // Flux expressions in terms of pressure correction pc:
  //   /  grad(pc) * area / k + v, inner
  //   \  a, boundary
  // ffk: diag coeff [i]
  // ffv: addition to flux [i]
  // Output:
  // ffe: result [i], Vect v stores expression: v[0]*cm + v[1]*cp + v[2]
  void GetFlux(
      const FieldFace<Scal>& ffk, const FieldFace<Scal>& ffv,
      FieldFace<ExprFace>& ffe) {
    ffe.Reinit(m);
    for (auto f : m.Faces()) {
      auto& e = ffe[f];
      IdxCell cp = m.GetCell(f, 1);
      Scal hr = m.GetArea(f) / m.GetVolume(cp);
      if (!ffbd_[f]) { // inner
        Scal a = -m.GetArea(f) * hr / ffk[f];
        e[0] = -a;
        e[1] = a;
      } else { // boundary
        e[0] = 0;
        e[1] = 0;
      }
      e[2] = ffv[f];
    }
  }
  // Expressions for sum of fluxes and source:
  //   sum(v) - sv * vol
  // ffv: fluxes [i]
  // fcsv: volume source [i]
  // Output:
  // fce: result [i]
  void GetFluxSum(
      const FieldFace<ExprFace>& ffv, const FieldCell<Scal>& fcsv,
      FieldCell<Expr>& fce) {
    fce.Reinit(m, Expr(0));
    for (auto c : m.Cells()) {
      auto& e = fce[c];
      for (auto q : m.Nci(c)) {
        IdxFace f = m.GetFace(c, q);
        ExprFace v = ffv[f] * m.GetOutwardFactor(c, q);
        e[0] += v[1 - q % 2];
        e[1 + q] += v[q % 2];
        e[Expr::dim - 1] += v[2];
      }
      e[Expr::dim - 1] -= fcsv[c] * m.GetVolume(c);
    }
  }
  // Expressions for sum of fluxes.
  //   sum(v)
  // ffv: fluxes [i]
  // Output:
  // fce: result [i]
  void GetFluxSum(const FieldFace<ExprFace>& ffv, FieldCell<Expr>& fce) {
    fce.Reinit(m, Expr(0));
    for (auto c : m.Cells()) {
      auto& e = fce[c];
      for (auto q : m.Nci(c)) {
        IdxFace f = m.GetFace(c, q);
        ExprFace v = ffv[f] * m.GetOutwardFactor(c, q);
        e[0] += v[1 - q % 2];
        e[1 + q] += v[q % 2];
        e[Expr::dim - 1] += v[2];
      }
    }
  }
  // Solve linear system fce = 0
  // fce: expressions [i]
  // Output:
  // fc: result [a]
  // m.GetSolveTmp(): modified temporary fields
  void Solve(const FieldCell<Expr>& fce, FieldCell<Scal>& fc) {
    auto sem = m.GetSem("solve");
    if (sem("solve")) {
      std::vector<Scal>* lsa;
      std::vector<Scal>* lsb;
      std::vector<Scal>* lsx;
      m.GetSolveTmp(lsa, lsb, lsx);
      lsx->resize(m.GetInBlockCells().size());
      size_t i = 0;
      for (auto c : m.Cells()) {
        (void)c;
        (*lsx)[i++] = 0;
      }
      auto l = ConvertLsCompact(fce, *lsa, *lsb, *lsx, m);
      using T = typename M::LS::T;
      l.t = T::symm; // solver type
      m.Solve(l);
    }
    if (sem("copy")) {
      std::vector<Scal>* lsa;
      std::vector<Scal>* lsb;
      std::vector<Scal>* lsx;
      m.GetSolveTmp(lsa, lsb, lsx);

      fc.Reinit(m);
      size_t i = 0;
      for (auto c : m.Cells()) {
        fc[c] = (*lsx)[i++];
      }
      CHECKNAN(fc, m.CN());
      m.Comm(&fc);
      if (par.linreport && m.IsRoot()) {
        std::cout << "pcorr:"
                  << " res=" << m.GetResidual() << " iter=" << m.GetIter()
                  << std::endl;
      }
    }
  }
  // Get diagcoeff from current convdiff equations
  void GetDiagCoeff(FieldCell<Scal>& fck, FieldFace<Scal>& ffk) {
    auto sem = m.GetSem("diag");
    if (sem("local")) {
      fck.Reinit(m, 0);
      for (auto d : dr_) {
        auto fct = cd_->GetDiag(d);
        for (auto c : m.Cells()) {
          fck[c] += fct[c];
        }
      }
      for (auto c : m.Cells()) {
        fck[c] /= dr_.size();
      }

      CHECKNAN(fck, m.CN())

      m.Comm(&fck);
    }
    if (sem("interp")) {
      ffk.Reinit(m);
      InterpolateI(fck, ffk, m);
    }
  }
  // Append explicit part of viscous force.
  // fcw: velocity [a]
  // Output:
  // fcf += viscous term [i]
  void AppendExplViscous(const FieldCell<Vect>& fcw, FieldCell<Vect>& fcf) {
    auto wf = Interpolate(fcw, mfcw_, m);
    for (auto d : dr_) {
      auto wfo = GetComponent(wf, d);
      auto gc = Gradient(wfo, m);
      auto gf = Interpolate(gc, mfcf_, m); // XXX adhoc zero-deriv cond
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
  // fcw: given velocity
  // ffv: given volume flux
  // fcpp: previous pressure
  // Output:
  // fcp: output pressure (may be aliased with fcpp)
  // fctv_: modified tmp
  void CalcPressure(
      const FieldCell<Vect>& fcw, const FieldFace<Scal>& ffv,
      const FieldCell<Scal>& fcpp, FieldCell<Scal>& fcp) {
    auto sem = m.GetSem("calcpressure");
    auto& fcl = fctv_; // evaluation of velocity equations
    if (sem.Nested("cd-asm")) {
      cd_->Assemble(fcw, ffv);
    }
    if (sem.Nested("diag")) {
      GetDiagCoeff(fck_, ffk_);
    }
    if (sem("eval")) {
      fcl.Reinit(m);
      for (auto d : dr_) {
        SetComponent(fcl, d, cd_->GetConst(d));
      }
      for (auto d : drr_) {
        SetComponent(fcl, d, 0);
      }
      m.Comm(&fcl);
    }
    if (sem("assemble")) {
      const Scal rh = par.rhie; // rhie factor
      FieldFace<Scal> ffa(m); // addition to flux TODO revise comment
      auto& ffbp = *owner_->ffbp_;
      for (auto f : m.Faces()) {
        IdxCell cm = m.GetCell(f, 0);
        IdxCell cp = m.GetCell(f, 1);

        if (!ffbd_[f]) { // if not boundary
          auto s = m.GetSurface(f);
          auto sa = m.GetArea(f);
          auto kf = rh * sa / ffk_[f];
          Vect bm = fcw[cm] - (fcl[cm] - fcgp_[cm] + fcb_[cm]) / fck_[cm] * rh;
          Vect bp = fcw[cp] - (fcl[cp] - fcgp_[cp] + fcb_[cp]) / fck_[cp] * rh;
          ffa[f] = (bm + bp).dot(s) * 0.5 + ffbp[f] * kf - ffv[f];
        } else { // if boundary
          ffa[f] = 0.;
        }
      }
      fcl.Free();

      GetFlux(fcpp, ffk_, ffa, ffvc_);

      GetFluxSum(ffvc_, fcpcs_);

      ApplyPcCond(fcpp, fcpcs_);
    }

    if (sem.Nested("pcorr-solve")) {
      Solve(fcpcs_, fcpc_);
    }

    if (sem("apply")) {
      fcgpc_ = Gradient(Interpolate(fcpc_, mfcpc_, m), m);

      // Correct pressure
      Scal pr = par.prelax; // pressure relaxation
      for (auto c : m.Cells()) {
        fcp[c] = fcpp[c] + fcpc_[c] * pr;
      }
      m.Comm(&fcp);
    }
  }
  void UpdateBc(typename M::Sem& sem) {
    if (sem.Nested("bc-inletflux")) {
      UFluid<M>::UpdateInletFlux(
          m, GetVelocity(Step::iter_curr), mfc_, par.inletflux_numid);
    }
    if (sem.Nested("bc-outlet")) {
      UFluid<M>::UpdateOutletBaseConditions(
          m, GetVelocity(Step::iter_curr), mfc_, *owner_->fcsv_);
    }
    if (sem("bc-derived")) {
      UpdateDerivedConditions();
    }
  }
  void CalcForce(typename M::Sem& sem) {
    if (sem.Nested("forceinit")) {
      CalcExtForce(*owner_->ffbp_, fcb_);
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
      ffd_ = Interpolate(*owner_->fcd_, mfcd_, m);

      // rotate layers
      fcp_prev = fcp_curr;
      ffv_.iter_prev = ffv_.iter_curr;
    }

    UpdateBc(sem);

    CalcForce(sem);

    if (sem("pgrad")) {
      auto ffp = Interpolate(fcp_curr, mfcp_, m);
      fcgp_ = Gradient(ffp, m);

      // append pressure gradient to force
      for (auto c : m.Cells()) {
        fcfcd_[c] += fcgp_[c] * (-1.);
      }
    }

    if (par.simpler) {
      if (sem.Nested("simpler")) {
        CalcPressure(
            cd_->GetVelocity(Step::iter_curr), ffv_.iter_curr.GetFieldFace(),
            fcp_curr, fcp_curr);
      }

      if (sem("pgrad")) {
        auto ffp = Interpolate(fcp_curr, mfcp_, m);
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

      GetFlux(ffk_, ffve_, ffvc_);
      CHECKNAN(ffvc_, m.CN())

      GetFluxSum(ffvc_, *owner_->fcsv_, fcpcs_);
      CHECKNAN(fcpcs_, m.CN())

      ApplyPcCond(fcp_curr, fcpcs_);
    }

    if (sem.Nested("pcorr-solve")) {
      Solve(fcpcs_, fcpc_);
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
        ffv_.iter_curr[f] = e[0] * fcpc_[cm] + e[1] * fcpc_[cp] + e[2];
      }
      CHECKNAN(ffv_.iter_curr, m.CN())

      auto ffpc = Interpolate(fcpc_, mfcpc_, m);
      CHECKNAN(ffpc, m.CN())
      fcgpc_ = Gradient(ffpc, m);
      CHECKNAN(fcgpc_, m.CN())

      // Calc velocity correction
      fcwc_.Reinit(m);
      for (auto c : m.Cells()) {
        fcwc_[c] = fcgpc_[c] / (-fck_[c]);
      }
      CHECKNAN(fcwc_, m.CN())
    }

    if (sem.Nested("convdiff-corr")) {
      // Correct velocity and comm
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
      ffv_.time_prev = ffv_.time_curr;
      fcp_.time_curr = fcp_.iter_curr;
      ffv_.time_curr = ffv_.iter_curr;
      CHECKNAN(fcp_.time_curr, m.CN())
      owner_->IncTime();
    }
    if (sem.Nested("convdiff-finish")) {
      cd_->FinishStep();
    }
  }
  double GetAutoTimeStep() {
    double dt = 1e10;
    auto& flux = ffv_.time_curr;
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
  M& m; // mesh
  GRange<size_t> dr_; // effective dimension range
  GRange<size_t> drr_; // remaining dimensions

  // Face conditions
  const MapEmbed<BCondFluid<Vect>>& mebc_;
  MapCondFaceFluid& mfc_; // fluid cond
  MapCondFace mfcw_; // velocity cond
  MapCondFace mfcp_; // pressure cond
  MapCondFace mfcf_; // force cond
  MapCondFace mfcpc_; // pressure corr cond
  MapCondFace mfcd_; // dynamic viscosity cond

  // Cell conditions
  MapCell<std::shared_ptr<CondCellFluid>> mcc_; // fluid cell cond
  MapCell<std::shared_ptr<CondCell>> mccp_; // pressure cell cond

  StepData<FieldEmbed<Scal>> ffv_; // volume flux
  StepData<FieldCell<Scal>> fcp_; // pressure

  std::unique_ptr<CD> cd_;

  // TODO: Const specifier for CondFace*

  FieldFace<bool> ffbd_; // is boundary

  // Cell fields:
  FieldCell<Vect> fcgp_; // gradient of pressure
  FieldCell<Scal> fck_; // diag coeff of velocity equation
  FieldCell<Expr> fcpcs_; // pressure correction linear system [i]
  FieldCell<Scal> fcpc_; // pressure correction
  FieldCell<Vect> fcgpc_; // gradient of pressure correction
  FieldCell<Vect> fcwc_; // velocity correction
  FieldCell<Vect> fcb_; // restored balanced force [s]
  FieldCell<Vect> fcfcd_; // force for convdiff [i]

  // tmp
  FieldCell<Scal> fct_;
  FieldCell<Vect> fctv_;

  // Face fields:
  FieldFace<Scal> ffd_; // dynamic viscosity
  FieldFace<Scal> ffve_; // predicted volume flux [i]
  FieldFace<Scal> ffk_; // diag coeff of velocity equation
  FieldFace<ExprFace> ffvc_; // expression for corrected volume flux [i]
};

template <class M_>
Simple<M_>::Simple(
    M& m, const FieldCell<Vect>& fcw, const MapEmbed<BCondFluid<Vect>>& mebc,
    MapCondFaceFluid& mfc, const MapCell<std::shared_ptr<CondCellFluid>>& mcc,
    FieldCell<Scal>* fcr, FieldCell<Scal>* fcd, FieldCell<Vect>* fcf,
    FieldFace<Scal>* ffbp, FieldCell<Scal>* fcsv, FieldCell<Scal>* fcsm,
    double t, double dt, Par par)
    : FluidSolver<M>(t, dt, m, fcr, fcd, fcf, ffbp, fcsv, fcsm)
    , imp(new Imp(this, fcw, mebc, mfc, mcc, par)) {}

template <class M_>
Simple<M_>::~Simple() = default;

template <class M_>
auto Simple<M_>::GetPar() const -> const Par& {
  return imp->par;
}

template <class M_>
void Simple<M_>::SetPar(Par par) {
  imp->par = par;
}

template <class M_>
void Simple<M_>::StartStep() {
  return imp->StartStep();
}

template <class M_>
void Simple<M_>::MakeIteration() {
  return imp->MakeIteration();
}

template <class M_>
void Simple<M_>::FinishStep() {
  return imp->FinishStep();
}

template <class M_>
auto Simple<M_>::GetVelocity(Step l) const -> const FieldCell<Vect>& {
  return imp->GetVelocity(l);
}

template <class M_>
auto Simple<M_>::GetPressure(Step l) const -> const FieldCell<Scal>& {
  return imp->fcp_.Get(l);
}

template <class M_>
auto Simple<M_>::GetVolumeFlux(Step l) const -> const FieldEmbed<Scal>& {
  return imp->ffv_.Get(l);
}

template <class M_>
double Simple<M_>::GetAutoTimeStep() const {
  return imp->GetAutoTimeStep();
}

template <class M_>
double Simple<M_>::GetError() const {
  return imp->cd_->GetError();
}

template <class M_>
auto Simple<M_>::GetVelocityCond() const -> const MapCondFace& {
  return imp->mfcw_;
}
