// Created by Petr Karnakov on 30.07.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <array>
#include <exception>
#include <fstream>
#include <limits>
#include <memory>

#include "approx.h"
#include "debug/isnan.h"
#include "geom/block.h"
#include "normal.h"
#include "reconst.h"
#include "trackerm.h"
#include "util/vof.h"

#include "vof.h"

template <class M_>
struct Vof<M_>::Imp {
  using Owner = Vof<M_>;
  using R = Reconst<Scal>;
  using TRM = Trackerm<M>;
  static constexpr size_t dim = M::dim;
  using Vect2 = GVect<Scal, 2>;
  using Sem = typename M::Sem;
  using MIdx = typename M::MIdx;

  Imp(Owner* owner, const FieldCell<Scal>& fcu, const FieldCell<Scal>& fccl,
      const MapCondFaceAdvection<Scal>& mfc, Par par)
      : owner_(owner)
      , par(par)
      , m(owner_->m)
      , layers(1)
      , mfc_(mfc)
      , fccl_(fccl)
      , fcim_(m, TRM::Pack(MIdx(0)))
      , fca_(m, GetNan<Scal>())
      , fcn_(m, GetNan<Vect>()) {
    fcu_.time_curr = fcu;

    UpdateBc(mfc_);

    for (auto c : m.AllCells()) {
      if (fcu[c] == 0) {
        fccl_[c] = kClNone;
      }
    }
  }
  void UpdateBc(const MapCondFaceAdvection<Scal>& mfc) {
    UVof<M>::GetAdvectionFaceCond(
        m, mfc, mfc_vf_, mfc_cl_, mfc_im_, mfc_n_, mfc_a_);
  }
  // reconstruct interface
  // uc: volume fraction [a]
  void Rec(const FieldCell<Scal>& uc) {
    DetectInterface(uc);
    // Compute fcn_ [s]
    UNormal<M>::CalcNormal(m, uc, fci_, par.dim, fcn_);
    auto h = m.GetCellSize();
    // Compute fca_ [s]
    for (auto c : m.SuCells()) {
      if (fci_[c]) {
        fca_[c] = R::GetLineA(fcn_[c], uc[c], h);
      } else {
        fca_[c] = GetNan<Scal>();
      }
    }
  }
  void DetectInterface(const FieldCell<Scal>& uc) {
    fci_.Reinit(m, false);
    // cell is 0<u<1
    for (auto c : m.AllCells()) {
      Scal u = uc[c];
      if (u > 0. && u < 1.) {
        fci_[c] = true;
      }
    }
    // cell is u=1 and neighbour is u=0
    for (auto c : m.SuCells()) {
      if (uc[c] == 1) {
        for (auto q : m.Nci(c)) {
          IdxCell cn = m.GetCell(c, q);
          if (uc[cn] == 0) {
            fci_[c] = true;
          }
        }
      }
    }
  }
  void StartStep() {
    auto sem = m.GetSem("start");
    if (sem("rotate")) {
      owner_->ClearIter();
      fcu_.time_prev = fcu_.time_curr;
      fcu_.iter_curr = fcu_.time_prev;
      UpdateBc(mfc_);
    }

    if (owner_->GetTime() == 0.) {
      if (sem("reconst")) {
        Rec(fcu_.time_curr);
      }
    }
  }
  // Makes advection sweep in one direction, updates uc [i]
  // uc: volume fraction [s]
  // d: direction
  // ffv: mixture flux [i]
  // fcn,fca: normal and plane constant [s]
  // mfc: face conditions, nullptr to keep boundary fluxes
  // type: 0: plain,
  //       1: Euler Explicit,
  //       2: Lagrange Explicit,
  //       3: Weymouth 2010
  // fcfm,fcfp: upwind mixture flux, required if type=2 [s]
  // fcuu: volume fraction for Weymouth div term
  // dt: time step
  // clipth: threshold for clipping, values outside [th,1-th] are clipped
  static void Sweep(
      FieldCell<Scal>& uc, size_t d, const FieldFace<Scal>& ffv,
      FieldCell<Scal>& fccl, FieldCell<Scal>& fcim, const FieldCell<Vect>& fcn,
      const FieldCell<Scal>& fca, const MapCondFace* mfc, int type,
      const FieldCell<Scal>* fcfm, const FieldCell<Scal>* fcfp,
      const FieldCell<Scal>* fcuu, Scal dt, Scal clipth, const M& m) {
    using Dir = typename M::Dir;
    using MIdx = typename M::MIdx;
    Dir md(d); // direction as Dir
    MIdx wd(md); // offset in direction d
    auto& bc = m.GetIndexCells();
    auto& bf = m.GetIndexFaces();
    MIdx gs = m.GetGlobalSize();
    auto h = m.GetCellSize();

    FieldFace<Scal> ffvu(m); // phase 2 flux

    if (!(type >= 0 && type <= 3)) {
      throw std::runtime_error(
          "Sweep(): unknown type '" + std::to_string(type) + "'");
    }

    // compute fluxes [i] and propagate color to downwind cells
    for (auto f : m.Faces()) {
      auto p = bf.GetMIdxDir(f);
      Dir df = p.second;

      if (df != md) {
        continue;
      }

      const Scal v = ffv[f]; // mixture flux
      IdxCell c = m.GetCell(f, v > 0. ? 0 : 1); // upwind cell
      if (uc[c] > 0 && uc[c] < 1) { // interfacial cell
        if (type == 0 || type == 1 || type == 3) {
          ffvu[f] = R::GetLineFlux(fcn[c], fca[c], h, v, dt, d);
        } else if (type == 2) { // Lagrange Explicit
          Scal vu = (v > 0. ? (*fcfm)[c] : (*fcfp)[c]);
          ffvu[f] = R::GetLineFluxStr(fcn[c], fca[c], h, v, vu, dt, d);
        }
      } else { // pure cell
        ffvu[f] = v * uc[c];
      }

      // propagate color to downwind cell if empty
      if (fccl[c] != kClNone) {
        IdxCell cd = m.GetCell(f, v > 0. ? 1 : 0); // downwind cell
        if (fccl[cd] == kClNone) {
          fccl[cd] = fccl[c];
          MIdx w = bc.GetMIdx(c);
          MIdx im = TRM::Unpack(fcim[c]);
          if (w[d] < 0) im[d] += 1;
          if (w[d] >= gs[d]) im[d] -= 1;
          fcim[cd] = TRM::Pack(im);
        }
      }
    }

    if (mfc) {
      // override flux in upwind boundaries
      FieldFace<Scal> ffu(m);
      InterpolateB(uc, *mfc, ffu, m);
      for (const auto& it : *mfc) {
        IdxFace f = it.first;
        Scal v = ffv[f];
        if ((it.second->GetNci() == 0) != (v > 0.)) {
          ffvu[f] = v * ffu[f];
        }
      }
    }

    // update volume fraction [i]
    for (auto c : m.Cells()) {
      auto w = bc.GetMIdx(c);
      const Scal lc = m.GetVolume(c);
      IdxFace fm = bf.GetIdx(w, md);
      IdxFace fp = bf.GetIdx(w + wd, md);
      // mixture cfl
      const Scal ds = (ffv[fp] - ffv[fm]) * dt / lc;
      // phase 2 cfl
      const Scal dl = (ffvu[fp] - ffvu[fm]) * dt / lc;
      auto& u = uc[c];
      if (type == 0) { // plain
        u += -dl;
      } else if (type == 1) { // Euler Implicit
        u = (u - dl) / (1. - ds);
      } else if (type == 2) { // Lagrange Explicit
        u += u * ds - dl;
      } else if (type == 3) { // Weymouth
        u += (*fcuu)[c] * ds - dl;
      }
      // clip
      if (!(u >= clipth)) {
        u = 0;
      } else if (!(u <= 1 - clipth)) {
        u = 1;
      }
      // clear color
      if (u == 0) {
        fccl[c] = kClNone;
        fcim[c] = TRM::Pack(MIdx(0));
      }
    }
  }
  void CommRec(
      Sem& sem, FieldCell<Scal>& uc, FieldCell<Scal>& fccl,
      FieldCell<Scal>& fcim) {
    if (sem("comm")) {
      m.Comm(&uc);
      m.Comm(&fccl);
      m.Comm(&fcim);
    }
    if (sem("bcreflect")) {
      BcApply(uc, mfc_vf_, m);
      BcApply(fccl, mfc_cl_, m);
      BcApply(fcim, mfc_im_, m);
    }
    if (sem("reconst")) {
      Rec(uc);
    }
  }
  void AdvAulisa(Sem& sem) {
    // directions, format: {dir EI, dir LE, ...}
    std::vector<size_t> dd;
    Scal vsc; // scaling factor for time step
    if (par.dim == 3) { // 3d
      if (count_ % 3 == 0) {
        dd = {0, 1, 1, 2, 2, 0};
      } else if (count_ % 3 == 1) {
        dd = {1, 2, 2, 0, 0, 1};
      } else {
        dd = {2, 0, 0, 1, 1, 2};
      }
      vsc = 0.5;
    } else { // 2d
      if (count_ % 2 == 0) {
        dd = {0, 1};
      } else {
        dd = {1, 0};
      }
      vsc = 1.0;
    }
    for (size_t id = 0; id < dd.size(); ++id) {
      size_t d = dd[id]; // direction as index
      if (sem("copyface")) {
        if (id % 2 == 1) { // copy fluxes for Lagrange Explicit step
          auto& ffv = *owner_->ffv_; // [f]ield [f]ace [v]olume flux
          fcfm_.Reinit(m);
          fcfp_.Reinit(m);
          for (auto c : m.Cells()) {
            fcfm_[c] = ffv[m.GetFace(c, 2 * d)];
            fcfp_[c] = ffv[m.GetFace(c, 2 * d + 1)];
          }
          m.Comm(&fcfm_);
          m.Comm(&fcfp_);
        }
      }
      auto& uc = fcu_.iter_curr;
      if (sem("sweep")) {
        Sweep(
            uc, d, *owner_->ffv_, fccl_, fcim_, fcn_, fca_, &mfc_vf_,
            id % 2 == 0 ? 1 : 2, &fcfm_, &fcfp_, nullptr,
            owner_->GetTimeStep() * vsc, par.clipth, m);
      }
      CommRec(sem, uc, fccl_, fcim_);
    }
  }
  static void BcMarchFill(FieldCell<Scal>& fcu, Scal fill, const M& m) {
    for (auto c : m.AllCells()) {
      auto x = m.GetCenter(c);
      if (!(Vect(0) <= x && x <= m.GetGlobalLength())) {
        fcu[c] = fill;
      }
    }
  }
  void AdvPlain(Sem& sem, int type) {
    std::vector<size_t> dd; // sweep directions
    if (par.dim == 3) { // 3d
      if (count_ % 3 == 0) {
        dd = {0, 1, 2};
      } else if (count_ % 3 == 1) {
        dd = {1, 2, 0};
      } else {
        dd = {2, 0, 1};
      }
    } else { // 2d
      if (count_ % 2 == 0) {
        dd = {0, 1};
      } else {
        dd = {1, 0};
      }
    }
    for (size_t id = 0; id < dd.size(); ++id) {
      auto& uc = fcu_.iter_curr;
      if (sem("sweep")) {
        Sweep(
            uc, dd[id], *owner_->ffv_, fccl_, fcim_, fcn_, fca_, &mfc_vf_, type,
            nullptr, nullptr, &fcuu_, owner_->GetTimeStep(), par.clipth, m);
      }
      CommRec(sem, uc, fccl_, fcim_);
    }
  }
  void Sharpen() {
    auto sem = m.GetSem("sharp");
    std::vector<size_t> dd; // sweep directions
    if (par.dim == 3) { // 3d
      if (count_ % 3 == 0) {
        dd = {0, 0, 1, 1, 2, 2};
      } else if (count_ % 3 == 1) {
        dd = {1, 1, 2, 2, 0, 0};
      } else {
        dd = {2, 2, 0, 0, 1, 1};
      }
    } else { // 2d
      if (count_ % 2 == 0) {
        dd = {0, 0, 1, 1};
      } else {
        dd = {1, 1, 0, 0};
      }
    }
    for (size_t id = 0; id < dd.size(); ++id) {
      size_t d = dd[id]; // direction as index
      auto& uc = fcu_.iter_curr;
      if (sem("sweep")) {
        const Scal sgn = (id % 2 == count_ / par.dim % 2 ? -1 : 1);
        FieldFace<Scal> ffv(m, m.GetCellSize().prod() * sgn * par.sharpen_cfl);
        // zero flux on boundaries
        for (const auto& it : mfc_vf_) {
          ffv[it.first] = 0;
        }
        Sweep(
            uc, d, ffv, fccl_, fcim_, fcn_, fca_, &mfc_vf_, 3, nullptr, nullptr,
            &fcuu_, 1., par.clipth, m);
      }
      CommRec(sem, uc, fccl_, fcim_);
    }
  }
  void MakeIteration() {
    auto sem = m.GetSem("iter");
    struct {
      FieldCell<Scal> fcclm; // previous color
    } * ctx(sem);
    auto& fcclm = ctx->fcclm;
    if (sem("init")) {
      auto& uc = fcu_.iter_curr;
      const Scal dt = owner_->GetTimeStep();
      auto& fcs = *owner_->fcs_;
      fcuu_.Reinit(m);
      for (auto c : m.Cells()) {
        uc[c] = fcu_.time_prev[c] + dt * fcs[c];
        fcuu_[c] = (uc[c] < 0.5 ? 0 : 1);
      }
    }

    using Scheme = typename Par::Scheme;
    switch (par.scheme) {
      case Scheme::plain:
        AdvPlain(sem, 0);
        break;
      case Scheme::aulisa:
        AdvAulisa(sem);
        break;
      case Scheme::weymouth:
        AdvPlain(sem, 3);
        break;
    }

    if (par.sharpen && sem.Nested("sharpen")) {
      Sharpen();
    }
    if (sem("bcc_clear")) {
      UVof<M>::BcClear(fcu_.iter_curr, mfc_, m);
    }
    if (sem("resetcolor")) {
      auto& fcu = fcu_.iter_curr;
      fcclm = fccl_;
      for (auto c : m.AllCells()) {
        fccl_[c] = (fcu[c] > 0.5 ? 0 : kClNone);
      }
      BcApply(fccl_, mfc_cl_, m);
    }
    if (sem.Nested()) {
      uvof_.Recolor(
          layers, &fcu_.iter_curr, &fccl_, &fcclm, par.clfixed, par.clfixed_x,
          par.coalth, mfc_cl_, par.verb, par.recolor_unionfind,
          par.recolor_reduce, par.recolor_grid, m);
    }
    if (sem("propagate")) {
      auto& u = fcu_.iter_curr;
      auto& cl = fccl_;
      for (auto f : m.Faces()) {
        for (size_t q : {0, 1}) {
          auto c = m.GetCell(f, q);
          auto cn = m.GetCell(f, 1 - q);
          if (cl[c] == kClNone && u[c] > 0 && u[cn] > u[c] &&
              cl[cn] != kClNone) {
            cl[c] = cl[cn];
          }
        }
      }
      m.Comm(&cl);
    }
    if (sem("stat")) {
      owner_->IncIter();
      ++count_;
    }
  }
  void FinishStep() {
    fcu_.time_curr = fcu_.iter_curr;
    owner_->IncTime();
  }
  void PostStep() {
    auto sem = m.GetSem("iter");
    if (sem("comm")) {
      // --> fcu [a]
      // --> fca [s], fcn [s]
      m.Comm(&fca_);
      m.Comm(&fcn_);
    }
    if (sem("reflect")) {
      // --> fca [a], fcn [a]
      BcApply(fca_, mfc_a_, m);
      BcApply(fcn_, mfc_n_, m);
      // --> reflected fca [a], fcn [a]
    }
  }
  void DumpInterface(std::string filename) {
    uvof_.DumpPoly(
        fcu_.time_curr, fcn_, fca_, fci_, filename, owner_->GetTime(),
        par.vtkbin, par.vtkmerge, m);
  }
  void DumpInterfaceMarch(std::string filename) {
    auto sem = m.GetSem("dump-interface-march");
    struct {
      FieldCell<Scal> fcut; // volume fraction
      FieldCell<Scal> fcclt; // color
    } * ctx(sem);
    auto& fcut = ctx->fcut;
    auto& fcclt = ctx->fcclt;
    if (sem("copy")) {
      fcut = fcu_.time_curr;
      fcclt = fccl_;
      if (par.bcc_reflectpoly) {
        BcReflectAll(fcut, mfc_vf_, m);
        BcReflectAll(fcclt, mfc_cl_, m);
      }
      if (par.dumppolymarch_fill >= 0) {
        BcMarchFill(fcut, par.dumppolymarch_fill, m);
      }
    }
    if (sem.Nested()) {
      uvof_.DumpPolyMarch(
          layers, &fcut, &fcclt, &fcn_, &fca_, &fci_, filename,
          owner_->GetTime(), par.vtkbin, par.vtkmerge, par.vtkiso,
          par.dumppolymarch_fill >= 0 ? &fcut : nullptr, m);
    }
  }

  Owner* owner_;
  Par par;
  M& m;
  GRange<size_t> layers;

  StepData<FieldCell<Scal>> fcu_;
  FieldCell<Scal> fcuu_; // volume fraction for Weymouth div term
  const MapCondFaceAdvection<Scal>& mfc_; // conditions on advection
  MapCondFace mfc_vf_; // conditions on vf
  MapCondFace mfc_cl_; // conditions on cl
  MapCondFace mfc_im_; // conditions on im
  MapCondFace mfc_n_; // conditions on n
  MapCondFace mfc_a_; // conditions on a

  FieldCell<Scal> fccl_; // color
  FieldCell<Scal> fcim_; // image
  FieldCell<Scal> fca_; // alpha (plane constant)
  FieldCell<Vect> fcn_; // n (normal to plane)
  FieldCell<bool> fci_; // interface mask (1: contains interface)
  size_t count_ = 0; // number of MakeIter() calls, used for splitting

  // tmp for MakeIteration, volume flux copied to cells
  FieldCell<Scal> fcfm_, fcfp_;
  UVof<M> uvof_;
};

template <class M_>
Vof<M_>::Vof(
    M& m, const FieldCell<Scal>& fcu, const FieldCell<Scal>& fccl,
    const MapCondFaceAdvection<Scal>& mfc, const FieldFace<Scal>* ffv,
    const FieldCell<Scal>* fcs, double t, double dt, Par par)
    : AdvectionSolver<M>(t, dt, m, ffv, fcs)
    , imp(new Imp(this, fcu, fccl, mfc, par)) {}

template <class M_>
Vof<M_>::~Vof() = default;

template <class M_>
auto Vof<M_>::GetPar() const -> const Par& {
  return imp->par;
}

template <class M_>
void Vof<M_>::SetPar(Par par) {
  imp->par = par;
}

template <class M_>
void Vof<M_>::StartStep() {
  imp->StartStep();
}

template <class M_>
void Vof<M_>::MakeIteration() {
  imp->MakeIteration();
}

template <class M_>
void Vof<M_>::FinishStep() {
  imp->FinishStep();
}

template <class M_>
auto Vof<M_>::GetField(Step l) const -> const FieldCell<Scal>& {
  return imp->fcu_.Get(l);
}

template <class M_>
auto Vof<M_>::GetColor() const -> const FieldCell<Scal>& {
  return imp->fccl_;
}

template <class M_>
auto Vof<M_>::GetImage(IdxCell c) const -> MIdx {
  return Trackerm<M>::Unpack(imp->fcim_[c]);
}

template <class M_>
auto Vof<M_>::GetMask() const -> const FieldCell<bool>& {
  return imp->fci_;
}

template <class M_>
auto Vof<M_>::GetAlpha() const -> const FieldCell<Scal>& {
  return imp->fca_;
}

template <class M_>
auto Vof<M_>::GetNormal() const -> const FieldCell<Vect>& {
  return imp->fcn_;
}

template <class M_>
void Vof<M_>::PostStep() {
  return imp->PostStep();
}

template <class M_>
void Vof<M_>::DumpInterface(std::string fn) const {
  return imp->DumpInterface(fn);
}

template <class M_>
void Vof<M_>::DumpInterfaceMarch(std::string fn) const {
  return imp->DumpInterfaceMarch(fn);
}
