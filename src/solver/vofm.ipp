// Created by Petr Karnakov on 01.09.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <array>
#include <exception>
#include <fstream>
#include <limits>
#include <memory>
#include <set>

#include "approx.h"
#include "approx_eb.h"
#include "debug/isnan.h"
#include "geom/block.h"
#include "multi.h"
#include "normal.h"
#include "reconst.h"
#include "trackerm.h"
#include "util/convdiff.h"
#include "util/vof.h"
#include "vofm.h"

template <class EB_>
struct Vofm<EB_>::Imp {
  using Owner = Vofm<EB_>;
  using R = Reconst<Scal>;
  using TRM = Trackerm<M>;
  static constexpr size_t dim = M::dim;
  using Vect2 = generic::Vect<Scal, 2>;
  using Sem = typename M::Sem;
  using UEB = UEmbed<M>;

  Imp(Owner* owner, const EB& eb0, const GRange<size_t> layers0,
      const Multi<const FieldCell<Scal>*>& fcu0,
      const Multi<const FieldCell<Scal>*>& fccl0,
      const MapCondFaceAdvection<Scal>& mfc, Par par)
      : owner_(owner)
      , par(par)
      , m(owner_->m)
      , eb(eb0)
      , layers(layers0)
      , fcuu_(layers, m)
      , fccls_(m, kClNone)
      , fcn_(layers, m, GetNan<Vect>())
      , fca_(layers, m, GetNan<Scal>())
      , fci_(layers, m, false)
      , fccl_(fccl0)
      , fcim_(layers, m, TRM::Pack(MIdx(0)))
      , ffvu_(layers, m, 0)
      , ffcl_(layers, m, kClNone)
      , ffi_(layers, m, false)
      , mfc_(mfc) {
    fcu0.assert_size(layers);
    fccl0.assert_size(layers);
    fcu_.time_curr = fcu0;

    CalcSum(
        layers, fcu_.Get(Step::time_curr), fccl_, fcus_.Get(Step::time_curr),
        fccls_, m);

    UpdateBc(mfc_);
  }
  // Computes combined volume fraction and color from layers
  static void CalcSum(
      const GRange<size_t>& layers, const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<const FieldCell<Scal>*>& fccl, FieldCell<Scal>& fcus,
      FieldCell<Scal>& fccls, const M& m) {
    fcus.Reinit(m, 0);
    fccls.Reinit(m, kClNone);
    for (auto l : layers) {
      for (auto c : m.AllCells()) {
        const Scal cl = (*fccl[l])[c];
        if (cl != kClNone) {
          const Scal u = (*fcu[l])[c];
          fcus[c] += u;
          fccls[c] = cl;
        }
      }
    }
    for (auto c : m.AllCells()) {
      auto& u = fcus[c];
      if (!(u >= 0)) {
        u = 0;
      } else if (!(u <= 1.)) {
        u = 1;
      }
    }
  }
  void UpdateBc(const MapCondFaceAdvection<Scal>& mfc) {
    std::tie(me_vf_, me_cl_, me_im_, me_n_, me_a_) =
        UVof<M>::GetAdvectionBc(m, mfc);
    mfc_cl_ = GetCond<Scal>(me_cl_);
  }
  // reconstruct interface
  void Rec(const Multi<FieldCell<Scal>*>& uc) {
    auto sem = m.GetSem("rec");
    if (sem("detect")) {
      DetectInterface(uc);
    }
    if (sem("local")) {
      // Compute fcn_ [s]
      for (auto i : layers) {
        auto& fcn = fcn_[i];
        auto& fci = fci_[i];
        auto& fccl = fccl_[i];

        const int sw = 1; // stencil halfwidth
        using MIdx = typename M::MIdx;
        GBlock<IdxCell, dim> bo(MIdx(-sw), MIdx(sw * 2 + 1));
        for (auto c : eb.SuCells()) {
          if (fci[c]) {
            auto uu = GetStencil<M, 1>{}(layers, uc, fccl_, c, fccl[c], m);
            fcn[c] = UNormal<M>::GetNormalYoungs(uu);
            UNormal<M>::GetNormalHeight(uu, fcn[c]);
          } else {
            fcn[c] = GetNan<Vect>();
          }
        }
      }

      // Override with average fcn_ [s]
      for (auto c : eb.SuCells()) {
        Vect na(0); // sum of normals
        Scal w = 0; // number of normals
        Scal us = 0; // sum of volume fraction
        for (auto i : layers) {
          auto& n = fcn_[i][c];
          if (!IsNan(n)) {
            na += (n.dot(na) < 0 ? -n : n);
            w += 1;
            us += (*uc[i])[c];
          }
        }
        if (w > 0 && us >= par.avgnorm0) {
          na /= w;
          for (auto i : layers) {
            auto& n = fcn_[i][c];
            if (!IsNan(n)) {
              // average normal oriented to have acute angle with current
              auto nal = (n.dot(na) < 0 ? -na : na);
              Scal u0 = par.avgnorm0;
              Scal u1 = par.avgnorm1;
              if (u0 < u1) {
                Scal a = std::min(1., std::max(0., (us - u0) / (u1 - u0)));
                n = nal * a + n * (1 - a);
              } else {
                n = nal;
              }
            }
          }
        }
      }

      // Compute fca_ [s]
      for (auto i : layers) {
        auto& fcn = fcn_[i];
        auto& fci = fci_[i];
        auto& fcu = *uc[i];
        auto& fca = fca_[i];

        auto h = eb.GetCellSize();
        for (auto c : eb.SuCells()) {
          if (fci[c]) {
            fca[c] = R::GetLineA(fcn[c], fcu[c], h);
          } else {
            fca[c] = GetNan<Scal>();
          }
        }
      }
    }
  }
  void DetectInterface(const Multi<const FieldCell<Scal>*>& uc) {
    // cell is 0<u<1
    for (auto i : layers) {
      auto& fci = fci_[i];
      auto& fcu = *uc[i];
      fci.Reinit(m, false);
      for (auto c : eb.AllCells()) {
        Scal u = fcu[c];
        if (u > 0. && u < 1.) {
          fci[c] = true;
        }
      }
      // cell is u=1 and neighbour is u=0
      for (auto c : eb.SuCells()) {
        if (fcu[c] == 1) {
          for (auto q : eb.Nci(c)) {
            bool b = false;
            for (auto j : layers) {
              IdxCell cn = eb.GetCell(c, q);
              if (fccl_[j][cn] == fccl_[i][c]) {
                if ((*uc[j])[cn] == 0) {
                  fci[c] = true;
                }
                b = true;
                break;
              }
            }
            if (!b) {
              fci[c] = true;
            }
          }
        }
      }
    }
  }
  void StartStep() {
    auto sem = m.GetSem("start");
    if (owner_->GetTime() == 0.) {
      if (sem.Nested("reconst")) {
        Rec(fcu_.time_curr);
      }
    }
    if (sem("rotate")) {
      owner_->ClearIter();
      fcu_.time_prev = fcu_.time_curr;
      fcu_.iter_curr = fcu_.time_prev;
      fcus_.time_prev = fcus_.time_curr;
      fcus_.iter_curr = fcus_.time_prev;
      UpdateBc(mfc_);
    }
  }
  enum class SweepType {
    plain, // sum of fluxes
    EI, // Euler Implicit (aulisa2009)
    LE, // Lagrange Explicit (aulisa2009)
    weymouth, // sum of fluxes and divergence (weymouth2010)
  };
  void DumpInterface(std::string filename) {
    uvof_.DumpPoly(
        layers, fcu_.time_curr, fccl_, fcn_, fca_, fci_, filename,
        owner_->GetTime(), par.vtkbin, par.vtkmerge, m);
  }
  void DumpInterfaceMarch(std::string filename) {
    auto sem = m.GetSem("dump-interface-march");
    struct {
      Multi<FieldCell<Scal>> fcut; // volume fraction
      Multi<FieldCell<Scal>> fcclt; // color
      FieldCell<Scal> fcust; // sum of  volume fraction
    } * ctx(sem);
    auto& fcut = ctx->fcut;
    auto& fcclt = ctx->fcclt;
    auto& fcust = ctx->fcust;
    if (sem("copy")) {
      fcust = fcus_.time_curr;
      fcclt = fccl_;
      fcut = fcu_.time_curr;
      for (auto i : layers) {
        if (par.bcc_reflectpoly) {
          BcReflectAll(fcut[i], me_vf_, m);
          BcReflectAll(fcclt[i], me_cl_, m);
        }
      }
      if (par.dumppolymarch_fill >= 0) {
        BcMarchFill(fcust, par.dumppolymarch_fill, m);
      }
    }
    if (sem.Nested()) {
      uvof_.DumpPolyMarch(
          layers, fcut, fcclt, fcn_, fca_, fci_, filename, owner_->GetTime(),
          par.vtkbin, par.vtkmerge, par.vtkiso,
          par.dumppolymarch_fill >= 0 ? &fcust : nullptr, m);
    }
  }
  void Sharpen(const Multi<FieldCell<Scal>*>& mfcu) {
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
      if (sem("sweep")) {
        const Scal sgn = (id % 2 == count_ / par.dim % 2 ? -1 : 1);
        FieldFace<Scal> ffv(m, 0);
        for (auto f : eb.Faces()) {
          const IdxCell cm = m.GetCell(f, 0);
          const IdxCell cp = m.GetCell(f, 1);
          ffv[f] = std::min(eb.GetVolume(cm), eb.GetVolume(cp)) * sgn *
                   par.sharpen_cfl;
        }
        // zero flux on boundaries
        for (const auto& it : mfc_.GetMapFace()) {
          ffv[it.first] = 0;
        }
        Sweep(
            mfcu, d, layers, ffv, fccl_, fcim_, fcn_, fca_, me_vf_,
            SweepType::weymouth, nullptr, nullptr, fcuu_, 1., par.clipth, eb);
      }
      CommRec(sem, mfcu, fccl_, fcim_);
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
  void AdvPlain(Sem& sem, const Multi<FieldCell<Scal>*>& mfcu, SweepType type) {
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
      size_t d = dd[id]; // direction as index
      if (sem("sweep")) {
        Sweep(
            mfcu, d, layers, owner_->fev_->GetFieldFace(), fccl_, fcim_, fcn_,
            fca_, me_vf_, type, nullptr, nullptr, fcuu_, owner_->GetTimeStep(),
            par.clipth, eb);
      }
      CommRec(sem, mfcu, fccl_, fcim_);
    }
  }
  // Makes advection sweep in one direction, updates uc [i] and fccl [i]
  // uc: volume fraction [s]
  // d: direction
  // ffv: mixture flux [i]
  // mfccl: color [s]
  // mfcim: image [s]
  // mfcn,mfca: normal and plane constant [s]
  // mfc: face conditions
  // type: sweep type
  // fcfm,fcfp: upwind mixture flux, required if type=2 [s]
  // fcuu: volume fraction for Weymouth div term
  // dt: time step
  // clipth: threshold for clipping, values outside [th,1-th] are clipped
  static void Sweep(
      const Multi<FieldCell<Scal>*>& mfcu, size_t d,
      const GRange<size_t>& layers, const FieldFace<Scal>& ffv,
      const Multi<FieldCell<Scal>*>& mfccl,
      const Multi<FieldCell<Scal>*>& mfcim,
      const Multi<const FieldCell<Vect>*>& mfcn,
      const Multi<const FieldCell<Scal>*>& mfca,
      const MapEmbed<BCond<Scal>>& mebc, SweepType type,
      const FieldCell<Scal>* fcfm, const FieldCell<Scal>* fcfp,
      const Multi<const FieldCell<Scal>*>& mfcuu, Scal dt, Scal clipth,
      const EB& eb) {
    using Dir = typename M::Dir;
    using MIdx = typename M::MIdx;
    const Dir md(d); // direction as Dir
    const MIdx wd(md); // offset in direction d
    const auto& m = eb.GetMesh();
    const auto& bc = m.GetIndexCells();
    const auto& bf = m.GetIndexFaces();
    const MIdx gs = m.GetGlobalSize();
    const auto h = m.GetCellSize();

    Multi<FieldFace<Scal>> mffvu(layers); // phase 2 flux
    Multi<FieldFace<Scal>> mffcl(layers); // face color

    for (auto i : layers) {
      auto& fcu = *mfcu[i];
      auto& fcn = *mfcn[i];
      auto& fca = *mfca[i];
      auto& fccl = *mfccl[i];
      auto& ffvu = mffvu[i];
      auto& ffcl = mffcl[i];

      // compute fluxes [i] and propagate color to downwind cells
      ffvu.Reinit(m, 0);
      ffcl.Reinit(m, kClNone);
      for (auto f : eb.Faces()) {
        auto p = bf.GetMIdxDir(f);
        Dir df = p.second;

        if (df != md) {
          continue;
        }

        // flux through face (maybe cut)
        const Scal v = ffv[f];
        // flux through full face that would give the same velocity
        const Scal v0 = v / eb.GetAreaFraction(f);
        const IdxCell c = m.GetCell(f, v > 0. ? 0 : 1); // upwind cell
        if (fccl[c] != kClNone) {
          ffcl[f] = fccl[c];
          if (fcu[c] > 0 && fcu[c] < 1) {
            switch (type) {
              case SweepType::plain:
              case SweepType::EI:
              case SweepType::weymouth: {
                const Scal vu0 = R::GetLineFlux(fcn[c], fca[c], h, v0, dt, d);
                ffvu[f] = (v >= 0 ? std::min(vu0, v) : std::max(vu0, v));
                break;
              }
              case SweepType::LE: {
                const Scal vc = (v > 0. ? (*fcfm)[c] : (*fcfp)[c]);
                const Scal vc0 = vc / eb.GetAreaFraction(f);
                const Scal vu0 =
                    R::GetLineFluxStr(fcn[c], fca[c], h, v0, vc0, dt, d);
                ffvu[f] = (v >= 0 ? std::min(vu0, v) : std::max(vu0, v));
                break;
              }
            }
          } else { // pure cell
            ffvu[f] = v * fcu[c];
          }

          // propagate to downwind cell if empty
          IdxCell cd = eb.GetCell(f, v > 0. ? 1 : 0); // downwind cell
          bool found = false; // found same color downwind
          for (auto j : layers) {
            if ((*mfccl[j])[cd] == fccl[c]) {
              found = true;
              break;
            }
          }
          if (!found) {
            for (auto j : layers) {
              if ((*mfccl[j])[cd] == kClNone) {
                (*mfccl[j])[cd] = fccl[c];
                MIdx w = bc.GetMIdx(c);
                MIdx im = TRM::Unpack((*mfcim[i])[c]);
                if (w[d] < 0) im[d] += 1;
                if (w[d] >= gs[d]) im[d] -= 1;
                (*mfcim[j])[cd] = TRM::Pack(im);
                break;
              }
            }
          }
        }
      }

      // override boundary flux
      const FieldFace<Scal> ffu = UEB::Interpolate(fcu, mebc, m);
      for (const auto& p : mebc.GetMapFace()) {
        const IdxFace f = p.first;
        const auto& bc = p.second;
        if (ffcl[f] != kClNone) {
          const Scal v = ffv[f];
          if ((bc.nci == 0) != (v > 0.)) {
            ffvu[f] = v * ffu[f];
          }
        }
      }
    }

    // update volume fraction [i]
    for (auto i : layers) {
      auto& fcu = *mfcu[i];
      auto& fcuu = *mfcuu[i];
      auto& fccl = *mfccl[i];
      for (auto c : eb.Cells()) {
        if (fccl[c] != kClNone) {
          auto w = bc.GetMIdx(c);
          const Scal vol = m.GetVolume(c);
          IdxFace fm = bf.GetIdx(w, md);
          IdxFace fp = bf.GetIdx(w + wd, md);
          // mixture cfl
          const Scal ds = (ffv[fp] - ffv[fm]) * dt / vol;
          // phase 2 cfl
          Scal vm = 0;
          Scal vp = 0;
          for (auto j : layers) {
            if (mffcl[j][fm] == fccl[c]) {
              vm = mffvu[j][fm];
              break;
            }
          }
          for (auto j : layers) {
            if (mffcl[j][fp] == fccl[c]) {
              vp = mffvu[j][fp];
              break;
            }
          }
          const Scal dl = (vp - vm) * dt / vol;
          auto& u = fcu[c];
          switch (type) {
            case SweepType::plain: {
              u += -dl;
              break;
            }
            case SweepType::EI: {
              u = (u - dl) / (1. - ds);
              break;
            }
            case SweepType::LE: {
              u += u * ds - dl;
              break;
            }
            case SweepType::weymouth: {
              u += fcuu[c] * ds - dl;
              break;
            }
          }

          // clip
          if (!(u >= clipth)) { // u < clipth or nan
            u = 0;
          } else if (!(u <= 1 - clipth)) {
            u = 1;
          }
          // clear color
          if (u == 0) {
            fccl[c] = kClNone;
            (*mfcim[i])[c] = TRM::Pack(MIdx(0));
          }
        }
      }
    }
  }
  void CommRec(
      Sem& sem, const Multi<FieldCell<Scal>*>& mfcu,
      const Multi<FieldCell<Scal>*>& mfccl,
      const Multi<FieldCell<Scal>*>& mfcim) {
    if (layers.size() * 3 <= m.GetMaxComm()) { // cubismnc or local
      if (sem("comm")) {
        for (auto i : layers) {
          m.Comm(mfcu[i]);
          m.Comm(mfccl[i]);
          m.Comm(mfcim[i]);
        }
      }
    } else { // legacy cubism
      for (auto i : layers) {
        if (sem("comm")) {
          m.Comm(mfcu[i]);
          m.Comm(mfccl[i]);
          m.Comm(mfcim[i]);
        }
      }
    }
    if (sem("bcreflect")) {
      for (auto i : layers) {
        BcApply(*mfcu[i], me_vf_, m);
        BcApply(*mfccl[i], me_cl_, m);
        BcApply(*mfcim[i], me_im_, m);
      }
    }
    if (sem.Nested("reconst")) {
      Rec(mfcu);
    }
  }
  void AdvAulisa(Sem& sem, const Multi<FieldCell<Scal>*>& mfcu) {
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
          auto& ffv =
              owner_->fev_->GetFieldFace(); // [f]ield [f]ace [v]olume flux
          fcfm_.Reinit(m);
          fcfp_.Reinit(m);
          for (auto c : eb.Cells()) {
            fcfm_[c] = ffv[eb.GetFace(c, 2 * d)];
            fcfp_[c] = ffv[eb.GetFace(c, 2 * d + 1)];
          }
          m.Comm(&fcfm_);
          m.Comm(&fcfp_);
        }
      }
      if (sem("sweep")) {
        Sweep(
            mfcu, d, layers, owner_->fev_->GetFieldFace(), fccl_, fcim_, fcn_,
            fca_, me_vf_, id % 2 == 0 ? SweepType::EI : SweepType::LE, &fcfm_,
            &fcfp_, nullptr, owner_->GetTimeStep() * vsc, par.clipth, eb);
      }
      CommRec(sem, mfcu, fccl_, fcim_);
    }
  }
  void MakeIteration() {
    auto sem = m.GetSem("iter");
    struct {
      Multi<FieldCell<Scal>> fcclm; // previous color
    } * ctx(sem);
    if (sem("init")) {
      for (auto l : layers) {
        auto& uc = fcu_.iter_curr[l];
        auto& ucm = fcu_.time_prev[l];
        auto& fcuu = fcuu_[l];
        fcuu.Reinit(m);
        const Scal dt = owner_->GetTimeStep();
        auto& fcs = *owner_->fcs_;
        for (auto c : eb.Cells()) {
          uc[c] = ucm[c] + dt * fcs[c];
          fcuu[c] = (uc[c] < 0.5 ? 0 : 1);
        }
      }
      ctx->fcclm = fccl_;
    }

    using Scheme = typename Par::Scheme;
    Multi<FieldCell<Scal>*> mfcu = fcu_.iter_curr;
    switch (par.scheme) {
      case Scheme::plain:
        AdvPlain(sem, mfcu, SweepType::plain);
        break;
      case Scheme::aulisa:
        AdvAulisa(sem, mfcu);
        break;
      case Scheme::weymouth:
        AdvPlain(sem, mfcu, SweepType::weymouth);
        break;
    }

    if (par.sharpen && sem.Nested("sharpen")) {
      Sharpen(mfcu);
    }
    if (sem("bcc_clear")) {
      if (par.cloverride) {
        for (auto i : layers) {
          UVof<M>::BcClearOverrideColor(
              fcu_.iter_curr[i], fccl_[i], 0., mfc_, m);
        }
      } else {
        for (auto i : layers) {
          UVof<M>::BcClear(fcu_.iter_curr[i], mfc_, m);
        }
      }
    }
    if (sem.Nested()) {
      uvof_.Recolor(
          layers, fcu_.iter_curr, fccl_, fccl_, par.clfixed, par.clfixed_x,
          par.coalth, mfc_cl_, par.verb, par.recolor_unionfind,
          par.recolor_reduce, par.recolor_grid, m);
    }

    if (sem("sum")) {
      CalcSum(
          layers, fcu_.Get(Step::iter_curr), fccl_, fcus_.Get(Step::iter_curr),
          fccls_, m);
    }
    if (sem("stat")) {
      owner_->IncIter();
      ++count_;
    }
  }
  void FinishStep() {
    fcu_.time_curr = fcu_.iter_curr;
    fcus_.time_curr = fcus_.iter_curr;
    owner_->IncTime();
  }
  void PostStep() {
    auto sem = m.GetSem("iter");
    // --> fcu [a], fca [s], fcn [s]
    if (layers.size() * 4 <= m.GetMaxComm()) { // cubismnc or local
      if (sem("comm")) {
        for (auto i : layers) {
          m.Comm(&fcn_[i]);
          m.Comm(&fca_[i]);
        }
      }
    } else { // legacy cubism
      for (auto i : layers) {
        if (sem("comm")) {
          m.Comm(&fcn_[i]);
          m.Comm(&fca_[i]);
        }
      }
    }
    if (sem("reflect")) {
      // --> fca [a], fcn [a]
      for (auto i : layers) {
        BcApply(fcn_[i], me_n_, m);
        BcApply(fca_[i], me_a_, m);
      }
      // --> reflected fca [a], fcn [a]
    }
  }

  Owner* owner_;
  Par par;
  M& m;
  const EB& eb;

  GRange<size_t> layers;
  StepData<Multi<FieldCell<Scal>>> fcu_;
  Multi<FieldCell<Scal>> fcuu_; // volume fraction for Weymouth div term
  StepData<FieldCell<Scal>> fcus_;
  FieldCell<Scal> fccls_;

  Multi<FieldCell<Vect>> fcn_; // n (normal to plane)
  Multi<FieldCell<Scal>> fca_; // alpha (plane constant)
  Multi<FieldCell<bool>> fci_; // interface mask (1: contains interface)
  Multi<FieldCell<Scal>> fccl_; // color
  Multi<FieldCell<Scal>> fcim_; // image
  Multi<FieldFace<Scal>> ffvu_; // flux: volume flux * field
  Multi<FieldFace<Scal>> ffcl_; // flux color (from upwind cell)
  Multi<FieldFace<bool>> ffi_; // interface mask
                               // (1: upwind cell contains interface)
  size_t count_ = 0; // number of MakeIter() calls, used for splitting

  // boundary conditions
  const MapCondFaceAdvection<Scal>& mfc_; // conditions on advection
  MapEmbed<BCond<Scal>> me_vf_; // volume fraction
  MapEmbed<BCond<Scal>> me_cl_; // color
  MapEmbed<BCond<Scal>> me_im_; // image
  MapEmbed<BCond<Vect>> me_n_; // normal
  MapEmbed<BCond<Scal>> me_a_; // plane constant
  MapCondFace mfc_cl_; // color

  // tmp for MakeIteration, volume flux copied to cells
  FieldCell<Scal> fcfm_, fcfp_;
  UVof<M> uvof_;
};

template <class EB_>
constexpr typename EB_::M::Scal Vofm<EB_>::kClNone;

template <class EB_>
Vofm<EB_>::Vofm(
    M& m, const EB& eb, const FieldCell<Scal>& fcu0,
    const FieldCell<Scal>& fccl0, const MapCondFaceAdvection<Scal>& mfc,
    const FieldEmbed<Scal>* fev, const FieldCell<Scal>* fcs, double t,
    double dt, Par par)
    : AdvectionSolver<M>(t, dt, m, fev, fcs) {
  const GRange<size_t> layers(par.layers);
  Multi<FieldCell<Scal>> fcu(layers, m, 0);
  Multi<FieldCell<Scal>> fccl(layers, m, kClNone);
  fcu[0] = fcu0;
  fccl[0] = fccl0;
  imp.reset(new Imp(this, eb, layers, fcu, fccl, mfc, par));
}

template <class EB_>
Vofm<EB_>::Vofm(
    M& m, const EB& eb, const Multi<const FieldCell<Scal>*>& fcu0,
    const Multi<const FieldCell<Scal>*>& fccl0,
    const MapCondFaceAdvection<Scal>& mfc, const FieldEmbed<Scal>* fev,
    const FieldCell<Scal>* fcs, double t, double dt, Par par)
    : AdvectionSolver<M>(t, dt, m, fev, fcs) {
  const GRange<size_t> layers(par.layers);
  imp.reset(new Imp(this, eb, layers, fcu0, fccl0, mfc, par));
}

template <class EB_>
Vofm<EB_>::~Vofm() = default;

template <class EB_>
auto Vofm<EB_>::GetEmbed() const -> const EB& {
  return imp->eb;
}

template <class EB_>
auto Vofm<EB_>::GetPar() const -> const Par& {
  return imp->par;
}

template <class EB_>
void Vofm<EB_>::SetPar(Par par) {
  imp->par = par;
}

template <class EB_>
void Vofm<EB_>::StartStep() {
  imp->StartStep();
}

template <class EB_>
void Vofm<EB_>::MakeIteration() {
  imp->MakeIteration();
}

template <class EB_>
void Vofm<EB_>::FinishStep() {
  imp->FinishStep();
}

template <class EB_>
auto Vofm<EB_>::GetField(Step l) const -> const FieldCell<Scal>& {
  return imp->fcus_.Get(l);
}

template <class EB_>
auto Vofm<EB_>::GetField(Step l, size_t i) const -> const FieldCell<Scal>& {
  return imp->fcu_.Get(l)[i];
}

template <class EB_>
auto Vofm<EB_>::GetFieldM() const -> Multi<const FieldCell<Scal>*> {
  return imp->fcu_.Get(Step::time_curr);
}

template <class EB_>
auto Vofm<EB_>::GetAlpha() const -> Multi<const FieldCell<Scal>*> {
  return imp->fca_;
}

template <class EB_>
auto Vofm<EB_>::GetImage(size_t l, IdxCell c) const -> MIdx {
  return Trackerm<M>::Unpack(imp->fcim_[l][c]);
}

template <class EB_>
auto Vofm<EB_>::GetMask() const -> Multi<const FieldCell<bool>*> {
  return imp->fci_;
}

template <class EB_>
auto Vofm<EB_>::GetColor() const -> Multi<const FieldCell<Scal>*> {
  return imp->fccl_;
}

template <class EB_>
auto Vofm<EB_>::GetColorSum() const -> const FieldCell<Scal>& {
  return imp->fccls_;
}

template <class EB_>
size_t Vofm<EB_>::GetNumLayers() const {
  return imp->layers.size();
}

template <class EB_>
auto Vofm<EB_>::GetNormal() const -> Multi<const FieldCell<Vect>*> {
  return imp->fcn_;
}

template <class EB_>
void Vofm<EB_>::PostStep() {
  return imp->PostStep();
}

template <class EB_>
void Vofm<EB_>::DumpInterface(std::string fn) const {
  return imp->DumpInterface(fn);
}

template <class EB_>
void Vofm<EB_>::DumpInterfaceMarch(std::string fn) const {
  return imp->DumpInterfaceMarch(fn);
}
