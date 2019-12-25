#pragma once

#include <array>
#include <exception>
#include <fstream>
#include <limits>
#include <memory>
#include <set>

#include "approx.h"
#include "debug/isnan.h"
#include "geom/block.h"
#include "multi.h"
#include "normal.h"
#include "partstrmeshm.h"
#include "reconst.h"
#include "trackerm.h"
#include "util/vof.h"
#include "vofm.h"

namespace solver {

template <class M_>
struct Vofm<M_>::Imp {
  using Owner = Vofm<M_>;
  using R = Reconst<Scal>;
  using PS = PartStr<Scal>;
  using PSM = PartStrMeshM<M>;
  using TRM = Trackerm<M>;
  static constexpr Scal kClNone = -1;
  static constexpr size_t dim = M::dim;
  using Vect2 = GVect<Scal, 2>;
  using Sem = typename M::Sem;

  Imp(Owner* owner, const FieldCell<Scal>& fcu0, const FieldCell<Scal>& fccl0,
      const MapCondFaceAdvection<Scal>& mfc, std::shared_ptr<Par> par)
      : owner_(owner)
      , par(par)
      , m(owner_->m)
      , mfc_(mfc)
      , layers(0, par->layers) {
    fcu_.time_curr.resize(layers);
    fcuu_.resize(layers);
    fcn_.resize(layers);
    fca_.resize(layers);
    fci_.resize(layers);
    fck_.resize(layers);
    fccl_.resize(layers);
    fcim_.resize(layers);
    ffvu_.resize(layers);
    ffcl_.resize(layers);
    ffi_.resize(layers);

    fcn_.InitAll(FieldCell<Vect>(m, GetNan<Vect>()));
    fca_.InitAll(FieldCell<Scal>(m, GetNan<Scal>()));
    fck_.InitAll(FieldCell<Scal>(m, GetNan<Scal>()));
    fci_.InitAll(FieldCell<bool>(m, false));

    fcus_.time_curr = fcu0;

    fcu_.time_curr[0].Reinit(m, 0);
    fcu_.time_curr.InitAll(fcu_.time_curr[0]);

    fccl_.InitAll(FieldCell<Scal>(m, kClNone));
    fcim_.InitAll(FieldCell<Scal>(m, TRM::Pack(MIdx(0))));
    ffvu_.InitAll(FieldFace<Scal>(m, 0));
    ffcl_.InitAll(FieldFace<Scal>(m, kClNone));
    ffi_.InitAll(FieldFace<bool>(m, false));

    fcu_.time_curr[0] = fcu0;
    fccl_[0] = fccl0;
    fccls_ = fccl0;

    UpdateBc(mfc_);

    // particles
    auto ps = std::make_shared<typename PS::Par>();
    Update(ps.get());
    auto psm = std::make_shared<typename PSM::Par>();
    psm->ps = ps;
    psm->ns = par->part_ns;
    psm->tol = par->part_tol;
    psm->itermax = par->part_itermax;
    psm->verb = par->part_verb;
    psm->dim = par->dim;
    psm->dump_fr = par->part_dump_fr;
    psm->vtkbin = par->vtkbin;
    psm->vtkmerge = par->vtkmerge;
    psm_ = std::unique_ptr<PSM>(new PSM(m, psm, layers));
  }
  void UpdateBc(const MapCondFaceAdvection<Scal>& mfc) {
    UVof<M>::GetAdvectionFaceCond(
        m, mfc, mfc_vf_, mfc_cl_, mfc_im_, mfc_n_, mfc_a_);
  }
  void Update(typename PS::Par* p) const {
    Scal hc = m.GetCellSize().norminf(); // cell size

    p->leq = par->part_h;
    p->relax = par->part_relax;
    p->npmax = par->part_np;
    p->segcirc = par->part_segcirc;
    p->hc = hc;
    p->dn = par->part_dn;
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
        for (auto c : m.SuCells()) {
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
      for (auto c : m.SuCells()) {
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
        if (w > 0 && us >= par->avgnorm0) {
          na /= w;
          for (auto i : layers) {
            auto& n = fcn_[i][c];
            if (!IsNan(n)) {
              // average normal oriented to have acute angle with current
              auto nal = (n.dot(na) < 0 ? -na : na);
              Scal u0 = par->avgnorm0;
              Scal u1 = par->avgnorm1;
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

        auto h = m.GetCellSize();
        for (auto c : m.SuCells()) {
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
      for (auto c : m.AllCells()) {
        Scal u = fcu[c];
        if (u > 0. && u < 1.) {
          fci[c] = true;
        }
      }
      // cell is u=1 and neighbour is u=0
      for (auto c : m.SuCells()) {
        if (fcu[c] == 1) {
          for (auto q : m.Nci(c)) {
            bool b = false;
            for (auto j : layers) {
              IdxCell cn = m.GetCell(c, q);
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
  void Dump() {
    auto sem = m.GetSem("iter");
    struct {
      // tmp fields
      Multi<FieldCell<Scal>> fcut; // volume fraction
      Multi<FieldCell<Scal>> fcclt; // color
      FieldCell<Scal> fcust; // sum of  volume fraction
    } * ctx(sem);
    bool dm = par->dmp->Try(owner_->GetTime(), owner_->GetTimeStep());
    if (par->dumppoly && dm && sem.Nested()) {
      uvof_.DumpPoly(
          layers, fcu_.iter_curr, fccl_, fcn_, fca_, fci_,
          GetDumpName("s", ".vtk", par->dmp->GetN()),
          owner_->GetTime() + owner_->GetTimeStep(),
          par->vtkbin, par->vtkmerge, m);
    }
    if (par->dumppolymarch && dm) {
      auto& fcut = ctx->fcut; // tmp volume fraction
      auto& fcclt = ctx->fcclt; // tmp color;
      auto& fcust = ctx->fcust; // tmp total volume
      if (sem("copy")) {
        fcust = fcus_.iter_curr;
        fcclt = fccl_;
        fcut = fcu_.iter_curr;
        for (auto i : layers) {
          if (par->bcc_reflectpoly) {
            BcReflectAll(fcut[i], mfc_vf_, m);
            BcReflectAll(fcclt[i], mfc_cl_, m);
          }
        }
        if (par->dumppolymarch_fill >= 0) {
          BcMarchFill(fcust, par->dumppolymarch_fill, m);
        }
      }
      if (sem.Nested()) {
        uvof_.DumpPolyMarch(
            layers, fcut, fcclt, fcn_, fca_, fci_,
            GetDumpName("sm", ".vtk", par->dmp->GetN()),
            owner_->GetTime() + owner_->GetTimeStep(),
            par->vtkbin, par->vtkmerge, par->vtkiso,
            par->dumppolymarch_fill >= 0 ? &fcust : nullptr, m);
      }
    }
    if (par->dumppart && dm && sem.Nested("part-dump")) {
      psm_->DumpParticles(fca_, fcn_, par->dmp->GetN(), owner_->GetTime());
    }
    if (par->dumppartinter && dm && sem.Nested("partinter-dump")) {
      psm_->DumpPartInter(fca_, fcn_, par->dmp->GetN(), owner_->GetTime());
    }
  }
  void Sharpen(const Multi<FieldCell<Scal>*>& mfcu) {
    auto sem = m.GetSem("sharp");
    std::vector<size_t> dd; // sweep directions
    if (par->dim == 3) { // 3d
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
        const Scal sgn = (id % 2 == count_ / par->dim % 2 ? -1 : 1);
        FieldFace<Scal> ffv(m, m.GetCellSize().prod() * sgn * par->sharpen_cfl);
        // zero flux on boundaries
        for (const auto& it : mfc_) {
          IdxFace f = it.GetIdx();
          ffv[f] = 0;
        }
        Sweep(
            mfcu, d, layers, ffv, fccl_, fcim_, fcn_, fca_, mfc_vf_, 3, nullptr,
            nullptr, fcuu_, 1., par->clipth, m);
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
  void AdvPlain(Sem& sem, const Multi<FieldCell<Scal>*>& mfcu, int type) {
    std::vector<size_t> dd; // sweep directions
    if (par->dim == 3) { // 3d
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
            mfcu, d, layers, *owner_->ffv_, fccl_, fcim_, fcn_, fca_, mfc_vf_,
            type, nullptr, nullptr, fcuu_, owner_->GetTimeStep(), par->clipth,
            m);
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
  // type: 0: plain,
  //       1: Euler Explicit,
  //       2: Lagrange Explicit,
  //       3: Weymouth 2010
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
      const Multi<const FieldCell<Scal>*>& mfca, const MapCondFace& mfc,
      int type, const FieldCell<Scal>* fcfm, const FieldCell<Scal>* fcfp,
      const Multi<const FieldCell<Scal>*>& mfcuu, Scal dt, Scal clipth,
      const M& m) {
    using Dir = typename M::Dir;
    using MIdx = typename M::MIdx;
    Dir md(d); // direction as Dir
    MIdx wd(md); // offset in direction d
    auto& bc = m.GetIndexCells();
    auto& bf = m.GetIndexFaces();
    MIdx gs = m.GetGlobalSize();
    auto h = m.GetCellSize();

    Multi<FieldFace<Scal>> mffvu(layers); // phase 2 flux
    Multi<FieldFace<Scal>> mffcl(layers); // face color

    if (!(type >= 0 && type <= 3)) {
      throw std::runtime_error(
          "Sweep(): unknown type '" + std::to_string(type) + "'");
    }

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
      for (auto f : m.Faces()) {
        auto p = bf.GetMIdxDir(f);
        Dir df = p.second;

        if (df != md) {
          continue;
        }

        const Scal v = ffv[f]; // mixture flux
        IdxCell c = m.GetCell(f, v > 0. ? 0 : 1); // upwind cell
        if (fccl[c] != kClNone) {
          ffcl[f] = fccl[c];
          if (fcu[c] > 0 && fcu[c] < 1) {
            if (type == 0 || type == 1 || type == 3) {
              ffvu[f] = R::GetLineFlux(fcn[c], fca[c], h, v, dt, d);
            } else if (type == 2) {
              Scal vu = (v > 0. ? (*fcfm)[c] : (*fcfp)[c]);
              ffvu[f] = R::GetLineFluxStr(fcn[c], fca[c], h, v, vu, dt, d);
            }
          } else { // pure cell
            ffvu[f] = v * fcu[c];
          }

          // propagate to downwind cell if empty
          IdxCell cd = m.GetCell(f, v > 0. ? 1 : 0); // downwind cell
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

      FieldFace<Scal> ffu(m);
      InterpolateB(fcu, mfc, ffu, m);
      // override boundary flux
      for (const auto& it : mfc) {
        IdxFace f = it.GetIdx();
        if (ffcl[f] != kClNone) {
          Scal v = ffv[f];
          if ((it.GetValue()->GetNci() == 0) != (v > 0.)) {
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
      for (auto c : m.Cells()) {
        if (fccl[c] != kClNone) {
          auto w = bc.GetMIdx(c);
          const Scal lc = m.GetVolume(c);
          IdxFace fm = bf.GetIdx(w, md);
          IdxFace fp = bf.GetIdx(w + wd, md);
          // mixture cfl
          const Scal ds = (ffv[fp] - ffv[fm]) * dt / lc;
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
          const Scal dl = (vp - vm) * dt / lc;
          auto& u = fcu[c];
          if (type == 0) { // plain
            u += -dl;
          } else if (type == 1) { // Euler Implicit
            u = (u - dl) / (1. - ds);
          } else if (type == 2) { // Lagrange Explicit
            u += u * ds - dl;
          } else if (type == 3) { // Weymouth
            u += fcuu[c] * ds - dl;
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
        BcApply(*mfcu[i], mfc_vf_, m);
        BcApply(*mfccl[i], mfc_cl_, m);
        BcApply(*mfcim[i], mfc_im_, m);
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
    if (par->dim == 3) { // 3d
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
      if (sem("sweep")) {
        Sweep(
            mfcu, d, layers, *owner_->ffv_, fccl_, fcim_, fcn_, fca_, mfc_vf_,
            id % 2 == 0 ? 1 : 2, &fcfm_, &fcfp_, nullptr,
            owner_->GetTimeStep() * vsc, par->clipth, m);
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
        for (auto c : m.Cells()) {
          uc[c] = ucm[c] + dt * fcs[c];
          fcuu[c] = (uc[c] < 0.5 ? 0 : 1);
        }
      }
      ctx->fcclm = fccl_;
    }

    using Scheme = typename Par::Scheme;
    Multi<FieldCell<Scal>*> mfcu = fcu_.iter_curr;
    switch (par->scheme) {
      case Scheme::plain:
        AdvPlain(sem, mfcu, 0);
        break;
      case Scheme::aulisa:
        AdvAulisa(sem, mfcu);
        break;
      case Scheme::weymouth:
        AdvPlain(sem, mfcu, 3);
        break;
    }

    if (par->sharpen && sem.Nested("sharpen")) {
      Sharpen(mfcu);
    }
    if (sem("bcc_clear")) {
      if (par->cloverride) {
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
          layers, fcu_.iter_curr, fccl_, fccl_, par->clfixed,
          par->clfixed_x, par->coalth, mfc_cl_, par->verb,
          par->recolor_unionfind, par->recolor_reduce, par->recolor_grid, m);
    }

    if (sem("sum")) {
      auto l = Layers::iter_curr;
      auto& fcus = fcus_.Get(l);
      fcus.Reinit(m, 0);
      fccls_.Reinit(m, kClNone);
      for (auto i : layers) {
        auto& fccl = fccl_[i];
        auto& fcu = fcu_.Get(l)[i];
        for (auto c : m.AllCells()) {
          if (fccl[c] != kClNone) {
            fcus[c] += fcu[c];
            fccls_[c] = fccl[c];
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
        BcApply(fcn_[i], mfc_n_, m);
        BcApply(fca_[i], mfc_a_, m);
      }
      // --> reflected fca [a], fcn [a]
    }
    if (par->part && sem.Nested("part")) {
      psm_->Part(fca_, fcn_, fci_, fccl_);
    }
    if (sem.Nested("dump")) {
      Dump();
    }
  }

  Owner* owner_;
  std::shared_ptr<Par> par;
  M& m;

  LayersData<Multi<FieldCell<Scal>>> fcu_;
  Multi<FieldCell<Scal>> fcuu_;
  LayersData<FieldCell<Scal>> fcus_;
  FieldCell<Scal> fccls_;
  const MapCondFaceAdvection<Scal>& mfc_; // conditions on advection
  MapCondFace mfc_vf_; // conditions on vf
  MapCondFace mfc_cl_; // conditions on cl
  MapCondFace mfc_im_; // conditions on cl
  MapCondFace mfc_n_; // conditions on n
  MapCondFace mfc_a_; // conditions on a

  Multi<FieldCell<Scal>> fca_; // alpha (plane constant)
  Multi<FieldCell<Vect>> fcn_; // n (normal to plane)
  Multi<FieldCell<Scal>> fck_; // curvature from height functions
  Multi<FieldCell<bool>> fci_; // interface mask (1: contains interface)
  Multi<FieldCell<Scal>> fccl_; // color
  Multi<FieldCell<Scal>> fcim_; // image
  Multi<FieldFace<Scal>> ffvu_; // flux: volume flux * field
  Multi<FieldFace<Scal>> ffcl_; // flux color (from upwind cell)
  Multi<FieldFace<bool>>
      ffi_; // interface mask (1: upwind cell contains interface)
  size_t count_ = 0; // number of MakeIter() calls, used for splitting

  std::unique_ptr<PSM> psm_;
  // tmp for MakeIteration, volume flux copied to cells
  FieldCell<Scal> fcfm_, fcfp_;
  GRange<size_t> layers;
  UVof<M> uvof_;
};

template <class M>
constexpr typename M::Scal Vofm<M>::Imp::kClNone;

template <class M_>
Vofm<M_>::Vofm(
    M& m, const FieldCell<Scal>& fcu, const FieldCell<Scal>& fccl,
    const MapCondFaceAdvection<Scal>& mfc, const FieldFace<Scal>* ffv,
    const FieldCell<Scal>* fcs, double t, double dt, std::shared_ptr<Par> par)
    : AdvectionSolver<M>(t, dt, m, ffv, fcs)
    , imp(new Imp(this, fcu, fccl, mfc, par)) {}

template <class M_>
Vofm<M_>::~Vofm() = default;

template <class M_>
auto Vofm<M_>::GetPar() -> Par* {
  return imp->par.get();
}

template <class M_>
void Vofm<M_>::StartStep() {
  imp->StartStep();
}

template <class M_>
void Vofm<M_>::MakeIteration() {
  imp->MakeIteration();
}

template <class M_>
void Vofm<M_>::FinishStep() {
  imp->FinishStep();
}

template <class M_>
auto Vofm<M_>::GetField(Layers l) const -> const FieldCell<Scal>& {
  return imp->fcus_.Get(l);
}

template <class M_>
auto Vofm<M_>::GetField(Layers l, size_t i) const -> const FieldCell<Scal>& {
  return imp->fcu_.Get(l)[i];
}

template <class M_>
auto Vofm<M_>::GetFieldM() const -> Multi<const FieldCell<Scal>*> {
  return imp->fcu_.Get(Layers::time_curr);
}

template <class M_>
auto Vofm<M_>::GetAlpha() const -> Multi<const FieldCell<Scal>*> {
  return imp->fca_;
}

template <class M_>
auto Vofm<M_>::GetImage(size_t l, IdxCell c) const -> MIdx {
  return Trackerm<M>::Unpack(imp->fcim_[l][c]);
}

template <class M_>
auto Vofm<M_>::GetColor() const -> Multi<const FieldCell<Scal>*> {
  return imp->fccl_;
}

template <class M_>
auto Vofm<M_>::GetColorSum() const -> const FieldCell<Scal>& {
  return imp->fccls_;
}

template <class M_>
size_t Vofm<M_>::GetNumLayers() const {
  return imp->layers.size();
}

template <class M_>
auto Vofm<M_>::GetNormal() const -> Multi<const FieldCell<Vect>*> {
  return imp->fcn_;
}

template <class M_>
auto Vofm<M_>::GetCurv() const -> Multi<const FieldCell<Scal>*> {
  return imp->psm_->GetCurv();
}

template <class M_>
auto Vofm<M_>::GetCurvSum() const -> const FieldCell<Scal>& {
  return *imp->psm_->GetCurv()[0];
}

template <class M_>
void Vofm<M_>::PostStep() {
  return imp->PostStep();
}

} // namespace solver
