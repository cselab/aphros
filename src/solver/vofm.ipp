#pragma once

#include <exception>
#include <fstream>
#include <array>
#include <memory>
#include <limits>
#include <set>

#include "vofm.h"
#include "geom/block.h"
#include "reconst.h"
#include "normal.h"
#include "debug/isnan.h"
#include "partstrmeshm.h"
#include "multi.h"
#include "util/vof.h"

namespace solver {

template <class T>
Multi<FieldCell<T>*> GetLayer(Multi<LayersData<FieldCell<T>>>& u, Layers l) {
  Multi<FieldCell<T>*> r(u.size());
  for (size_t i = 0; i < u.size(); ++i) {
    r[i] = &u[i].Get(l);
  }
  return r;
}


template <class M_>
struct Vofm<M_>::Imp {
  using Owner = Vofm<M_>;
  using R = Reconst<Scal>;
  using PS = PartStr<Scal>;
  using PSM = PartStrMeshM<M>;
  static constexpr Scal kClNone = -1;
  static constexpr size_t dim = M::dim;
  using Vect2 = GVect<Scal, 2>;
  using Sem = typename M::Sem;

  Imp(Owner* owner, const FieldCell<Scal>& fcu0, const FieldCell<Scal>& fccl0,
      const MapFace<std::shared_ptr<CondFace>>& mfc, 
      std::shared_ptr<Par> par)
      : owner_(owner), par(par), m(owner_->m)
      , mfc_(mfc), layers(0, 4)
  {
    fcu_.resize(layers.size());
    fcuu_.resize(layers.size());
    fcn_.resize(layers.size());
    fca_.resize(layers.size());
    fci_.resize(layers.size());
    fck_.resize(layers.size());
    fccl_.resize(layers.size());
    fcclt_.resize(layers.size());
    ffvu_.resize(layers.size());
    ffcl_.resize(layers.size());
    ffi_.resize(layers.size());


    fcn_.InitAll(FieldCell<Vect>(m, GetNan<Vect>()));
    fca_.InitAll(FieldCell<Scal>(m, GetNan<Scal>()));
    fck_.InitAll(FieldCell<Scal>(m, GetNan<Scal>()));
    fci_.InitAll(FieldCell<bool>(m, false));

    fcus_.time_curr = fcu0;

    fcu_[0].time_curr.Reinit(m, 0);
    fcu_.InitAll(fcu_[0]);

    fccl_.InitAll(FieldCell<Scal>(m, kClNone));
    ffvu_.InitAll(FieldFace<Scal>(m, 0));
    ffcl_.InitAll(FieldFace<Scal>(m, kClNone));
    ffi_.InitAll(FieldFace<bool>(m, false));


    fcu_[0].time_curr = fcu0;
    fccl_[0] = fccl0;
    fccls_ = fccl0;

    for (auto it : mfc_) {
      IdxFace f = it.GetIdx();
      mfvz_[f] = std::make_shared<
          CondFaceGradFixed<Vect>>(Vect(0), it.GetValue()->GetNci());
    }

    // particles
    auto ps = std::make_shared<typename PS::Par>();
    Update(ps.get());
    auto psm = std::make_shared<typename PSM::Par>();
    psm->ps = ps;
    psm->intth = par->part_intth;
    psm->ns = par->part_ns;
    psm->tol = par->part_tol;
    psm->itermax = par->part_itermax;
    psm->verb = par->part_verb;
    psm->dim = par->dim;
    psm->bcc_reflect = par->bcc_reflect;
    psm->dump_fr = par->part_dump_fr;
    psm->maxr = par->part_maxr;
    psm->vtkbin = par->vtkbin;
    psm->vtkmerge = par->vtkmerge;
    psm_ = std::unique_ptr<PSM>(new PSM(m, psm, layers));
  }
  void Update(typename PS::Par* p) const {
    Scal hc = m.GetCellSize().norminf(); // cell size

    p->leq = par->part_h;
    p->relax = par->part_relax;
    p->constr = par->part_constr;
    p->npmax = par->part_np;
    p->segcirc = par->part_segcirc;
    p->hc = hc;

    p->kstr = par->part_kstr;
    p->kbend = par->part_kbend;
    p->kattr = par->part_kattr;
    p->bendmean = par->part_bendmean;
    p->ka = par->part_kattr;
    p->kt = par->part_kbend;
    p->kx = par->part_kstr;
    p->tmax = par->part_tmax;
    p->dtmax = par->part_dtmax;
    p->anglim = par->part_anglim;
    p->dn = par->part_dn;

    {
      using FT = typename PS::FT;
      switch (par->part_attrforce) { 
        case Par::AF::line:
          p->forcetype = FT::line;
          break;
        case Par::AF::center:
          p->forcetype = FT::center;
          break;
        case Par::AF::volume:
          p->forcetype = FT::volume;
          break;
        default:
          throw std::runtime_error("Update(): Unknown part_attrforce");
      }
    }
  }
  // reconstruct interface
  void Rec(const Multi<FieldCell<Scal>*>& uc) {
    auto sem = m.GetSem("rec");
    if (sem("detect")) {
      DetectInterface(uc);
    }
    if (sem("local")) {
      for (auto i : layers) {
        auto& fcn = fcn_[i];
        auto& fci = fci_[i];
        auto& fcu = *uc[i];
        auto& fca = fca_[i];
        auto& fccl = fccl_[i];

        const int sw = 1; // stencil halfwidth
        using MIdx = typename M::MIdx;
        GBlock<IdxCell, dim> bo(MIdx(-sw), MIdx(sw * 2 + 1));
        auto h = m.GetCellSize();
        // Compute fcn_, fca_ [s]
        for (auto c : m.SuCells()) {
          if (fci[c]) {
            auto uu = GetStencil<M, 1>{}(layers, uc, fccl_, c, fccl[c], m);
            fcn[c] = UNormal<M>::GetNormalYoungs(uu);
            UNormal<M>::GetNormalHeight(uu, fcn[c]);
            fca[c] = R::GetLineA(fcn[c], fcu[c], h);
          } else {
            fcn[c] = GetNan<Vect>();
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
              IdxCell cn = m.GetNeighbourCell(c, q);
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
        Rec(GetLayer(fcu_, Layers::time_curr));
      }
    }
    if (sem("rotate")) {
      owner_->ClearIter();
      for (auto& u : fcu_.data()) {
        u.time_prev = u.time_curr;
        u.iter_curr = u.time_prev;
      }
      fcus_.time_prev = fcus_.time_curr;
      fcus_.iter_curr = fcus_.time_prev;
    }
  }
  void Dump() {
    auto sem = m.GetSem("iter");
    bool dm = par->dmp->Try(owner_->GetTime(),
                            owner_->GetTimeStep());
    if (par->dumppoly && dm && sem.Nested()) {
      uvof_.DumpPoly(
          layers, GetLayer(fcu_, Layers::iter_curr), fccl_, fcn_, fca_, fci_,
          GetDumpName("s", ".vtk", par->dmp->GetN()),
          owner_->GetTime() + owner_->GetTimeStep(),
          par->poly_intth, par->vtkbin, par->vtkmerge, m);
    }
    if (par->dumppolymarch && dm && sem.Nested()) {
      uvof_.DumpPolyMarch(
          layers, GetLayer(fcu_, Layers::iter_curr), fccl_, fcn_, fca_, fci_,
          GetDumpName("sm", ".vtk", par->dmp->GetN()),
          owner_->GetTime() + owner_->GetTimeStep(),
          par->poly_intth, par->vtkbin, par->vtkmerge, par->vtkiso, m);
    }
    if (par->dumppart && dm && sem.Nested("part-dump")) {
      psm_->DumpParticles(fca_, fcn_, par->dmp->GetN(), owner_->GetTime());
    }
    if (par->dumppartinter && dm && sem.Nested("partinter-dump")) {
      psm_->DumpPartInter(fca_, fcn_, par->dmp->GetN(), owner_->GetTime());
    }
  }
  void Recolor() {
    auto sem = m.GetSem("recolor");
    if (sem("init")) {
      fcclt_.InitAll(FieldCell<Scal>(m, kClNone));
      // initial unique color
      Scal q = m.GetId() * m.GetInBlockCells().size() * layers.size();
      for (auto i : layers) {
        for (auto c : m.Cells()) {
          if (fccl_[i][c] != kClNone) {
            fcclt_[i][c] = (q += 1);
          }
        }
        m.Comm(&fcclt_[i]);
      }
    }
    sem.LoopBegin();
    if (sem("min")) {
      const int sw = 1; // stencil halfwidth
      using MIdx = typename M::MIdx;
      auto& bc = m.GetIndexCells();
      GBlock<IdxCell, dim> bo(MIdx(-sw), MIdx(sw * 2 + 1));
      size_t tries = 0;
      size_t cells = 0;
      while (true) {
        bool chg = false;
        for (auto i : layers) {
          for (auto c : m.Cells()) {
            if (fccl_[i][c] != kClNone) {
              MIdx w = bc.GetMIdx(c);
              for (MIdx wo : bo) {
                IdxCell cn = bc.GetIdx(w + wo);
                for (auto j : layers) {
                  if (fccl_[j][cn] == fccl_[i][c]) {
                    if (fcclt_[j][cn] < fcclt_[i][c]) {
                      chg = true;
                      ++cells;
                      fcclt_[i][c] = fcclt_[j][cn];
                    }
                  }
                }
              }
            }
          }
        }
        if (!chg) {
          break;
        }
        ++tries;
      }
      for (auto i : layers) {
        m.Comm(&fcclt_[i]);
      }
      recolor_tries_ = tries;
      m.Reduce(&recolor_tries_, "max");
    }
    if (sem("check")) {
      if (par->verb && m.IsRoot()) {
        std::cerr << "recolor:"
          << " max tries: " << recolor_tries_ 
          << std::endl;
      }
      if (!recolor_tries_) {
        sem.LoopBreak();
      }
    }
    sem.LoopEnd();
    if (sem("free")) {
      for (auto i : layers) {
        for (auto c : m.AllCells()) {
          fccl_[i][c] = fcclt_[i][c];
        }
        fcclt_[i].Free();
      }
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
        Sweep(mfcu, d, layers, ffv, fccl_, fcn_, fca_, mfc_,
            3, nullptr, nullptr, fcuu_, 1., par->clipth, m);
        for (auto i : layers) {
          m.Comm(mfcu[i]);
          m.Comm(&fccl_[i]);
        }
      }
      if (par->bcc_reflect && sem("reflect")) {
        for (auto i : layers) {
          BcReflect(*mfcu[i], mfc_, par->bcc_fill, m);
        }
      }
      if (sem.Nested("reconst")) {
        Rec(mfcu);
      }
    }
  }
  // set volume fraction to 0 or 1 near wall
  static void BcClear(FieldCell<Scal>& uc,
                      const MapFace<std::shared_ptr<CondFace>>& mfc, 
                      const M& m) {
    for (const auto& it : mfc) {
      CondFace* cb = it.GetValue().get();
      if (dynamic_cast<CondFaceReflect*>(cb)) {
        IdxFace f = it.GetIdx();
        CondFace* cb = it.GetValue().get(); 
        size_t nci = cb->GetNci();
        IdxCell c = m.GetNeighbourCell(f, nci);
        uc[c] = (uc[c] > 0.5 ? 1. : 0.);
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
        Sweep(mfcu, d, layers, *owner_->ffv_, fccl_, fcn_, fca_, mfc_,
              type, nullptr, nullptr, fcuu_,
              owner_->GetTimeStep(), par->clipth, m);
        for (auto i : layers) {
          m.Comm(mfcu[i]);
          m.Comm(&fccl_[i]);
        }
      }
      if (par->bcc_reflect && sem("reflect")) {
        for (auto i : layers) {
          BcReflect(*mfcu[i], mfc_, par->bcc_fill, m);
        }
      }
      if (sem.Nested("reconst")) {
        Rec(mfcu);
      }
    }
  }
  // Makes advection sweep in one direction, updates uc [i] and fccl [s]
  // uc: volume fraction [s]
  // d: direction
  // ffv: mixture flux [i]
  // fcn,fca: normal and plane constant [s]
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
      const GRange<size_t>& layers,
      const FieldFace<Scal>& ffv,
      const Multi<FieldCell<Scal>*>& mfccl,
      const Multi<const FieldCell<Vect>*>& mfcn,
      const Multi<const FieldCell<Scal>*>& mfca,
      const MapFace<std::shared_ptr<CondFace>>& mfc, int type,
      const FieldCell<Scal>* fcfm, const FieldCell<Scal>* fcfp,
      const Multi<const FieldCell<Scal>*>& mfcuu,
      Scal dt, Scal clipth, const M& m) {
    using Dir = typename M::Dir;
    using MIdx = typename M::MIdx;
    Dir md(d); // direction as Dir
    MIdx wd(md); // offset in direction d
    auto& bc = m.GetIndexCells();
    auto& bf = m.GetIndexFaces();
    auto h = m.GetCellSize();

    Multi<FieldFace<Scal>> mffvu(layers.size()); // phase 2 flux
    Multi<FieldFace<Scal>> mffcl(layers.size()); // face color

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
        IdxCell c = m.GetNeighbourCell(f, v > 0. ? 0 : 1); // upwind cell
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
          IdxCell cd = m.GetNeighbourCell(f, v > 0. ? 1 : 0); // downwind cell
          bool found = false;
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
            if (mffcl[j][fm] == fccl[c]) { vm = mffvu[j][fm]; break; }
          }
          for (auto j : layers) {
            if (mffcl[j][fp] == fccl[c]) { vp = mffvu[j][fp]; break; }
          }
          const Scal dl = (vp - vm) * dt / lc;
          auto& u = fcu[c];
          if (type == 0) {        // plain
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
          }
        }
      }
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
            fcfm_[c] = ffv[m.GetNeighbourFace(c, 2 * d)];
            fcfp_[c] = ffv[m.GetNeighbourFace(c, 2 * d + 1)];
          }
          m.Comm(&fcfm_);
          m.Comm(&fcfp_);
        }
      }
      if (sem("sweep")) {
        Sweep(mfcu, d, layers, *owner_->ffv_, fccl_, fcn_, fca_, mfc_,
              id % 2 == 0 ? 1 : 2, &fcfm_, &fcfp_, nullptr,
              owner_->GetTimeStep() * vsc, par->clipth, m);
        for (auto i : layers) {
          m.Comm(mfcu[i]);
          m.Comm(&fccl_[i]);
        }
      }
      if (par->bcc_reflect && sem("reflect")) {
        for (auto i : layers) {
          BcReflect(*mfcu[i], mfc_, par->bcc_fill, m);
        }
      }
      if (sem.Nested("reconst")) {
        Rec(mfcu);
      }
    }
  }
  void MakeIteration() {
    auto sem = m.GetSem("iter");
    if (sem("init")) {
      for (auto i : layers) {
        auto& fcu = fcu_[i];
        auto& uc = fcu.iter_curr;
        auto& fcuu = fcuu_[i];
        fcuu.Reinit(m);
        const Scal dt = owner_->GetTimeStep();
        auto& fcs = *owner_->fcs_;
        for (auto c : m.Cells()) {
          uc[c] = fcu.time_prev[c] + dt * fcs[c];
          fcuu[c] = (uc[c] < 0.5 ? 0 : 1);
        }
      }
    }

    using Scheme = typename Par::Scheme;
    auto mfcu = GetLayer(fcu_, Layers::iter_curr);
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
    if (par->bcc_clear && sem("clear")) {
      for (auto i : layers) {
        BcClear(fcu_[i].iter_curr, mfc_, m);
      }
    }
    if (sem.Nested("recolor")) {
      Recolor();
    }

    if (sem("sum")) {
      auto l = Layers::iter_curr;
      auto& fcus = fcus_.Get(l);
      fcus.Reinit(m, 0);
      fccls_.Reinit(m, kClNone);
      for (auto i : layers) {
        auto& fccl = fccl_[i];
        auto& fcu = fcu_[i].Get(l);
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
    for (auto& u : fcu_.data()) {
      u.time_curr = u.iter_curr;
    }
    fcus_.time_curr = fcus_.iter_curr;
    owner_->IncTime();
  }
  void PostStep() {
    auto sem = m.GetSem("iter");
    for (auto i : layers) {
      if (sem("comm")) {
        // --> fcu [a]
        // --> fca [s], fcn [s]
        m.Comm(&fcn_[i]);
        m.Comm(&fca_[i]);
      }
    }
    if (par->bcc_reflect && sem("reflect")) {
      // --> fca [a], fcn [a]
      for (auto i : layers) {
        BcReflect(fcn_[i], mfc_, Vect(0), m);
        BcReflect(fca_[i], mfc_, Scal(0), m);
      }
      // --> reflected fca [a], fcn [a]
    }
    if (par->part && sem.Nested("part")) {
      psm_->Part(GetLayer(fcu_, Layers::iter_curr),
                 fca_, fcn_, fci_, fccl_, nullptr, mfc_);
    }
    if (sem.Nested("dump")) {
      Dump();
    }
  }


  Owner* owner_;
  std::shared_ptr<Par> par;
  M& m;

  Multi<LayersData<FieldCell<Scal>>> fcu_;
  Multi<FieldCell<Scal>> fcuu_;
  LayersData<FieldCell<Scal>> fcus_;
  FieldCell<Scal> fccls_;
  MapFace<std::shared_ptr<CondFace>> mfc_;
  MapFace<std::shared_ptr<CondFace>> mfvz_; // zero-derivative bc for Vect

  Multi<FieldCell<Scal>> fca_; // alpha (plane constant)
  Multi<FieldCell<Vect>> fcn_; // n (normal to plane)
  Multi<FieldCell<Scal>> fck_; // curvature from height functions
  Multi<FieldCell<bool>> fci_; // interface mask (1: contains interface)
  Multi<FieldCell<Scal>> fccl_;  // color
  Multi<FieldCell<Scal>> fcclt_;  // tmp color
  Multi<FieldFace<Scal>> ffvu_;  // flux: volume flux * field
  Multi<FieldFace<Scal>> ffcl_;  // flux color (from upwind cell)
  Multi<FieldFace<bool>> ffi_;   // interface mask (1: upwind cell contains interface)
  size_t count_ = 0; // number of MakeIter() calls, used for splitting

  std::unique_ptr<PSM> psm_;
  // tmp for MakeIteration, volume flux copied to cells
  FieldCell<Scal> fcfm_, fcfp_;
  GRange<size_t> layers;
  Scal recolor_tries_;
  UVof<M> uvof_;
};

template <class M>
constexpr typename M::Scal Vofm<M>::Imp::kClNone;

template <class M_>
Vofm<M_>::Vofm(
    M& m, const FieldCell<Scal>& fcu, const FieldCell<Scal>& fccl,
    const MapFace<std::shared_ptr<CondFace>>& mfc,
    const FieldFace<Scal>* ffv, const FieldCell<Scal>* fcs,
    double t, double dt, std::shared_ptr<Par> par)
    : AdvectionSolver<M>(t, dt, m, ffv, fcs)
    , imp(new Imp(this, fcu, fccl, mfc, par))
{}

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
  return imp->fcu_[i].Get(l);
}

template <class M_>
auto Vofm<M_>::GetField(size_t i) const -> const FieldCell<Scal>& {
  return imp->fcu_[i].Get(Layers::time_curr);
}

template <class M_>
auto Vofm<M_>::GetAlpha(size_t i) const -> const FieldCell<Scal>& {
  (void) i;
  return imp->fca_[i];
}

template <class M_>
auto Vofm<M_>::GetColor(size_t i) const -> const FieldCell<Scal>& {
  return imp->fccl_[i];
}

template <class M_>
auto Vofm<M_>::GetColor() const -> const FieldCell<Scal>& {
  return imp->fccls_;
}

template <class M_>
size_t Vofm<M_>::GetNumLayers() const {
  return imp->layers.size();
}

template <class M_>
auto Vofm<M_>::GetNormal(size_t i) const -> const FieldCell<Vect>& {
  return imp->fcn_[i];
}

template <class M_>
auto Vofm<M_>::GetCurv() const -> const FieldCell<Scal>& {
  return imp->psm_->GetCurv(0);
}

template <class M_>
auto Vofm<M_>::GetCurv(size_t i) const -> const FieldCell<Scal>& {
  (void) i;
  return imp->psm_->GetCurv(i);
}

template <class M_>
void Vofm<M_>::PostStep() {
  return imp->PostStep();
}

} // namespace solver
