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
#include "partstr.h"
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

  Imp(Owner* owner, const FieldCell<Scal>& fcu0, const FieldCell<Scal>& fccl0,
      const MapFace<std::shared_ptr<CondFace>>& mfc, 
      std::shared_ptr<Par> par)
      : owner_(owner), par(par), m(owner_->m)
      , mfc_(mfc), layers(0, 4)
  {
    fcu_.resize(layers.size());
    fcn_.resize(layers.size());
    fca_.resize(layers.size());
    fci_.resize(layers.size());
    fck_.resize(layers.size());
    fccl_.resize(layers.size());
    fcclt_.resize(layers.size());
    ffvu_.resize(layers.size());
    ffcl_.resize(layers.size());
    ffi_.resize(layers.size());


    fcn_.InitAll(FieldCell<Vect>(m, Vect(0)));
    fca_.InitAll(FieldCell<Scal>(m, 0));
    fck_.InitAll(FieldCell<Scal>(m, 0));
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
    psm_ = std::unique_ptr<PSM>(new PSM(m, psm));
    psm->vtkbin = par->vtkbin;
    psm->vtkmerge = par->vtkmerge;
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
        auto& fccl = fccl_[i];
        if (par->bcc_reflect) {
          BcReflect(fcu, mfc_, par->bcc_fill, m);
        }

        const int sw = 1; // stencil halfwidth
        using MIdx = typename M::MIdx;
        GBlock<IdxCell, dim> bo(MIdx(-sw), MIdx(sw * 2 + 1));
        for (auto c : m.Cells()) {
          if (fci[c]) {
            auto uu = GetStencil<M, 1>{}(layers, uc, fccl_, c, fccl[c], m);
            fcn[c] = UNormal<M>::GetNormalYoungs(uu);
          }
        }
      }
      for (auto i : layers) {
        if (par->bcc_reflect) {
          auto& fcn = fcn_[i];
          BcReflect(fcn, mfc_, Vect(0), m);
        }
      }
      for (auto i : layers) {
        auto& fcn = fcn_[i];
        auto& fca = fca_[i];
        auto& fci = fci_[i];
        auto& fcu = *uc[i];
        for (auto c : m.SuCells()) {
          if (fci[c]) {
            fca[c] = R::GetLineA(fcn[c], fcu[c], m.GetCellSize());
          }
        }
        if (par->bcc_reflect) {
          BcReflect(fca, mfc_, Scal(0), m);
        }
      }
    }
    for (auto i : layers) {
      if (sem("comm")) {
        auto& fcn = fcn_[i];
        auto& fca = fca_[i];
        m.Comm(&fca);
        m.Comm(&fcn);
      }
    }
  }
  void DetectInterface(const Multi<const FieldCell<Scal>*>& uc) {
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
  void MakeIteration() {
    auto sem = m.GetSem("iter");
    if (sem("init")) {
      for (auto i : layers) {
        auto& fcu = fcu_[i];
        auto& uc = fcu.iter_curr;
        const Scal dt = owner_->GetTimeStep();
        auto& fcs = *owner_->fcs_;
        for (auto c : m.Cells()) {
          uc[c] =
              fcu.time_prev[c] +  // previous time step
              dt * fcs[c]; // source
        }
      }
    }

    using Dir = typename M::Dir;
    using MIdx = typename M::MIdx;
    // directions, format: {dir LE, dir EI, ...}
    std::vector<size_t> dd; 
    Scal vsc; // scaling factor for ffv, used for splitting
    if (par->dim == 3) { // 3d
      if (count_ % 3 == 0) {
        dd = {0, 1, 1, 2, 2, 0};
      } else if (count_ % 3 == 1) {
        dd = {1, 2, 2, 0, 0, 1};
      } else {
        dd = {2, 0, 0, 1, 1, 2};
      }
      vsc = 0.5;
    } else {
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
        for (auto i : layers) {
          ffvu_[i].Reinit(m, 0);
          ffcl_[i].Reinit(m, kClNone);
          ffi_[i].Reinit(m, false);
        }
      }
      for (auto i : layers) {
        if (sem("flux")) {
          Dir md(d); // direction as Dir
          MIdx wd(md); // offset in direction d
          auto& bf = m.GetIndexFaces();
          auto h = m.GetCellSize();
          const Scal dt = owner_->GetTimeStep();
          auto& ffv = *owner_->ffv_; // [f]ield [f]ace [v]olume flux

          // Returns phase 2 flux
          auto F = [this,vsc,id,h,dt,d,&ffv](
              IdxFace f, IdxCell cu, Scal u, Vect n, Scal a, bool in) {
            const Scal v = ffv[f] * vsc; // mixture flux
            if (in) {
              if (id % 2 == 0) { // Euler Implicit
                return R::GetLineFlux(n, a, h, v, dt, d);
              } else {           // Lagrange Explicit
                // upwind mixture flux
                Scal vu = (v > 0. ? fcfm_[cu] : fcfp_[cu]) * vsc;
                return R::GetLineFluxStr(n, a, h, v, vu, dt, d);
              }
            }
            return v * u;
          };

          auto& fcu = fcu_[i].iter_curr;
          auto& fcn = fcn_[i];
          auto& fca = fca_[i];
          auto& fci = fci_[i];
          auto& fccl = fccl_[i];
          auto& ffvu = ffvu_[i];
          auto& ffcl = ffcl_[i];
          auto& ffi = ffi_[i];

          for (auto f : m.Faces()) {
            if (bf.GetDir(f) != md) {
              continue;
            }

            // upwind cell
            IdxCell cu = m.GetNeighbourCell(f, ffv[f] > 0. ? 0 : 1);
            if (fccl[cu] != kClNone) {
              ffvu[f] = F(f, cu, fcu[cu], fcn[cu], fca[cu], fci[cu]);
              ffi[f] = fci[cu];
              ffcl[f] = fccl[cu];
            }
          }

          FieldFace<Scal> ffu(m);
          // interpolate field value to boundaries
          InterpolateB(fcu, mfc_, ffu, m);
          // override boundary upwind flux
          for (const auto& it : mfc_) {
            IdxFace f = it.GetIdx();
            CondFace* cb = it.GetValue().get(); 
            Scal v = ffv[f];
            if ((cb->GetNci() == 0) != (v > 0.)) {
              ffvu[f] = v * ffu[f];
              ffi[f] = true;
            }
          }
        }
      }
      if (sem("cell")) {
        for (auto i : layers) {
          Dir md(d); // direction as Dir
          MIdx wd(md); // offset in direction d
          const Scal dt = owner_->GetTimeStep();
          auto& bc = m.GetIndexCells();
          auto& bf = m.GetIndexFaces();
          auto& ffv = *owner_->ffv_; // mixture flux
          auto& fcu = fcu_[i].iter_curr;
          auto& fccl = fccl_[i];

          for (auto c : m.Cells()) {
            auto w = bc.GetMIdx(c);
            const Scal lc = m.GetVolume(c);
            Scal cl = fccl[c];
            if (cl == kClNone) {
              continue;
            }
            // faces
            IdxFace fm = bf.GetIdx(w, md);
            IdxFace fp = bf.GetIdx(w + wd, md);
            // mixture fluxes
            const Scal vm = ffv[fm] * vsc;
            const Scal vp = ffv[fp] * vsc;
            // mixture cfl
            const Scal sm = vm * dt / lc;
            const Scal sp = vp * dt / lc;
            const Scal ds = sp - sm;
            // phase 2 fluxes
            Scal qm = 0;
            Scal qp = 0;
            // interface
            bool im = false;
            bool ip = false;
            for (auto j : layers) {
              if (ffcl_[j][fm] == cl) {
                qm = ffvu_[j][fm];
                im = ffi_[j][fm];
                break;
              }
            }
            for (auto j : layers) {
              if (ffcl_[j][fp] == cl) {
                qp = ffvu_[j][fp];
                ip = ffi_[j][fp];
                break;
              }
            }
            if (!im && !ip) {
              //continue; // XXX
            }
            // phase 2 cfl
            const Scal lm = qm * dt / lc;
            const Scal lp = qp * dt / lc;
            const Scal dl = lp - lm;
            auto& u = fcu[c];
            if (id % 2 == 0) { // Euler Implicit
              u = (u - dl) / (1. - ds);
            } else { // Lagrange Explicit
              u = u * (1. + ds) - dl;
            }
            // clip
            const Scal th = par->clipth;
            if (u < th) {
              u = 0.;
            } else if (u > 1. - th) {
              u = 1.;
            } else if (IsNan(u)) {
              u = 0.;
            }
            // update color
            if (u == 0) {
              fccl[c] = kClNone;
            }
          }
        }

        for (auto i : layers) {
          Dir md(d); // direction as Dir
          MIdx wd(md); // offset in direction d
          const Scal dt = owner_->GetTimeStep();
          auto& bc = m.GetIndexCells();
          auto& bf = m.GetIndexFaces();
          auto& ffv = *owner_->ffv_; // mixture flux
          auto& ffcl = ffcl_[i];
          for (auto f : m.Faces()) {
            auto p = bf.GetMIdxDir(f);
            MIdx w = p.first;
            Dir df = p.second;
            if (df != md) {
              continue;
            }
            Scal cl = ffcl[f];
            if (cl == kClNone) {
              continue;
            }

            const Scal v = ffv[f] * vsc;

            IdxFace fm, fp;
            IdxCell c;
            if (v > 0) {
              fm = f;
              fp = bf.GetIdx(w + wd, md);
              c = bc.GetIdx(w);
            } else {
              fm = bf.GetIdx(w - wd, md);
              fp = f;
              c = bc.GetIdx(w - wd);
            }

            const Scal lc = m.GetVolume(c);

            // mixture fluxes
            const Scal vm = ffv[fm] * vsc;
            const Scal vp = ffv[fp] * vsc;
            // mixture cfl
            const Scal sm = vm * dt / lc;
            const Scal sp = vp * dt / lc;
            const Scal ds = sp - sm;
            // phase 2 fluxes
            Scal qm = 0;
            Scal qp = 0;
            // interface
            bool im = false;
            bool ip = false;
            // qm,im from layer with same color
            for (auto j : layers) {
              if (ffcl_[j][fm] == cl) {
                qm = ffvu_[j][fm];
                im = ffi_[j][fm];
                break;
              }
            }
            // qp,ip from layer with same color
            for (auto j : layers) {
              if (ffcl_[j][fp] == cl) {
                qp = ffvu_[j][fp];
                ip = ffi_[j][fp];
                break;
              }
            }
            if (!im && !ip) {
              //continue; // XXX
            }
            // phase 2 cfl
            const Scal lm = qm * dt / lc;
            const Scal lp = qp * dt / lc;
            const Scal dl = lp - lm;

            // find cell with same color
            const size_t jnone = -1;
            size_t j = jnone;
            for (auto jj : layers) {
              if (fccl_[jj][c] == cl) {
                j = jj;
                break;
              }
            }
            if (j == jnone) { // if not found, find first empty
              for (auto jj : layers) {
                if (fccl_[jj][c] == kClNone) {
                  j = jj;
                  break;
                }
              }
            }
            if (j == jnone) {
              //std::cerr << "warn: no layer i=" << i << " w=" << w << "\n";
              j = 0;
            }

            if (fccl_[j][c] != kClNone) {
              continue;
            }

            auto& u = fcu_[j].iter_curr[c];
            if (id % 2 == 0) { // Euler Implicit
              u = (u - dl) / (1. - ds);
            } else { // Lagrange Explicit
              u = u * (1. + ds) - dl;
            }

            // clip
            const Scal th = par->clipth;
            if (u < th) {
              u = 0.;
            } else if (u > 1. - th) {
              u = 1.;
            } else if (IsNan(u)) {
              u = 0.;
            }
            // update color
            auto& clc = fccl_[j][c];
            if (u == 0) {
              clc = kClNone;
            } else {
              clc = cl;
            }
          }
        }
        for (auto i : layers) {
          m.Comm(&fcu_[i].iter_curr);
          m.Comm(&fccl_[i]);
        }
      }
      if (sem.Nested("reconst")) {
        Rec(GetLayer(fcu_, Layers::iter_curr));
      }
    }
    if (sem.Nested("recolor")) {
      Recolor();
    }
    if (sem("stat")) {
      owner_->IncIter();
      ++count_;
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
    // Curvature from gradient of volume fraction
    if (par->curvgrad && sem("curv")) {
      for (auto i : layers) {
        auto ffu = Interpolate(fcu_[i].iter_curr, mfc_, m); // [s]
        auto fcg = Gradient(ffu, m); // [s]
        auto ffg = Interpolate(fcg, mfvz_, m); // [i]

        auto& fck = fck_[i];
        fck.Reinit(m, GetNan<Scal>()); // curvature [i]

        auto& fci = fci_[i];
        for (auto c : m.Cells()) {
          if (!fci[c]) {
            continue;
          }
          Scal s = 0.;
          for (auto q : m.Nci(c)) {
            IdxFace f = m.GetNeighbourFace(c, q);
            auto& g = ffg[f];
            auto n = g / g.norm();  // inner normal
            s += -n.dot(m.GetOutwardSurface(c, q));
          }
          fck[c] = s / m.GetVolume(c);
        }
      }
    }
    if (par->part && sem.Nested("part")) {
      psm_->Part(GetLayer(fcu_, Layers::iter_curr),
                 fca_, fcn_, fci_, fccl_, mfc_);
    }
    if (sem.Nested("dump")) {
      Dump();
    }
  }


  Owner* owner_;
  std::shared_ptr<Par> par;
  M& m;

  Multi<LayersData<FieldCell<Scal>>> fcu_;
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
