#pragma once

#include <exception>
#include <fstream>
#include <array>
#include <memory>
#include <limits>

#include "vof.h"
#include "geom/block.h"
#include "reconst.h"
#include "normal.h"
#include "debug/isnan.h"
#include "partstrmeshm.h"
#include "util/vof.h"

namespace solver {

template <class M_>
struct Vof<M_>::Imp {
  using Owner = Vof<M_>;
  using R = Reconst<Scal>;
  using PS = PartStr<Scal>;
  using PSM = PartStrMeshM<M>;
  static constexpr size_t dim = M::dim;
  using Vect2 = GVect<Scal, 2>;
  using Sem = typename M::Sem;

  Imp(Owner* owner, const FieldCell<Scal>& fcu,
      const MapFace<std::shared_ptr<CondFace>>& mfc, 
      std::shared_ptr<Par> par)
      : owner_(owner), par(par), m(owner_->m)
      , mfc_(mfc), fca_(m, GetNan<Scal>()), fcn_(m, GetNan<Vect>())
      , fck_(m, 0), fch_(m, Vect(0))
  {
    fcu_.time_curr = fcu;
    for (auto it : mfc_) {
      IdxFace f = it.GetIdx();
      mfvz_[f] = std::make_shared<
          CondFaceGradFixed<Vect>>(Vect(0), it.GetValue()->GetNci());
    }

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
    psm_ = std::unique_ptr<PSM>(new PSM(m, psm));
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

  // fcu: volume fraction [a]
  // fcud2: volume fraction difference double (xpp-xmm, ypp-ymm, zpp-zmm) [i]
  void CalcDiff2(const FieldCell<Scal>& fcu, FieldCell<Vect>& fcud2) {
    fcud2.Reinit(m);
    for (auto c : m.Cells()) {
      for (size_t d = 0; d < dim; ++d) {
        auto cm = m.GetNeighbourCell(c, 2 * d);
        auto cmm = m.GetNeighbourCell(cm, 2 * d);
        auto cp = m.GetNeighbourCell(c, 2 * d + 1);
        auto cpp = m.GetNeighbourCell(cp, 2 * d + 1);
        fcud2[c][d] = fcu[cpp] - fcu[cmm];
      }
    }
  }
  // fcu: volume fraction [a]
  // fcud2: volume fraction difference double (xpp-xmm, ...) [a]
  // Output:
  // fcud4: volume fraction difference quad (xp4-xm4, ...) [i]
  void CalcDiff4(const FieldCell<Scal>& fcu, 
      const FieldCell<Vect>& fcud2, FieldCell<Vect>& fcud4) {
    fcud4.Reinit(m);
    for (auto c : m.Cells()) {
      for (size_t d = 0; d < dim; ++d) {
        auto cm = m.GetNeighbourCell(c, 2 * d);
        auto cmm = m.GetNeighbourCell(cm, 2 * d);
        auto cp = m.GetNeighbourCell(c, 2 * d + 1);
        auto cpp = m.GetNeighbourCell(cp, 2 * d + 1);
        Scal um4 = fcu[c] - fcud2[cmm][d];
        Scal up4 = fcud2[cpp][d] + fcu[c];
        fcud4[c][d] = up4 - um4;
      }
    }
  }
  // fcu: volume fraction [a]
  // fcud4: volume fraction difference double (xp4-xm4, ...) [a]
  // Output:
  // fcud6: volume fraction difference quad (xp6-xm6, ...) [i]
  void CalcDiff6(const FieldCell<Scal>& fcu, 
      const FieldCell<Vect>& fcud4, FieldCell<Vect>& fcud6) {
    fcud6.Reinit(m);
    for (auto c : m.Cells()) {
      for (size_t d = 0; d < dim; ++d) {
        auto cm = m.GetNeighbourCell(c, 2 * d);
        auto cmm = m.GetNeighbourCell(cm, 2 * d);
        auto cp = m.GetNeighbourCell(c, 2 * d + 1);
        auto cpp = m.GetNeighbourCell(cp, 2 * d + 1);
        Scal um6 = fcu[cpp] - fcud4[cmm][d];
        Scal up6 = fcud4[cpp][d] + fcu[cmm];
        fcud6[c][d] = up6 - um6;
      }
    }
  }
  // reconstruct interface
  // uc: volume fraction [a]
  void Rec(const FieldCell<Scal>& uc) {
    auto sem = m.GetSem("rec");
    if (sem("local")) {
      DetectInterface(uc);
      // Compute normal fcn_ [s]
      UNormal<M>::CalcNormal(m, uc, fci_, par->dim, fcn_);
      auto h = m.GetCellSize();
      // Reconstruct interface fca_ [s]
      for (auto c : m.SuCells()) {
        if (fci_[c]) {
          fca_[c] = R::GetLineA(fcn_[c], uc[c], h);
        } else {
          fca_[c] = GetNan<Scal>();
        }
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
          IdxCell cn = m.GetNeighbourCell(c, q);
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
    }

    if (owner_->GetTime() == 0.) {
      if (sem.Nested("reconst")) {
        Rec(fcu_.time_curr);
      }
    }
  }
  void Dump() {
    auto sem = m.GetSem("iter");
    bool dm = par->dmp->Try(owner_->GetTime(),
                            owner_->GetTimeStep());
    if (par->dumppoly && dm && sem.Nested("dumppoly")) {
      uvof_.DumpPoly(
          fcu_.iter_curr, fcn_, fca_, fci_,
          GetDumpName("s", ".vtk", par->dmp->GetN()),
          owner_->GetTime() + owner_->GetTimeStep(),
          par->poly_intth, par->vtkbin, par->vtkmerge, m);
    }
    if (par->dumppart && dm && sem.Nested("part-dump")) {
      psm_->DumpParticles(&fca_, &fcn_, par->dmp->GetN(), owner_->GetTime());
    }
    if (par->dumppartinter && dm && sem.Nested("partinter-dump")) {
      psm_->DumpPartInter(&fca_, &fcn_, par->dmp->GetN(), owner_->GetTime());
    }
  }
  // Makes advection sweep in one direction, updates uc [i]
  // uc: volume fraction [s]
  // d: direction
  // ffv: mixture flux [i]
  // fcn,fca: normal and plane constant [s]
  // mfc: face conditions
  // type: 0: plain,
  //       1: Euler Explicit,
  //       2: Lagrange Explicit,
  //       3: Weymouth 2010
  // vsc: scale for mixture flux
  // fcfm,fcfp: upwind mixture flux, required if le=1 [s]
  // dt: time step
  void Sweep(FieldCell<Scal>& uc, size_t d,
             const FieldFace<Scal>& ffv,
             const FieldCell<Vect>& fcn, const FieldCell<Scal>& fca,
             const MapFace<std::shared_ptr<CondFace>>& mfc,
             int type, Scal vsc,
             const FieldCell<Scal>& fcfm, const FieldCell<Scal>& fcfp,
             Scal dt, const M& m) {
    using Dir = typename M::Dir;
    using MIdx = typename M::MIdx;
    Dir md(d); // direction as Dir
    MIdx wd(md); // offset in direction d
    auto& bc = m.GetIndexCells();
    auto& bf = m.GetIndexFaces();
    auto h = m.GetCellSize();

    FieldFace<Scal> ffvu(m); // phase 2 flux

    if (!(type >= 0 && type <= 3)) {
      throw std::runtime_error(
          "Sweep(): unknown type '" + std::to_string(type) + "'");
    }

    // compute fluxes [i]
    for (auto f : m.Faces()) {
      auto p = bf.GetMIdxDir(f);
      Dir df = p.second;

      if (df != md) {
        continue;
      }

      const Scal v = ffv[f] * vsc; // mixture flux
      IdxCell cu = m.GetNeighbourCell(f, v > 0. ? 0 : 1); // upwind cell
      if (uc[cu] > 0 && uc[cu] < 1) { // interfacial cell
        if (type == 0 || type == 1 || type == 3) {
          ffvu[f] = R::GetLineFlux(fcn[cu], fca[cu], h, v, dt, d);
        } else if (type == 2) { // Lagrange Explicit
          Scal vu = (v > 0. ? fcfm[cu] : fcfp[cu]) * vsc;
          ffvu[f] = R::GetLineFluxStr(fcn[cu], fca[cu], h, v, vu, dt, d);
        }
      } else { // pure cell
        ffvu[f] = v * uc[cu];
      }
    }

    FieldFace<Scal> ffu(m); // volume fraction on boundaries
    InterpolateB(uc, mfc, ffu, m);

    // override boundary upwind flux
    for (const auto& it : mfc) {
      IdxFace f = it.GetIdx();
      CondFace* cb = it.GetValue().get(); 
      Scal v = ffv[f];
      if ((cb->GetNci() == 0) != (v > 0.)) {
        ffvu[f] = v * ffu[f];
      }
    }

    // update volume fraction [i]
    for (auto c : m.Cells()) {
      auto w = bc.GetMIdx(c);
      const Scal lc = m.GetVolume(c);
      IdxFace fm = bf.GetIdx(w, md);
      IdxFace fp = bf.GetIdx(w + wd, md);
      // mixture cfl
      const Scal ds = (ffv[fp] - ffv[fm]) * vsc * dt / lc;
      // phase 2 cfl
      const Scal dl = (ffvu[fp] - ffvu[fm]) * dt / lc;
      if (type == 0) {        // plain
        uc[c] += -dl;
      } else if (type == 1) { // Euler Implicit
        uc[c] = (uc[c] - dl) / (1. - ds);
      } else if (type == 2) { // Lagrange Explicit
        uc[c] += uc[c] * ds - dl;
      } else if (type == 3) { // Weymouth
        uc[c] += (uc[c] < 0.5 ? 0 : 1) * ds - dl;
      }
    }
  }
  void Clip(FieldCell<Scal>& uc, Scal th) {
    for (auto c : m.Cells()) {
      Scal& u = uc[c];
      if (u < th) {
        u = 0;
      } else if (u > 1 - th) {
        u = 1;
      } else if (IsNan(u)) {
        u = 0;
      }
    }
  }
  void AdvAulisa(Sem& sem) {
    // directions, format: {dir EI, dir LE, ...}
    std::vector<size_t> dd;
    Scal vsc; // scaling factor for ffv
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
        auto& uc = fcu_.iter_curr;
        Sweep(uc, d, *owner_->ffv_, fcn_, fca_, mfc_,
              id % 2 == 0 ? 1 : 2, vsc, fcfm_, fcfp_,
              owner_->GetTimeStep(), m);
        Clip(uc, par->clipth);
        m.Comm(&uc);
      }
      if (par->bcc_reflect && sem.Nested("reflect")) {
        BcReflect(fcu_.iter_curr, mfc_, par->bcc_fill, m);
      }
      if (sem.Nested("reconst")) {
        // Compute fcn_,fca_ [s] , fci_ [a]
        Rec(fcu_.iter_curr);
      }
    }
  }
  void MakeIteration() {
    auto sem = m.GetSem("iter");
    if (sem("init")) {
      auto& uc = fcu_.iter_curr;
      const Scal dt = owner_->GetTimeStep();
      auto& fcs = *owner_->fcs_;
      for (auto c : m.Cells()) {
        uc[c] = 
            fcu_.time_prev[c] +  // previous time step
            dt * fcs[c]; // source
      }
    }

    AdvAulisa(sem);

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
    if (par->bcc_reflect && sem("reflect")) {
      // --> fca [a], fcn [a]
      BcReflect(fca_, mfc_, Scal(0), m);
      BcReflect(fcn_, mfc_, Vect(0), m);
      // --> reflected fca [a], fcn [a]
    }
    if (par->curvgrad && sem("curv")) {
      // Curvature from gradient of volume fraction
      auto ffu = Interpolate(fcu_.iter_curr, mfc_, m); // [s]
      auto fcg = Gradient(ffu, m); // [s]
      auto ffg = Interpolate(fcg, mfvz_, m); // [i]

      fck_.Reinit(m, GetNan<Scal>()); // curvature [i]
      for (auto c : m.Cells()) {
        if (!fci_[c]) {
          continue;
        }
        Scal s = 0.;
        for (auto q : m.Nci(c)) {
          IdxFace f = m.GetNeighbourFace(c, q);
          auto& g = ffg[f];
          auto n = g / g.norm();  // inner normal
          s += -n.dot(m.GetOutwardSurface(c, q));
        }
        fck_[c] = s / m.GetVolume(c);
      }
    }
    if (!par->curvgrad) {
      auto& uc = fcu_.iter_curr;
      if (sem("diff2")) {
        CalcDiff2(uc, fcud2_);
        m.Comm(&fcud2_);
      }
      if (sem("diff4")) {
        CalcDiff4(uc, fcud2_, fcud4_);
        m.Comm(&fcud4_);
      }
      /*
      if (sem("diff6")) {
        CalcDiff6(uc, fcud4_, fcud6_);
        m.Comm(&fcud6_);
      }
      */
      if (sem("height")) {
        //UNormal<M>::CalcHeight(
        //    m, uc, fcud2_, fcud4_, fcud6_, par->dim, fch_);
        UNormal<M>::CalcHeight(m, uc, fcud2_, fcud4_, par->dim, fch_);
        //UNormal<M>::CalcHeight(m, uc, fcud2_, par->dim, fch_);
        m.Comm(&fch_);
      }
      if (sem("curvcomm")) {
        UNormal<M>::CalcCurvHeight(m, uc, fch_, fcn_, par->dim, fck_);
        m.Comm(&fck_);
      }
    }
    if (par->part && sem.Nested("part")) {
      const FieldCell<Scal>* fccl(nullptr);
      psm_->Part(&fcu_.iter_curr, &fca_, &fcn_, &fci_, fccl, &fck_, mfc_);
    }
    if (sem.Nested("dump")) {
      Dump();
    }
  }


  Owner* owner_;
  std::shared_ptr<Par> par;
  M& m;

  LayersData<FieldCell<Scal>> fcu_;
  MapFace<std::shared_ptr<CondFace>> mfc_;
  MapFace<std::shared_ptr<CondFace>> mfvz_; // zero-derivative bc for Vect

  FieldCell<Scal> fca_; // alpha (plane constant)
  FieldCell<Vect> fcn_; // n (normal to plane)
  FieldCell<Scal> fck_; // curvature from height functions
  FieldCell<bool> fci_; // interface mask (1: contains interface)
  FieldCell<Vect> fcud2_; // volume fraction difference double
  FieldCell<Vect> fcud4_; // volume fraction difference quad
  FieldCell<Vect> fcud6_; // volume fraction difference 6th
  FieldCell<Vect> fch_; // curvature from height functions
  size_t count_ = 0; // number of MakeIter() calls, used for splitting

  std::unique_ptr<PartStrMeshM<M>> psm_;
  // tmp for MakeIteration, volume flux copied to cells
  FieldCell<Scal> fcfm_, fcfp_;
  UVof<M> uvof_;
};

template <class M_>
Vof<M_>::Vof(
    M& m, const FieldCell<Scal>& fcu,
    const MapFace<std::shared_ptr<CondFace>>& mfc,
    const FieldFace<Scal>* ffv, const FieldCell<Scal>* fcs,
    double t, double dt, std::shared_ptr<Par> par)
    : AdvectionSolver<M>(t, dt, m, ffv, fcs)
    , imp(new Imp(this, fcu, mfc, par))
{}

template <class M_>
Vof<M_>::~Vof() = default;

template <class M_>
auto Vof<M_>::GetPar() -> Par* {
  return imp->par.get();
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
auto Vof<M_>::GetField(Layers l) const -> const FieldCell<Scal>& {
  return imp->fcu_.Get(l);
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
auto Vof<M_>::GetHeight() const -> const FieldCell<Vect>& {
  return imp->fch_;
}

template <class M_>
auto Vof<M_>::GetCurv() const -> const FieldCell<Scal>& {
  return imp->par->part_k ? imp->psm_->GetCurv(0) : imp->fck_;
}

template <class M_>
void Vof<M_>::PostStep() {
  return imp->PostStep();
}


// curvature from height function
template <class M_>
auto Vof<M_>::GetCurvH() const -> const FieldCell<Scal>& {
  return imp->fck_;
}

// curvature from particles
template <class M_>
auto Vof<M_>::GetCurvP() const -> const FieldCell<Scal>& {
  return imp->psm_->GetCurv(0);
}

} // namespace solver
