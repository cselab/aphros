#pragma once

#include <exception>
#include <fstream>
#include <array>
#include <memory>
#include <limits>

#include "vof.h"
#include "geom/block.h"
#include "dump/vtk.h"
#include "reconst.h"
#include "normal.h"
#include "debug/isnan.h"
#include "partstr.h"

namespace solver {

template <class T>
class Multi {
 public:
  Multi() : d_(1) {}
  Multi(size_t n) : d_(n) {}
  // cast to pointer
  template <class U>
  Multi(Multi<U>& u) : d_(u.size()) {
    for (size_t i = 0; i < u.size(); ++i) {
      d_[i] = &u[i];
    }
  }
  // cast to const
  template <class U>
  Multi(const Multi<U*>& u) : d_(u.size()) {
    for (size_t i = 0; i < u.size(); ++i) {
      d_[i] = const_cast<const U*>(u[i]);
    }
  }
  T& operator[](size_t i) {
    return d_[i];
  }
  const T& operator[](size_t i) const {
    return d_[i];
  }
  size_t size() const {
    return d_.size();
  }
  void resize(size_t n) {
    d_.resize(n);
  }
  std::vector<T>& data() {
    return d_;
  }
  const std::vector<T>& data() const {
    return d_;
  }
  void InitAll(const T& u) {
    for (auto& a : d_) {
      a = u;
    }
  }
  std::vector<T*> GetPtr() {
    std::vector<T*> r;
    for (auto& a : d_) {
      r.push_bask(&a);
    }
    return r;
  }
  std::vector<const T*> GetConstPtr() const {
    std::vector<const T*> r;
    for (auto& a : d_) {
      r.push_bask(&a);
    }
    return r;
  }

 private:
  std::vector<T> d_;
};

// Copies values from layer 0 to layers [1:] in cells with mask=0
// in stencil around each cell with mask=1.
// mfc: fields to propagate
// mask: if true, copy values from layer 0 in stencil 1.
template <class M>
void Propagate(Multi<FieldCell<typename M::Scal>*> mfc,
               const Multi<FieldCell<bool>*>& mask, M& m) {

  if (!mfc.size()) {
    return;
  }

  using MIdx = typename M::MIdx;
  auto& bc = m.GetIndexCells();
  const int sw = 2; // stencil halfwidth
  const int sn = sw * 2 + 1;
  // block of offsets
  GBlock<IdxCell, M::dim> bo(MIdx(-sw), MIdx(sn)); 

  auto& fc0 = *mfc[0];
  for (size_t i = 1; i < mfc.size(); ++i) {
    auto& fc = *mfc[i];
    auto& ms = *mask[i];
    for (auto c : m.Cells()) {
      if (ms[c]) {
        MIdx w = bc.GetMIdx(c);
        for (MIdx wo : bo) {
          IdxCell cn = bc.GetIdx(w + wo); 
          if (!ms[cn]) {
            fc[cn] = fc0[cn];
          }
        }
      }
    }
  }
}

// sw: stencil halfwidth
template <class M>
void Propagate(FieldCell<bool>& mask, size_t sw, M& m) {
  using MIdx = typename M::MIdx;
  auto& bc = m.GetIndexCells();
  const int sn = sw * 2 + 1;
  // block of offsets
  GBlock<IdxCell, M::dim> bo(MIdx(-sw), MIdx(sn)); 

  auto ms = mask;
  for (auto c : m.Cells()) {
    if (ms[c]) {
      MIdx w = bc.GetMIdx(c);
      for (MIdx wo : bo) {
        IdxCell cn = bc.GetIdx(w + wo);
        mask[cn] = true;
      }
    }
  }
}

template <class T>
Multi<FieldCell<T>*> GetLayer(Multi<LayersData<FieldCell<T>>>& u, Layers l) {
  Multi<FieldCell<T>*> r(u.size());
  for (size_t i = 0; i < u.size(); ++i) {
    r[i] = &u[i].Get(l);
  }
  return r;
}

template <class M>
FieldCell<bool> And(const FieldCell<bool>& u, 
                   const FieldCell<bool>& v, const M& m) {
  FieldCell<bool> r = u;
  for (auto c : m.AllCells()) {
    r[c] = r[c] && v[c];
  }
  return r;
}

template <class M_>
struct Vof<M_>::Imp {
  using Owner = Vof<M_>;
  using R = Reconst<Scal>;
  using PS = PartStr<Scal>;
  using PSM = PartStrMesh<M>;
  static constexpr size_t dim = M::dim;
  using Vect2 = GVect<Scal, 2>;

  Imp(Owner* owner, const FieldCell<Scal>& fcu,
      const MapFace<std::shared_ptr<CondFace>>& mfc, 
      std::shared_ptr<Par> par)
      : owner_(owner), par(par), m(owner_->m)
      , mfc_(mfc), fck_(m, 0)
      , fch_(m, Vect(0))
  {
    fcu_.resize(2);
    fcn_.resize(fcu_.size());
    fca_.resize(fcu_.size());
    fci_.resize(fcu_.size());
    fcm_.resize(fcu_.size());

    fcn_.InitAll(FieldCell<Vect>(m, Vect(0)));
    fca_.InitAll(FieldCell<Scal>(m, 0));
    fcus_.time_curr = fcu;
    fcu_[0].time_curr.Reinit(m, 0);
    fcu_[1].time_curr.Reinit(m, 0);
    fcm_.InitAll(FieldCell<bool>(m, false));
    for (auto c : m.AllCells()) {
      fcm_[0][c] = true;
      if (m.GetCenter(c)[1] < 0.5) {
        fcu_[0].time_curr[c] = fcu[c];
      } else {
        fcu_[1].time_curr[c] = fcu[c];
        fcm_[1][c] = true;
      }
    }
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

  // Dump cut polygons
  void DumpPoly() {
    auto sem = m.GetSem("dumppoly");
    if (sem("local")) {
      dl_.clear();
      dlc_.clear();
      dll_.clear();
      auto h = m.GetCellSize();
      for (size_t i = 0; i < fcu_.size(); ++i) {
        auto& fcn = fcn_[i];
        auto& fca = fca_[i];
        auto& fci = fci_[i];
        auto& fcu = fcu_[i].iter_curr;
        for (auto c : m.Cells()) {
          Scal u = fcu[c];
          if (IsNan(u) || IsNan(fcn[c]) || IsNan(fca[c])) {
            continue;
          }
          const Scal th = par->poly_intth;
          if (fci[c] && u > th && u < 1. - th) {
            dl_.push_back(R::GetCutPoly(m.GetCenter(c), fcn[c], fca[c], h));
            dlc_.push_back(m.GetHash(c));
            dll_.push_back(i);
          }
        }
      }
      using TV = typename M::template OpCatVT<Vect>;
      m.Reduce(std::make_shared<TV>(&dl_));
      using TS = typename M::template OpCatT<Scal>;
      m.Reduce(std::make_shared<TS>(&dlc_));
    }
    if (sem("write")) {
      if (m.IsRoot()) {
        std::string fn = GetDumpName("s", ".vtk", par->dmp->GetN());
        std::cout << std::fixed << std::setprecision(8)
            << "dump" 
            << " t=" << owner_->GetTime() + owner_->GetTimeStep()
            << " to " << fn << std::endl;
        WriteVtkPoly(fn, dl_, {&dlc_, &dll_}, {"c", "l"}, 
            "Reconstructed linear interface", true);
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
  void Rec(const Multi<FieldCell<Scal>*>& uc) {
    auto sem = m.GetSem("rec");
    if (sem("local")) {
      DetectInterface(uc);
      for (size_t i = 0; i < uc.size(); ++i) {
        auto& fcn = fcn_[i];
        auto& fca = fca_[i];
        auto& fci = fci_[i];
        auto& fcm = fcm_[i];
        auto& fcu = *uc[i];
        if (par->bcc_reflect) {
          BcReflect(fcu, mfc_, par->bcc_fill, m);
        }
        // Compute normal and curvature [s]
        //Propagate(GetLayer(fcu_, Layers::time_curr), fcm_, m);
        UNormal<M>::CalcNormal(m, fcu, And(fci, fcm, m), par->dim, fcn);
        auto h = m.GetCellSize();
        // Reconstruct interface [s]
        for (auto c : m.SuCells()) {
          if (fcm[c]) {
            fca[c] = R::GetLineA(fcn[c], fcu[c], h);
          }
        }
        if (par->bcc_reflect) {
          BcReflect(fca, mfc_, Scal(0), m);
          BcReflect(fcn, mfc_, Vect(0), m);
        }
      }
    }
  }
  void DetectInterface(const Multi<const FieldCell<Scal>*>& uc) {
    for (size_t i = 0; i < uc.size(); ++i) {
      auto& fci = fci_[i];
      auto& fcm = fcm_[i];
      auto& fcu = *uc[i];
      fci.Reinit(m, false);
      // volume fraction different from 0 or 1
      for (auto c : m.AllCells()) {
        if (fcm[c]) {
          Scal u = fcu[c];
          if (u > 0. && u < 1.) {
            fci[c] = true;
          }
        }
      }
      // neighbour cell has different value but both are 0 or 1
      for (auto f : m.SuFaces()) {
        IdxCell cm = m.GetNeighbourCell(f, 0);
        IdxCell cp = m.GetNeighbourCell(f, 1);
        if (fcm[cm] || fcm[cp]) {
          Scal um = fcu[cm];
          Scal up = fcu[cp];
          if ((um == 0. || um == 1.) && (up == 0. || up == 1.) && (um != up)) {
            fci[cm] = true;
            fci[cp] = true;
          }
        }
      }
      fcm.Reinit(m, false);
      for (auto c : m.AllCells()) {
        if (fcu[c] > 0) {
          fcm[c] = true;
        }
      }
      Propagate(fcm, 1, m);
    }
  }
  void StartStep() {
    auto sem = m.GetSem("start");
    if (sem("rotate")) {
      owner_->ClearIter();
      for (auto& u : fcu_.data()) {
        u.time_prev = u.time_curr;
        u.iter_curr = u.time_prev;
      }
      fcus_.time_prev = fcus_.time_curr;
      fcus_.iter_curr = fcus_.time_prev;
    }

    if (owner_->GetTime() == 0.) {
      if (sem.Nested("reconst")) {
        Rec(GetLayer(fcu_, Layers::time_curr));
      }
    }
  }
  void Dump() {
    auto sem = m.GetSem("iter");
    bool dm = par->dmp->Try(owner_->GetTime(),
                            owner_->GetTimeStep());
    if (par->dumppoly && dm && sem.Nested("dumppoly")) {
      DumpPoly();
    }
    if (par->dumppart && dm && sem.Nested("part-dump")) {
      psm_->DumpParticles(fca_[0], fcn_[0], par->dmp->GetN(), owner_->GetTime());
    }
    if (par->dumppartinter && dm && sem.Nested("partinter-dump")) {
      psm_->DumpPartInter(fca_[0], fcn_[0], par->dmp->GetN(), owner_->GetTime());
    }
  }
  void MakeIteration() {
    auto sem = m.GetSem("iter");
    if (sem("init")) {
      for (size_t i = 0; i < fcu_.size(); ++i) {
        auto& fcm = fcm_[i];
        auto& fcu = fcu_[i];
        auto& uc = fcu.iter_curr;
        const Scal dt = owner_->GetTimeStep();
        auto& fcs = *owner_->fcs_;
        for (auto c : m.Cells()) {
          if (fcm[c]) {
            uc[c] = 
                fcu.time_prev[c] +  // previous time step
                dt * fcs[c]; // source
          }
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
      }
      for (size_t i = 0; i < fcu_.size(); ++i) {
        if (sem("adv")) {
          Dir md(d); // direction as Dir
          MIdx wd(md); // offset in direction d
          auto& bc = m.GetIndexCells();
          auto& bf = m.GetIndexFaces();
          auto h = m.GetCellSize();
          auto& ffv = *owner_->ffv_; // [f]ield [f]ace [v]olume flux
          const Scal dt = owner_->GetTimeStep();

          FieldFace<Scal> ffvu(m); // flux: volume flux * field

          ffi_.Reinit(m);
          auto& fcu = fcu_[i].iter_curr;
          auto& fcn = fcn_[i];
          auto& fca = fca_[i];
          auto& fci = fci_[i];
          for (auto f : m.Faces()) {
            auto p = bf.GetMIdxDir(f);
            Dir df = p.second;

            if (df != md) {
              continue;
            }

            // mixture flux
            const Scal v = ffv[f] * vsc;
            // upwind cell
            IdxCell cu = m.GetNeighbourCell(f, v > 0. ? 0 : 1);
            if (fci[cu]) { // cell contains interface, flux from reconstruction
              if (id % 2 == 0) { // Euler Implicit
                // phase 2 flux
                ffvu[f] = R::GetLineFlux(fcn[cu], fca[cu], h, v, dt, d);
              } else { // Lagrange Explicit
                // XXX: if id % 2 == 1, fcfm_ and fcfp_ contain fluxes
                // upwind mixture flux
                Scal vu = (v > 0. ? fcfm_[cu] : fcfp_[cu]) * vsc;
                // phase 2 flux
                ffvu[f] = R::GetLineFluxStr(fcn[cu], fca[cu], h, v, vu, dt, d);
              }
              ffi_[f] = true;
            } else {
              ffvu[f] = v * fcu[cu];
              ffi_[f] = false;
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
              ffi_[f] = true;
            }
          }

          for (auto c : m.Cells()) {
            auto w = bc.GetMIdx(c);
            const Scal lc = m.GetVolume(c);
            // faces
            IdxFace fm = bf.GetIdx(w, md);
            IdxFace fp = bf.GetIdx(w + wd, md);
            if (!ffi_[fm] && !ffi_[fp]) {
              continue;
            }
            // mixture fluxes
            const Scal vm = ffv[fm] * vsc;
            const Scal vp = ffv[fp] * vsc;
            // mixture cfl
            const Scal sm = vm * dt / lc;
            const Scal sp = vp * dt / lc;
            const Scal ds = sp - sm;
            // phase 2 fluxes
            Scal qm = ffvu[fm];
            Scal qp = ffvu[fp];
            // phase 2 cfl
            const Scal lm = qm * dt / lc;
            const Scal lp = qp * dt / lc;
            const Scal dl = lp - lm;
            if (id % 2 == 0) { // Euler Implicit
              fcu[c] = (fcu[c] - dl) / (1. - ds);
            } else { // Lagrange Explicit
              fcu[c] = fcu[c] * (1. + ds) - dl;
            }
          }

          // clip
          const Scal th = par->clipth;
          for (auto c : m.Cells()) {
            Scal& u = fcu[c];
            if (u < th) {
              u = 0.;
            } else if (u > 1. - th) {
              u = 1.;
            } else if (IsNan(u)) {
              u = 0.;
            }
          }
          m.Comm(&fcu);
        }
      }
      if (sem.Nested("reconst")) {
        Rec(GetLayer(fcu_, Layers::iter_curr));
      }
    }
    if (sem("stat")) {
      owner_->IncIter();
      ++count_;
      auto l = Layers::iter_curr;
      auto& fcus = fcus_.Get(l);
      fcus.Reinit(m, 0);
      for (size_t i = 0; i < fcu_.size(); ++i) {
        auto& fcm = fcm_[i];
        auto& fcu = fcu_[i].Get(l);
        for (auto c : m.AllCells()) {
          if (fcm[c]) {
            fcus[c] += fcu[c];
          }
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
      auto ffu = Interpolate(fcu_[0].iter_curr, mfc_, m); // [s]
      auto fcg = Gradient(ffu, m); // [s]
      auto ffg = Interpolate(fcg, mfvz_, m); // [i]

      fck_.Reinit(m, GetNan<Scal>()); // curvature [i]
      auto& fci = fci_[0];
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
        fck_[c] = s / m.GetVolume(c);
      }
    }
    if (!par->curvgrad) {
      auto& uc = fcu_[0].iter_curr;
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
        UNormal<M>::CalcCurvHeight(m, uc, fch_, fcn_[0], par->dim, fck_);
        m.Comm(&fck_);
      }
    }
    if (par->part && sem.Nested("part")) {
      psm_->Part(fcu_[0].iter_curr, fca_[0], fcn_[0], fci_[0], fck_, mfc_);
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
  MapFace<std::shared_ptr<CondFace>> mfc_;
  MapFace<std::shared_ptr<CondFace>> mfvz_; // zero-derivative bc for Vect

  Multi<FieldCell<Scal>> fca_; // alpha (plane constant)
  Multi<FieldCell<Vect>> fcn_; // n (normal to plane)
  FieldCell<Scal> fck_; // curvature from height functions
  Multi<FieldCell<bool>> fci_; // interface mask (1: contains interface)
  Multi<FieldCell<bool>> fcm_; // layer mask
  FieldFace<bool> ffi_; // interface mask (1: upwind cell contains interface)
  FieldCell<Vect> fcud2_; // volume fraction difference double
  FieldCell<Vect> fcud4_; // volume fraction difference quad
  FieldCell<Vect> fcud6_; // volume fraction difference 6th
  FieldCell<Vect> fch_; // curvature from height functions
  size_t count_ = 0; // number of MakeIter() calls, used for splitting

  std::vector<std::vector<Vect>> dl_; // dump poly
  std::vector<Scal> dlc_; // dump poly cell
  std::vector<Scal> dll_; // dump poly layer

  std::unique_ptr<PartStrMesh<M>> psm_;
  // tmp for MakeIteration, volume flux copied to cells
  FieldCell<Scal> fcfm_, fcfp_;
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
  return imp->fcus_.Get(l);
}

template <class M_>
auto Vof<M_>::GetField(Layers l, size_t i) const -> const FieldCell<Scal>& {
  return imp->fcu_[i].Get(l);
}

template <class M_>
auto Vof<M_>::GetField(size_t i) const -> const FieldCell<Scal>& {
  return imp->fcu_[i].Get(Layers::time_curr);
}

template <class M_>
auto Vof<M_>::GetAlpha(size_t i) const -> const FieldCell<Scal>& {
  (void) i;
  return imp->fca_[i];
}

template <class M_>
auto Vof<M_>::GetMask(size_t i) const -> const FieldCell<bool>& {
  return imp->fcm_[i];
}

template <class M_>
size_t Vof<M_>::GetNumLayers() const {
  return imp->fcu_.size();
}

template <class M_>
auto Vof<M_>::GetNormal(size_t i) const -> const FieldCell<Vect>& {
  return imp->fcn_[i];
}

template <class M_>
auto Vof<M_>::GetHeight() const -> const FieldCell<Vect>& {
  return imp->fch_;
}

template <class M_>
auto Vof<M_>::GetCurv() const -> const FieldCell<Scal>& {
  return imp->par->part_k ? imp->psm_->GetCurv() : imp->fck_;
}

template <class M_>
auto Vof<M_>::GetCurv(size_t i) const -> const FieldCell<Scal>& {
  (void) i;
  return imp->par->part_k ? imp->psm_->GetCurv() : imp->fck_;
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
  return imp->psm_->GetCurv();
}

} // namespace solver
