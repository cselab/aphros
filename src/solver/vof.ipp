// Created by Petr Karnakov on 30.07.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <array>
#include <exception>
#include <fstream>
#include <limits>
#include <memory>

#include "approx.h"
#include "approx_eb.h"
#include "debug/isnan.h"
#include "geom/block.h"
#include "normal.h"
#include "reconst.h"
#include "trackerm.h"
#include "util/convdiff.h"
#include "util/vof.h"

#include "vof.h"

template <class EB_>
struct Vof<EB_>::Imp {
  using Owner = Vof<EB_>;
  using R = Reconst<Scal>;
  using TRM = Trackerm<M>;
  static constexpr size_t dim = M::dim;
  using Sem = typename M::Sem;
  using MIdx = typename M::MIdx;
  using UEB = UEmbed<M>;

  Imp(Owner* owner, const EB& eb0, const FieldCell<Scal>& fcu,
      const FieldCell<Scal>& fccl, Par par0)
      : owner_(owner)
      , par(par0)
      , m(owner_->m)
      , eb(eb0)
      , layers(1)
      , mebc_(owner->mebc_)
      , fccl_(fccl)
      , fcim_(m, TRM::Pack(MIdx(0)))
      , fcim_unpack_(m, MIdx(0))
      , fca_(m, GetNan<Scal>())
      , fcn_(m, GetNan<Vect>())
      , fci_(m, false) {
    par.dim = std::min(par.dim, M::dim);

    fcu_.time_curr = fcu;

    UpdateBc(mebc_);

    for (auto c : m.AllCells()) {
      if (fcu[c] == 0) {
        fccl_[c] = kClNone;
      }
    }
  }
  void UpdateBc(const MapEmbed<BCondAdvection<Scal>>& mebc) {
    std::tie(me_vf_, me_cl_, me_im_, me_n_, me_a_) =
        UVof<M>::GetAdvectionBc(m, mebc);
  }
  // Computes normal and plane constant in interfacial cells.
  // uc: volume fraction to compute plane constant [a]
  // Output:
  // fci_: interface mask
  // fca_: plane constant [s]
  // fcn_: normal [s]
  void ReconstPlanes(const FieldCell<Scal>& uc) {
    DetectInterface(uc);
    UNormal<M>::CalcNormal(m, uc, fci_, par.dim, fcn_);
    auto h = m.GetCellSize();
    for (auto c : eb.SuCells()) {
      if (fci_[c]) {
        fca_[c] = R::GetLineA(fcn_[c], uc[c], h);
      } else {
        fca_[c] = GetNan<Scal>();
      }
    }
  }
  void DetectInterface(const FieldCell<Scal>& uc) {
    fci_.Reinit(m, false);
    for (auto c : eb.AllCells()) {
      const Scal u = uc[c];
      if (u > 0 && u < 1) {
        fci_[c] = true;
      }
    }
  }
  // Extrapolates volume fraction to halo and excluded cells
  // with a linear least-squares fit.
  // Output:
  // fcu: updated halo ad excluded cells.
  void ExtrapolateLinear(Sem& sem, FieldCell<Scal>& fcu) {
    if (sem("extrap-linear")) {
      // Fill volume fraction in halo cells from quadratic extrapolation
      for (auto p : mebc_.GetMapFace()) {
        const IdxFace f = p.first;
        auto& bc = mebc_.at(f);
        auto cc = m.GetCellColumn(f, bc.nci);
        const IdxCell cm = cc[1];
        const IdxCell cp = cc[2];
        const IdxCell cpp = cc[3];
        fcu[cm] = 2 * fcu[cp] - fcu[cpp];
      }
      if (eb.kIsEmbed) {
        for (auto c : m.Cells()) {
          if (eb.IsExcluded(c)) {
            std::vector<Vect> xx;
            std::vector<Scal> uu;
            for (auto cn : m.Stencil5(c)) {
              if (!eb.IsExcluded(c)) {
                xx.push_back(m.GetCenter(cn));
                uu.push_back(fcu[cn]);
              }
            }
            if (xx.size()) {
              const auto p = ULinearFit<Vect>::FitLinear(xx, uu);
              fcu[c] =
                  ULinearFit<typename EB::Vect>::EvalLinear(p, m.GetCenter(c));
            }
          }
        }
        m.Comm(&fcu);
      }
    }
  }
  // Extrapolates volume fraction, normal, and plane constant
  // to excluded cells by PLIC plane from the nearest regular cell.
  void ExtrapolatePlic(Sem& sem, FieldCell<Scal>& uc) {
    if (sem("extrap-plic")) {
      // Fill halo cells from PLIC extrapolation.
      for (auto p : mebc_.GetMapFace()) {
        const IdxFace f = p.first;
        auto& bc = mebc_.at(f);
        auto cc = m.GetCellColumn(f, bc.nci);
        const IdxCell cm = cc[1];
        const IdxCell cp = cc[2];
        if (fci_[cp]) {
          fcn_[cm] = fcn_[cp];
          fca_[cm] =
              fca_[cp] - (m.GetCenter(cm) - m.GetCenter(cp)).dot(fcn_[cp]);
          uc[cm] = R::GetLineU(fcn_[cp], fca_[cp], m.GetCellSize());
          fccl_[cm] = fccl_[cp];
          fci_[cm] = (uc[cm] > 0 && uc[cm] < 1);
        } else {
          uc[cm] = uc[cp];
          fci_[cm] = false;
        }
      }
      if (eb.kIsEmbed) {
        for (auto c : eb.Cells()) {
          if (eb.IsCut(c)) {
            const IdxCell cn = eb.GetRegularNeighbor(c);
            if (!eb.IsRegular(cn)) {
              fci_[c] = false;
              uc[c] = 0;
            } else if (fci_[cn]) {
              fcn_[c] = fcn_[cn];
              fca_[c] =
                  fca_[cn] - (m.GetCenter(c) - m.GetCenter(cn)).dot(fcn_[cn]);
              fci_[c] = fci_[cn];
              uc[c] = R::GetLineU(fcn_[c], fca_[c], m.GetCellSize());
            } else {
              fci_[c] = false;
              uc[c] = uc[cn];
            }
          }
        }
        m.Comm(&uc);
        m.Comm(&fca_);
        m.Comm(&fcn_);
      }
    }
  }
  void StartStep() {
    auto sem = m.GetSem("start");
    if (sem("rotate")) {
      owner_->ClearIter();
      fcu_.time_prev = fcu_.time_curr;
      fcu_.iter_curr = fcu_.time_prev;
      UpdateBc(mebc_);
    }

    if (owner_->GetTime() == 0.) {
      auto& fcu = fcu_.time_curr;
      if (par.extrapolate_boundaries) {
        ExtrapolateLinear(sem, fcu);
      }
      if (sem("reconst")) {
        ReconstPlanes(fcu);
      }
    }
  }
  enum class SweepType {
    plain, // sum of fluxes
    EI, // Euler Implicit (aulisa2009)
    LE, // Lagrange Explicit (aulisa2009)
    weymouth, // sum of fluxes and divergence (weymouth2010)
  };
  // Makes advection sweep in one direction, updates uc [i]
  // uc: volume fraction [s]
  // dir: direction
  // ffv: mixture flux [i]
  // fcn,fca: normal and plane constant [s]
  // mebc: face conditions, nullptr to keep boundary fluxes
  // fcfm,fcfp: upwind mixture flux, required if type=LE [s]
  // fcuu: volume fraction for Weymouth div term
  // dt: time step
  // clipth: threshold for clipping, values outside [th,1-th] are clipped
  static void Sweep(
      FieldCell<Scal>& uc, size_t dir, const FieldFace<Scal>& ffv,
      FieldCell<Scal>& fccl, FieldCell<Scal>& fcim, const FieldCell<Vect>& fcn,
      const FieldCell<Scal>& fca, const MapEmbed<BCond<Scal>>* mebc,
      SweepType type, const FieldCell<Scal>* fcfm, const FieldCell<Scal>* fcfp,
      const FieldCell<Scal>* fcuu, Scal dt, Scal clipth, const EB& eb) {
    const auto& m = eb.GetMesh();
    const auto& indexc = m.GetIndexCells();
    const auto& indexf = m.GetIndexFaces();
    const MIdx globalsize = m.GetGlobalSize();
    const auto h = m.GetCellSize();
    const auto d = m.direction(dir);

    FieldFace<Scal> ffvu(m, 0); // phase 2 flux

    // compute fluxes [i] and propagate color to downwind cells
    for (auto f : eb.Faces()) {
      const auto p = indexf.GetMIdxDir(f);
      if (p.second.raw() != d) {
        continue;
      }

      // flux through face (maybe cut)
      const Scal v = ffv[f];
      // flux through full face that would give the same velocity
      const Scal v0 = v / eb.GetAreaFraction(f);
      const IdxCell c = m.GetCell(f, v > 0 ? 0 : 1); // upwind cell
      if (uc[c] > 0 && uc[c] < 1 && fcn[c].sqrnorm() > 0) { // interfacial cell
        switch (type) {
          case SweepType::plain:
          case SweepType::EI:
          case SweepType::weymouth: {
            ffvu[f] = R::GetLineFlux(fcn[c], fca[c], h, v0, dt, d);
            break;
          }
          case SweepType::LE: {
            const Scal vc = (v > 0 ? (*fcfm)[c] : (*fcfp)[c]);
            const Scal vc0 = vc / eb.GetAreaFraction(f);
            ffvu[f] = R::GetLineFluxStr(fcn[c], fca[c], h, v0, vc0, dt, d);
            break;
          }
        }
      } else { // pure cell or cell with undefined normal
        ffvu[f] = v0 * uc[c];
      }

      // propagate color to downwind cell if empty
      if (fccl[c] != kClNone) {
        const IdxCell cd = m.GetCell(f, v > 0 ? 1 : 0); // downwind cell
        if (fccl[cd] == kClNone) {
          fccl[cd] = fccl[c];
          const MIdx w = indexc.GetMIdx(c);
          MIdx im = TRM::Unpack(fcim[c]);
          if (w[d] < 0) im[d] += 1;
          if (w[d] >= globalsize[d]) im[d] -= 1;
          fcim[cd] = TRM::Pack(im);
        }
      }
    }

    if (mebc) {
      // override flux in upwind boundaries
      const FieldFace<Scal> ffu = UEB::Interpolate(uc, *mebc, m);
      for (const auto& p : mebc->GetMapFace()) {
        const IdxFace f = p.first;
        const auto& bc = p.second;
        const Scal v = ffv[f];
        if ((bc.nci == 0) != (v > 0)) {
          ffvu[f] = v * ffu[f];
        }
      }
    }

    // update volume fraction [i]
    for (auto c : eb.CellsM()) {
      const auto fm = c.face(-d);
      const auto fp = c.face(d);
      // mixture cfl
      // add 1e-16 to avoid zero denominator for excluded faces
      const Scal ds = (ffv[fp] / (eb.GetAreaFraction(fp) + 1e-16) -
                       ffv[fm] / (eb.GetAreaFraction(fm) + 1e-16)) *
                      dt / c.volume;
      // phase 2 cfl
      const Scal dl = (ffvu[fp] - ffvu[fm]) * dt / c.volume;
      auto& u = uc[c];
      switch (type) {
        case SweepType::plain: {
          u += -dl;
          break;
        }
        case SweepType::EI: {
          u = (u - dl) / (1 - ds);
          break;
        }
        case SweepType::LE: {
          u += u * ds - dl;
          break;
        }
        case SweepType::weymouth: {
          u += (*fcuu)[c] * ds - dl;
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
        fcim[c] = TRM::Pack(MIdx(0));
      }
    }
  }
  // Removes orphan fragments, those for which the volume fraction
  // in the 3x3x3 stencil does not exceed the threshold.
  // fcu: volume fractions
  // fcim: packed image vectors
  // filterth: threshold for detecting orphan fragments
  static void FilterOrphan(
      FieldCell<Scal>& fcu, FieldCell<Scal>& fccl, FieldCell<Scal>& fcim,
      Scal filterth, const EB& eb) {
    for (auto c : eb.Cells()) {
      bool orphan = true;
      for (auto cn : eb.Stencil(c)) {
        if (fcu[cn] >= filterth) {
          orphan = false;
          break;
        }
      }
      if (orphan) {
        fcu[c] = 0;
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
      BcApply(uc, me_vf_, m);
      BcApply(fccl, me_cl_, m);
      BcApply(fcim, me_im_, m);
    }
    if (par.extrapolate_boundaries) {
      ExtrapolateLinear(sem, uc);
    }
    if (sem("reconst")) {
      ReconstPlanes(uc);
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
      const auto d = m.direction(dd[id]);
      if (sem("copyface")) {
        if (id % 2 == 1) { // copy fluxes for Lagrange Explicit step
          auto& ffv =
              owner_->fev_->GetFieldFace(); // [f]ield [f]ace [v]olume flux
          fcfm_.Reinit(m);
          fcfp_.Reinit(m);
          for (auto c : eb.CellsM()) {
            fcfm_[c] = ffv[c.face(-d)];
            fcfp_[c] = ffv[c.face(d)];
          }
          m.Comm(&fcfm_);
          m.Comm(&fcfp_);
        }
      }
      auto& uc = fcu_.iter_curr;
      if (sem("sweep")) {
        Sweep(
            uc, d, owner_->fev_->GetFieldFace(), fccl_, fcim_, fcn_, fca_,
            &me_vf_, id % 2 == 0 ? SweepType::EI : SweepType::LE, &fcfm_,
            &fcfp_, nullptr, owner_->GetTimeStep() * vsc, par.clipth, eb);
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
  void AdvPlain(Sem& sem, SweepType type) {
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
            uc, dd[id], owner_->fev_->GetFieldFace(), fccl_, fcim_, fcn_, fca_,
            &me_vf_, type, nullptr, nullptr, &fcuu_, owner_->GetTimeStep(),
            par.clipth, eb);
      }
      CommRec(sem, uc, fccl_, fcim_);
      if (par.extrapolate_boundaries) {
        ExtrapolatePlic(sem, uc);
      }
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
        FieldFace<Scal> ffv(m, 0);
        for (auto f : eb.Faces()) {
          const IdxCell cm = m.GetCell(f, 0);
          const IdxCell cp = m.GetCell(f, 1);
          ffv[f] = std::min(eb.GetVolume(cm), eb.GetVolume(cp)) * sgn *
                   par.sharpen_cfl;
        }
        // zero flux on boundaries
        for (const auto& it : me_vf_.GetMapFace()) {
          ffv[it.first] = 0;
        }
        Sweep(
            uc, d, ffv, fccl_, fcim_, fcn_, fca_, &me_vf_, SweepType::weymouth,
            nullptr, nullptr, &fcuu_, 1., par.clipth, eb);
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
      for (auto c : eb.Cells()) {
        auto um = fcu_.time_prev[c];
        uc[c] = um + dt * fcs[c] * (1 - um);
        fcuu_[c] = (uc[c] < 0.5 ? 0 : 1);
      }
    }

    using Scheme = typename Par::Scheme;
    switch (par.scheme) {
      case Scheme::plain:
        AdvPlain(sem, SweepType::plain);
        break;
      case Scheme::aulisa:
        AdvAulisa(sem);
        break;
      case Scheme::weymouth:
        AdvPlain(sem, SweepType::weymouth);
        break;
    }
    if (par.sharpen && sem.Nested("sharpen")) {
      Sharpen();
    }
    if (par.filterth > 0) {
      if (sem("filterorphan")) {
        FilterOrphan(fcu_.iter_curr, fccl_, fcim_, par.filterth, eb);
      }
      CommRec(sem, fcu_.iter_curr, fccl_, fcim_);
    }
    if (modifier_) {
      if (sem("modify")) {
        modifier_(fcu_.iter_curr, fccl_, eb);
      }
      CommRec(sem, fcu_.iter_curr, fccl_, fcim_);
    }
    if (sem("bcc_clear")) {
      UVof<M>::BcClear(fcu_.iter_curr, mebc_, m);
    }
    if (sem("resetcolor")) {
      auto& fcu = fcu_.iter_curr;
      fcclm = fccl_;
      for (auto c : eb.AllCells()) {
        fccl_[c] = (fcu[c] > 0.5 ? 0 : kClNone);
      }
      BcApply(fccl_, me_cl_, m);
    }
    if (sem.Nested()) {
      if (par.labeling) {
        par.labeling->Recolor(
            layers, &fcu_.iter_curr, &fccl_, &fcclm, me_cl_, m);
      } else {
        uvof_.Recolor(
            layers, &fcu_.iter_curr, &fccl_, &fcclm, par.clfixed, par.clfixed_x,
            par.coalth, me_cl_, par.verb, par.recolor_unionfind,
            par.recolor_reduce, par.recolor_grid, m);
      }
    }
    if (sem("propagate")) {
      auto& u = fcu_.iter_curr;
      auto& cl = fccl_;
      for (auto f : eb.Faces()) {
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
    modifier_ = nullptr;
  }
  void PostStep() {
    auto sem = m.GetSem("iter");
    if (sem("comm")) {
      // --> fcu [a]
      // --> fca [s], fcn [s]
      m.Comm(&fca_);
      m.Comm(&fcn_);
    }
    if (sem("local")) {
      // --> fca [a], fcn [a]
      BcApply(fca_, me_a_, m);
      BcApply(fcn_, me_n_, m);
      // --> reflected fca [a], fcn [a]

      // unpack image vector
      for (auto c : m.AllCells()) {
        fcim_unpack_[c] = TRM::Unpack(fcim_[c]);
      }
    }
  }
  void DumpInterface(
      std::string filename,
      std::vector<Multi<const FieldCell<Scal>*>> extra_fields,
      std::vector<std::string> extra_names) {
    typename UVof<M>::DumpPolyArgs args;
    args.layers = GRange<size_t>(1);
    args.fcu = &fcu_.time_curr;
    args.fccl = (FieldCell<Scal>*)(nullptr);
    args.fcn = &fcn_;
    args.fca = &fca_;
    args.fci = &fci_;
    args.filename = filename;
    args.time = owner_->GetTime();
    args.binary = par.vtkbin;
    args.merge = par.vtkmerge;
    args.extra_fields = extra_fields;
    args.extra_names = extra_names;
    uvof_.DumpPoly(args, m);
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
        BcReflectAll(fcut, me_vf_, m);
        BcReflectAll(fcclt, me_cl_, m);
      }
      if (par.dumppolymarch_fill >= 0) {
        BcMarchFill(fcut, par.dumppolymarch_fill, m);
      }
    }
    if (sem.Nested()) {
      uvof_.DumpPolyMarch(
          layers, &fcut, &fcclt, &fcn_, filename, owner_->GetTime(),
          par.vtkpoly, par.vtkbin, par.vtkmerge, par.vtkiso,
          par.dumppolymarch_fill >= 0 ? &fcut : nullptr, m);
    }
  }

  Owner* owner_;
  Par par;
  M& m;
  const EB& eb;
  GRange<size_t> layers;

  StepData<FieldCell<Scal>> fcu_;
  FieldCell<Scal> fcuu_; // volume fraction for Weymouth div term

  // boundary conditions
  const MapEmbed<BCondAdvection<Scal>>& mebc_; // advection
  MapEmbed<BCond<Scal>> me_vf_; // volume fraction
  MapEmbed<BCond<Scal>> me_cl_; // color
  MapEmbed<BCond<Scal>> me_im_; // image
  MapEmbed<BCond<Vect>> me_n_; // normal
  MapEmbed<BCond<Scal>> me_a_; // plane constant

  FieldCell<Scal> fccl_; // color
  FieldCell<Scal> fcim_; // image
  FieldCell<MIdx> fcim_unpack_; // image unpacked
  FieldCell<Scal> fca_; // alpha (plane constant)
  FieldCell<Vect> fcn_; // n (normal to plane)
  FieldCell<bool> fci_; // interface mask (1: contains interface)
  size_t count_ = 0; // number of MakeIter() calls, used for splitting

  // tmp for MakeIteration, volume flux copied to cells
  FieldCell<Scal> fcfm_, fcfp_;
  UVof<M> uvof_;
  std::function<void(FieldCell<Scal>& fcu, FieldCell<Scal>& fccl, const EB&)>
      modifier_;
};

template <class EB_>
Vof<EB_>::Vof(
    M& m_, const EB& eb, const FieldCell<Scal>& fcu,
    const FieldCell<Scal>& fccl, const MapEmbed<BCondAdvection<Scal>>& mebc,
    const FieldEmbed<Scal>* fev, const FieldCell<Scal>* fcs, double t,
    double dt, Par par)
    : AdvectionSolver<M>(t, dt, m_, fev, fcs, mebc)
    , imp(new Imp(this, eb, fcu, fccl, par)) {}

template <class EB_>
Vof<EB_>::~Vof() = default;

template <class EB_>
auto Vof<EB_>::GetEmbed() const -> const EB& {
  return imp->eb;
}

template <class EB_>
auto Vof<EB_>::GetPar() const -> const Par& {
  return imp->par;
}

template <class EB_>
void Vof<EB_>::SetPar(Par par) {
  imp->par = par;
}

template <class EB_>
void Vof<EB_>::StartStep() {
  imp->StartStep();
}

template <class EB_>
void Vof<EB_>::Sharpen() {
  imp->Sharpen();
}

template <class EB_>
void Vof<EB_>::MakeIteration() {
  imp->MakeIteration();
}

template <class EB_>
void Vof<EB_>::FinishStep() {
  imp->FinishStep();
}

template <class EB_>
auto Vof<EB_>::GetField(Step l) const -> const FieldCell<Scal>& {
  return imp->fcu_.Get(l);
}

template <class EB_>
auto Vof<EB_>::GetColor() const -> const FieldCell<Scal>& {
  return imp->fccl_;
}

template <class EB_>
auto Vof<EB_>::GetPlic() const -> Plic {
  return {imp->layers, &GetField(), &GetAlpha(), &GetNormal(), &GetMask(),
          nullptr,     &GetColor(), &GetImage(), imp->mebc_};
}

template <class EB_>
auto Vof<EB_>::GetMask() const -> const FieldCell<bool>& {
  return imp->fci_;
}

template <class EB_>
auto Vof<EB_>::GetAlpha() const -> const FieldCell<Scal>& {
  return imp->fca_;
}

template <class EB_>
auto Vof<EB_>::GetNormal() const -> const FieldCell<Vect>& {
  return imp->fcn_;
}

template <class EB_>
auto Vof<EB_>::GetImage() const -> const FieldCell<MIdx>& {
  return imp->fcim_unpack_;
}

template <class EB_>
void Vof<EB_>::PostStep() {
  return imp->PostStep();
}

template <class EB_>
void Vof<EB_>::DumpInterface(
    std::string fn,
    std::vector<Multi<const FieldCell<Scal>*>> extra_fields,
    std::vector<std::string> extra_names) const {
  return imp->DumpInterface(fn, extra_fields, extra_names);
}

template <class EB_>
void Vof<EB_>::DumpInterfaceMarch(std::string fn) const {
  return imp->DumpInterfaceMarch(fn);
}

template <class EB_>
void Vof<EB_>::AddModifier(
    std::function<void(FieldCell<Scal>& fcu, FieldCell<Scal>& fccl, const EB&)>
        func) {
  imp->modifier_ = func;
}
