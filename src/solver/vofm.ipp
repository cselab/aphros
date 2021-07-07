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
#include "dump/hdf.h"
#include "geom/block.h"
#include "multi.h"
#include "normal.h"
#include "reconst.h"
#include "trackerm.h"
#include "util/convdiff.h"
#include "util/filesystem.h"
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

  Imp(Owner* owner, const EB& eb_, const GRange<size_t> layers0,
      const Multi<const FieldCell<Scal>*>& fcu0,
      const Multi<const FieldCell<Scal>*>& fccl0, Par par_)
      : owner_(owner)
      , par(par_)
      , m(owner_->m)
      , eb(eb_)
      , layers(layers0)
      , fcuu_(layers, m, 0)
      , fccls_(m, kClNone)
      , fcn_(layers, m, GetNan<Vect>())
      , fca_(layers, m, GetNan<Scal>())
      , fci_(layers, m, false)
      , fccl_(fccl0)
      , fcim_(layers, m, TRM::Pack(MIdx(0)))
      , fcim_unpack_(layers, m, MIdx(0))
      , mebc_(owner_->mebc_) {
    par.dim = std::min(par.dim, M::dim);

    fcu0.assert_size(layers);
    fccl0.assert_size(layers);
    fcu_.time_curr = fcu0;

    CalcSum(
        layers, fcu_.Get(Step::time_curr), fccl_, fcus_.Get(Step::time_curr),
        fccls_, m);

    UpdateBc(mebc_);
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
  void UpdateBc(const MapEmbed<BCondAdvection<Scal>>& mebc) {
    std::tie(me_vf_, me_cl_, me_im_, me_n_, me_a_) =
        UVof<M>::GetAdvectionBc(m, mebc);
  }
  // reconstruct interface
  void ReconstPlanes(const Multi<FieldCell<Scal>*>& uc) {
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
        for (auto c : eb.SuCells()) {
          if (fci[c]) {
            const auto uu =
                GetStencilValues<Scal>(layers, uc, fccl_, c, fccl[c], m);
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
        if (u > 0 && u < 1) {
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
  // Extrapolates volume fraction to halo and excluded cells
  // with a linear least-squares fit.
  // Output:
  // fcu: updated halo ad excluded cells.
  void ExtrapolateLinear(
      Sem& sem, const Multi<FieldCell<Scal>*>& fcu,
      const Multi<FieldCell<Scal>*>& fccl) {
    if (sem("extrap-linear")) {
      // Fill volume fraction in halo cells from quadratic extrapolation
      for (auto p : mebc_.GetMapFace()) {
        const IdxFace f = p.first;
        auto& bc = mebc_.at(f);
        auto cc = m.GetCellColumn(f, bc.nci);
        const IdxCell cm = cc[1];
        const IdxCell cp = cc[2];
        const IdxCell cpp = cc[3];
        for (auto lp : layers) {
          const Scal clp = (*fccl[lp])[cp];
          Scal upp = 0;
          for (auto lpp : layers) {
            if ((*fccl[lpp])[cpp] == clp) {
              upp = (*fcu[lpp])[cpp];
            }
          }
          (*fcu[lp])[cm] = 2 * (*fcu[lp])[cp] - upp;
          (*fccl[lp])[cm] = (*fccl[lp])[cp];
        }
      }
    }
  }
  // Extrapolates volume fraction, normal, and plane constant
  // to excluded cells by PLIC plane from the nearest regular cell.
  void ExtrapolatePlic(
      Sem& sem, const Multi<FieldCell<Scal>*>& fcu,
      const Multi<FieldCell<Scal>*>& fccl) {
    if (sem("extrap-plic")) {
      // Fill halo cells from PLIC extrapolation.
      for (auto p : mebc_.GetMapFace()) {
        const IdxFace f = p.first;
        auto& bc = mebc_.at(f);
        auto cc = m.GetCellColumn(f, bc.nci);
        const IdxCell cm = cc[1];
        const IdxCell cp = cc[2];
        for (auto lp : layers) {
          if (fci_[lp][cp]) {
            fcn_[lp][cm] = fcn_[lp][cp];
            fca_[lp][cm] =
                fca_[lp][cp] -
                (m.GetCenter(cm) - m.GetCenter(cp)).dot(fcn_[lp][cp]);
            (*fcu[lp])[cm] =
                R::GetLineU(fcn_[lp][cp], fca_[lp][cp], m.GetCellSize());
            (*fccl[lp])[cm] = (*fccl[lp])[cp];
            fci_[lp][cm] = ((*fcu[lp])[cp] > 0 && (*fcu[lp])[cm] < 1);
          } else {
            (*fcu[lp])[cm] = (*fcu[lp])[cp];
            fci_[lp][cm] = false;
          }
        }
      }
    }
  }
  void StartStep() {
    auto sem = m.GetSem("start");
    if (owner_->GetTime() == 0.) {
      if (par.extrapolate_boundaries) {
        ExtrapolateLinear(sem, fcu_.time_curr, fccl_);
      }
      if (sem.Nested("reconst")) {
        ReconstPlanes(fcu_.time_curr);
      }
    }
    if (sem("rotate")) {
      owner_->ClearIter();
      fcu_.time_prev = fcu_.time_curr;
      fcu_.iter_curr = fcu_.time_prev;
      fcus_.time_prev = fcus_.time_curr;
      fcus_.iter_curr = fcus_.time_prev;
      UpdateBc(mebc_);
    }
  }
  enum class SweepType {
    plain, // sum of fluxes
    EI, // Euler Implicit (aulisa2009)
    LE, // Lagrange Explicit (aulisa2009)
    weymouth, // sum of fluxes and divergence (weymouth2010)
  };
  void DumpInterface(
      std::string filename,
      std::vector<Multi<const FieldCell<Scal>*>> extra_fields,
      std::vector<std::string> extra_names) {
    typename UVof<M>::DumpPolyArgs args;
    args.layers = layers;
    args.fcu = fcu_.time_curr;
    args.fccl = fccl_;
    args.fcn = fcn_;
    args.fca = fca_;
    args.fci = fci_;
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
          layers, fcut, fcclt, fcn_, filename, owner_->GetTime(), par.vtkpoly,
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
        for (const auto& it : mebc_.GetMapFace()) {
          ffv[it.first] = 0;
        }
        Sweep(
            mfcu, d, layers, ffv, fccl_, fcim_, fcn_, fca_, me_vf_,
            SweepType::weymouth, nullptr, nullptr, fcuu_, 1, par.clipth, eb);
      }
      CommRec(sem, mfcu, fccl_, fcim_);
      if (par.extrapolate_boundaries) {
        ExtrapolatePlic(sem, mfcu, fccl_);
      }
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
    for (auto d : dd) {
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
  // mebc: face conditions
  // type: sweep type
  // fcfm,fcfp: upwind mixture flux, required if type=2 [s]
  // fcuu: volume fraction for Weymouth div term
  // dt: time step
  // clipth: threshold for clipping, values outside [th,1-th] are clipped
  static void Sweep(
      const Multi<FieldCell<Scal>*>& mfcu, size_t dir,
      const GRange<size_t>& layers, const FieldFace<Scal>& ffv,
      const Multi<FieldCell<Scal>*>& mfccl,
      const Multi<FieldCell<Scal>*>& mfcim,
      const Multi<const FieldCell<Vect>*>& mfcn,
      const Multi<const FieldCell<Scal>*>& mfca,
      const MapEmbed<BCond<Scal>>& mebc, SweepType type,
      const FieldCell<Scal>* fcfm, const FieldCell<Scal>* fcfp,
      const Multi<const FieldCell<Scal>*>& mfcuu, Scal dt, Scal clipth,
      const EB& eb) {
    using MIdx = typename M::MIdx;
    const auto& m = eb.GetMesh();
    const auto& indexc = m.GetIndexCells();
    const auto& indexf = m.GetIndexFaces();
    const MIdx globalsize = m.GetGlobalSize();
    const auto h = m.GetCellSize();
    const auto d = m.direction(dir);

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
        auto p = indexf.GetMIdxDir(f);
        if (p.second.raw() != d) {
          continue;
        }

        // flux through face (maybe cut)
        const Scal v = ffv[f];
        // flux through full face that would give the same velocity
        const Scal v0 = v / eb.GetAreaFraction(f);
        const IdxCell c = m.GetCell(f, v > 0 ? 0 : 1); // upwind cell
        if (v != 0 && fccl[c] != kClNone) {
          ffcl[f] = fccl[c];
          if (fcu[c] > 0 && fcu[c] < 1 && fcn[c].sqrnorm() > 0) {
            switch (type) {
              case SweepType::plain:
              case SweepType::EI:
              case SweepType::weymouth: {
                const Scal vu0 = R::GetLineFlux(fcn[c], fca[c], h, v0, dt, d);
                ffvu[f] = (v >= 0 ? std::min(vu0, v) : std::max(vu0, v));
                break;
              }
              case SweepType::LE: {
                const Scal vc = (v > 0 ? (*fcfm)[c] : (*fcfp)[c]);
                const Scal vc0 = vc / eb.GetAreaFraction(f);
                const Scal vu0 =
                    R::GetLineFluxStr(fcn[c], fca[c], h, v0, vc0, dt, d);
                ffvu[f] = (v >= 0 ? std::min(vu0, v) : std::max(vu0, v));
                break;
              }
            }
          } else { // pure cell or cell with undefined normal
            ffvu[f] = v * fcu[c];
          }

          // propagate to downwind cell if empty
          IdxCell cd = eb.GetCell(f, v > 0 ? 1 : 0); // downwind cell
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
                const MIdx w = indexc.GetMIdx(c);
                MIdx im = TRM::Unpack((*mfcim[i])[c]);
                if (w[d] < 0) im[d] += 1;
                if (w[d] >= globalsize[d]) im[d] -= 1;
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
          if ((bc.nci == 0) != (v > 0)) {
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
      for (auto c : eb.CellsM()) {
        if (fccl[c] != kClNone) {
          const auto fm = c.face(-d);
          const auto fp = c.face(d);
          // mixture cfl
          const Scal ds = (ffv[fp] - ffv[fm]) * dt / c.volume;
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
          const Scal dl = (vp - vm) * dt / c.volume;
          auto& u = fcu[c];
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
  // Removes orphan fragments, those for which the volume fraction
  // in the 3x3x3 stencil does not exceed the threshold.
  // fcu: volume fractions
  // fccl: colors
  // fcim: packed image vectors
  // filterth: threshold for detecting orphan fragments
  static void FilterOrphan(
      const Multi<FieldCell<Scal>*>& fcu, const GRange<size_t>& layers,
      const Multi<FieldCell<Scal>*>& fccl, const Multi<FieldCell<Scal>*>& fcim,
      Scal filterth, const EB& eb) {
    for (auto l : layers) {
      for (auto c : eb.Cells()) {
        if ((*fccl[l])[c] != kClNone) {
          bool orphan = true;
          for (auto ln : layers) {
            for (auto cn : eb.Stencil(c)) {
              if ((*fccl[l])[c] == (*fccl[ln])[cn] &&
                  (*fcu[ln])[cn] >= filterth) {
                orphan = false;
                break;
              }
            }
            if (!orphan) {
              break;
            }
          }
          if (orphan) {
            (*fcu[l])[c] = 0;
            (*fccl[l])[c] = kClNone;
            (*fcim[l])[c] = TRM::Pack(MIdx(0));
          }
        }
      }
    }
  }
  void CommRec(
      Sem& sem, const Multi<FieldCell<Scal>*>& mfcu,
      const Multi<FieldCell<Scal>*>& mfccl,
      const Multi<FieldCell<Scal>*>& mfcim) {
    if (sem("comm")) {
      for (auto i : layers) {
        m.Comm(mfcu[i]);
        m.Comm(mfccl[i]);
        m.Comm(mfcim[i]);
      }
    }
    if (sem("bcreflect")) {
      for (auto i : layers) {
        BcApply(*mfcu[i], me_vf_, m);
        BcApply(*mfccl[i], me_cl_, m);
        BcApply(*mfcim[i], me_im_, m);
      }
    }
    if (par.extrapolate_boundaries) {
      ExtrapolateLinear(sem, mfcu, mfccl);
    }
    if (sem.Nested("reconst")) {
      ReconstPlanes(mfcu);
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
    }

    using Scheme = typename Par::Scheme;
    const Multi<FieldCell<Scal>*> mfcu = fcu_.iter_curr;
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
    if (par.filterth > 0) {
      if (sem("filterorphan")) {
        FilterOrphan(mfcu, layers, fccl_, fcim_, par.filterth, eb);
      }
      CommRec(sem, mfcu, fccl_, fcim_);
    }
    if (modifier_) {
      if (sem("modify")) {
        modifier_(mfcu, fccl_, layers, eb);
      }
      CommRec(sem, mfcu, fccl_, fcim_);
    }
    if (sem("bcc_clear")) {
      if (par.cloverride) {
        for (auto i : layers) {
          UVof<M>::BcClearOverrideColor(
              fcu_.iter_curr[i], fccl_[i], 0., mebc_, m);
        }
      } else {
        for (auto i : layers) {
          UVof<M>::BcClear(fcu_.iter_curr[i], mebc_, m);
        }
      }
    }
    if (sem.Nested() && par.recolor) {
      uvof_.Recolor(
          layers, fcu_.iter_curr, fccl_, fccl_, par.clfixed, par.clfixed_x,
          par.coalth, me_cl_, par.verb, par.recolor_unionfind,
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
    modifier_ = nullptr;
  }
  void PostStep() {
    auto sem = m.GetSem("iter");
    // --> fcu [a], fca [s], fcn [s]
    if (sem("comm")) {
      for (auto l : layers) {
        m.Comm(&fcn_[l]);
        m.Comm(&fca_[l]);
      }
    }
    if (sem("reflect")) {
      // --> fca [a], fcn [a]
      for (auto l : layers) {
        BcApply(fcn_[l], me_n_, m);
        BcApply(fca_[l], me_a_, m);
      }
      // --> reflected fca [a], fcn [a]

      // unpack image vector
      for (auto c : m.AllCells()) {
        for (auto l : layers) {
          fcim_unpack_[l][c] = TRM::Unpack(fcim_[l][c]);
        }
      }
    }
  }
  static std::string GetFieldPath(
      std::string dirpath, std::string field, size_t layer) {
    return dirpath + "/" + field + "_" + std::to_string(layer) + ".h5";
  }
  void SaveState(std::string dirpath) const {
    auto sem = m.GetSem();
    if (sem()) {
      if (m.IsRoot()) {
        util::Makedir(dirpath);
      }
    }
    for (auto l : layers) {
      auto fieldpath = [dirpath, l](std::string field) {
        return GetFieldPath(dirpath, field, l);
      };
      if (sem.Nested()) {
        Hdf<M>::Write(fcu_.time_curr[l], fieldpath("vf"), m);
      }
      if (sem.Nested()) {
        Hdf<M>::Write(fccl_[l], fieldpath("cl"), m);
      }
    }
  }
  void LoadState(std::string dirpath) {
    auto sem = m.GetSem();
    if (sem()) {
      if (m.IsRoot()) {
        fassert(util::IsDir(dirpath), "Not a directory '" + dirpath + "'");
      }
    }
    for (auto l : layers) {
      auto fieldpath = [dirpath, l](std::string field) {
        return GetFieldPath(dirpath, field, l);
      };
      if (sem.Nested()) {
        Hdf<M>::Read(fcu_.time_curr[l], fieldpath("vf"), m);
      }
      if (sem.Nested()) {
        Hdf<M>::Read(fccl_[l], fieldpath("cl"), m);
      }
      if (sem()) {
        m.Comm(&fcu_.time_curr[l]);
        m.Comm(&fccl_[l]);
      }
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
  Multi<FieldCell<MIdx>> fcim_unpack_; // image unpacked
  size_t count_ = 0; // number of MakeIter() calls, used for splitting

  // boundary conditions
  const MapEmbed<BCondAdvection<Scal>>& mebc_; // conditions on advection
  MapEmbed<BCond<Scal>> me_vf_; // volume fraction
  MapEmbed<BCond<Scal>> me_cl_; // color
  MapEmbed<BCond<Scal>> me_im_; // image
  MapEmbed<BCond<Vect>> me_n_; // normal
  MapEmbed<BCond<Scal>> me_a_; // plane constant

  // tmp for MakeIteration, volume flux copied to cells
  FieldCell<Scal> fcfm_, fcfp_;
  UVof<M> uvof_;
  std::function<void(
      const Multi<FieldCell<Scal>*>& fcu, const Multi<FieldCell<Scal>*>& fccl,
      GRange<size_t> layers, const EB&)>
      modifier_;
};

template <class EB_>
constexpr typename EB_::M::Scal Vofm<EB_>::kClNone;

template <class EB_>
Vofm<EB_>::Vofm(
    M& m_, const EB& eb, const FieldCell<Scal>& fcu0,
    const FieldCell<Scal>& fccl0, const MapEmbed<BCondAdvection<Scal>>& mebc,
    const FieldEmbed<Scal>* fev, const FieldCell<Scal>* fcs, double t,
    double dt, Par par)
    : AdvectionSolver<M>(t, dt, m_, fev, fcs, mebc) {
  const GRange<size_t> layers(par.layers);
  Multi<FieldCell<Scal>> fcu(layers, m, 0);
  Multi<FieldCell<Scal>> fccl(layers, m, kClNone);
  fcu[0] = fcu0;
  fccl[0] = fccl0;
  imp.reset(new Imp(this, eb, layers, fcu, fccl, par));
}

template <class EB_>
Vofm<EB_>::Vofm(
    M& m_, const EB& eb, const Multi<const FieldCell<Scal>*>& fcu0,
    const Multi<const FieldCell<Scal>*>& fccl0,
    const MapEmbed<BCondAdvection<Scal>>& mebc, const FieldEmbed<Scal>* fev,
    const FieldCell<Scal>* fcs, double t, double dt, Par par)
    : AdvectionSolver<M>(t, dt, m_, fev, fcs, mebc) {
  const GRange<size_t> layers(par.layers);
  imp.reset(new Imp(this, eb, layers, fcu0, fccl0, par));
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
void Vofm<EB_>::Sharpen() {
  imp->Sharpen(imp->fcu_.iter_curr);
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
auto Vofm<EB_>::GetPlic() const -> Plic {
  return {imp->layers, GetFieldM(), GetAlpha(), GetNormal(), GetMask(),
          GetColor(),  GetColor(),  GetImage(), imp->mebc_};
}

template <class EB_>
size_t Vofm<EB_>::GetNumLayers() const {
  return imp->layers.size();
}

template <class EB_>
void Vofm<EB_>::SaveState(std::string dirpath) const {
  return imp->SaveState(dirpath);
}

template <class EB_>
void Vofm<EB_>::LoadState(std::string dirpath) {
  return imp->LoadState(dirpath);
}

template <class EB_>
auto Vofm<EB_>::GetNormal() const -> Multi<const FieldCell<Vect>*> {
  return imp->fcn_;
}

template <class EB_>
auto Vofm<EB_>::GetImage() const -> Multi<const FieldCell<MIdx>*> {
  return imp->fcim_unpack_;
}

template <class EB_>
void Vofm<EB_>::PostStep() {
  return imp->PostStep();
}

template <class EB_>
void Vofm<EB_>::DumpInterface(
    std::string filename,
    std::vector<Multi<const FieldCell<Scal>*>> extra_fields,
    std::vector<std::string> extra_names) const {
  return imp->DumpInterface(filename, extra_fields, extra_names);
}

template <class EB_>
void Vofm<EB_>::DumpInterfaceMarch(std::string filename) const {
  return imp->DumpInterfaceMarch(filename);
}

template <class EB_>
void Vofm<EB_>::AddModifier(
    std::function<void(
        const Multi<FieldCell<Scal>*>&, const Multi<FieldCell<Scal>*>&,
        GRange<size_t>, const EB&)>
        func) {
  imp->modifier_ = func;
}
