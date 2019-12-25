#pragma once

#include <array>
#include <exception>
#include <fstream>
#include <limits>
#include <memory>

#include "approx.h"
#include "tvd.h"
#include "util/vof.h"

namespace solver {

template <class M_>
struct Tvd<M_>::Imp {
  static constexpr size_t dim = M::dim;

  Imp(Tvd* owner, const FieldCell<Scal>& fcu,
      const MapCondFaceAdvection<Scal>& mfc, std::shared_ptr<Par> par)
      : owner_(owner), par(par), m(owner_->m), mfc_(mfc), fck_(m, 0) {
    fcu_.time_curr = fcu;
    for (auto it : mfc_) {
      IdxFace f = it.GetIdx();
      mfvz_[f].Set<CondFaceGradFixed<Vect>>(Vect(0), it.GetValue().GetNci());
    }
  }
  void StartStep() {
    owner_->ClearIter();
    fcu_.time_prev = fcu_.time_curr;
    fcu_.iter_curr = fcu_.time_prev;
    UVof<M>::GetAdvectionFaceCond(
        m, mfc_, mfc_vf_, mfc_cl_, mfc_im_, mfc_n_, mfc_a_);
  }
  void MakeIteration() {
    auto sem = m.GetSem("iter");
    auto& curr = fcu_.iter_curr;
    if (sem("init")) {
      const Scal dt = owner_->GetTimeStep();
      auto& fcs = *owner_->fcs_;
      for (auto c : m.Cells()) {
        curr[c] = fcu_.time_prev[c] + // previous time step
                  dt * fcs[c]; // source
      }
    }
    for (size_t d = 0; d < (par->split ? dim : 1); ++d) {
      if (sem("adv")) {
        ffu_ = Interpolate(curr, mfc_vf_, m);
        fcg_ = Gradient(ffu_, m);
        ffvu_ = InterpolateSuperbee(curr, fcg_, mfc_vf_, *(owner_->ffv_), m);

        for (auto f : m.Faces()) {
          ffvu_[f] *= (*owner_->ffv_)[f];
        }

        if (par->split) {
          // Filter volume flux in one direction
          Vect vd(0);
          vd[d] = 1.;

          // Weighten by n.dot(vd)
          for (auto f : m.Faces()) {
            auto n = m.GetNormal(f);
            ffvu_[f] *= (n.dot(vd)) / n.norm1();
          }
        }

        const Scal dt = owner_->GetTimeStep();
        for (auto c : m.Cells()) {
          Scal s = 0.; // sum of fluxes
          for (auto q : m.Nci(c)) {
            IdxFace f = m.GetFace(c, q);
            s += ffvu_[f] * m.GetOutwardFactor(c, q);
          }
          curr[c] += dt / m.GetVolume(c) * (-s); // fluxes
        }
        m.Comm(&curr);
      }
    }

    // Interface sharpening
    if (std::abs(par->sharp) != 0. && sem("sharp")) {
      auto& ac = fcu_.iter_curr;

      auto& af = sharp_af_;
      auto& gc = sharp_gc_;
      auto& gf = sharp_gf_;
      af = Interpolate(ac, mfc_vf_, m);
      gc = Gradient(af, m);
      gf = Interpolate(gc, mfvz_, m);
      const Scal kc = par->sharp;
      const Scal kd = par->sharp * par->sharpo;

      for (auto f : m.Faces()) {
        const Vect g = gf[f]; // gradient on face
        const Vect n = g / (g.norm() + 1e-6); // normal to interface
        const Vect s = m.GetSurface(f); // surface vector on face
        IdxCell cm = m.GetCell(f, 0);
        IdxCell cp = m.GetCell(f, 1);
        Scal am = ac[cm];
        Scal ap = ac[cp];
        ap = std::max(0., std::min(1., ap));
        am = std::max(0., std::min(1., am));
        ffvu_[f] = kc * (ap > am ? (1. - (ap - am)) * (1. - ap) * am
                                 : -(1. - (am - ap)) * (1. - am) * ap) -
                   kd * (ap - am);
        ffvu_[f] *= std::abs(n.dot(s));
      }

      for (auto c : m.Cells()) {
        Scal s = 0.; // sum of fluxes
        for (auto q : m.Nci(c)) {
          IdxFace f = m.GetFace(c, q);
          s += ffvu_[f] * m.GetOutwardFactor(c, q);
        }

        const Scal dt = owner_->GetTimeStep();
        ac[c] += dt / m.GetVolume(c) * (-s);
      }
      m.Comm(&ac);
    }

    if (sem("curv")) {
      ffu_ = Interpolate(curr, mfc_vf_, m); // [s]
      fcg_ = Gradient(ffu_, m); // [s]
      ffg_ = Interpolate(fcg_, mfvz_, m); // [i]

      fck_.Reinit(m); // curvature [i]
      for (auto c : m.Cells()) {
        Scal s = 0.;
        for (auto q : m.Nci(c)) {
          IdxFace f = m.GetFace(c, q);
          auto& g = ffg_[f];
          // TODO: revise 1e-6
          auto n = g / (g.norm() + 1e-6); // inner normal
          s += -n.dot(m.GetOutwardSurface(c, q));
        }
        fck_[c] = s / m.GetVolume(c);
      }
      m.Comm(&fck_);
    }

    if (sem("stat")) {
      owner_->IncIter();
    }
  }
  void FinishStep() {
    fcu_.time_curr = fcu_.iter_curr;
    owner_->IncTime();
  }
  const FieldCell<Scal>& GetField(Layers l) const {
    return fcu_.Get(l);
  }

  Tvd* owner_;
  std::shared_ptr<Par> par;
  M& m;

  // TODO: revise tmp fields
  LayersData<FieldCell<Scal>> fcu_;
  const MapCondFaceAdvection<Scal>& mfc_;
  MapCondFace mfvz_; // zero-derivative bc for Vect
  MapCondFace mfc_vf_; // conditions on vf
  MapCondFace mfc_cl_; // conditions on cl
  MapCondFace mfc_im_; // conditions on cl
  MapCondFace mfc_n_; // conditions on n
  MapCondFace mfc_a_; // conditions on a

  FieldFace<Scal> ffvu_; // flux: volume flux * field
  FieldFace<Scal> ffu_; // field on faces
  FieldCell<Vect> fcg_; // gradient
  FieldFace<Vect> ffg_; // gradient
  FieldFace<Scal> sharp_af_;
  FieldCell<Vect> sharp_gc_;
  FieldFace<Vect> sharp_gf_;
  FieldCell<Scal> fck_; // curvature
};

template <class M_>
Tvd<M_>::Tvd(
    M& m, const FieldCell<Scal>& fcu, const MapCondFaceAdvection<Scal>& mfc,
    const FieldFace<Scal>* ffv, const FieldCell<Scal>* fcs, double t, double dt,
    std::shared_ptr<Par> par)
    : AdvectionSolver<M>(t, dt, m, ffv, fcs)
    , imp(new Imp(this, fcu, mfc, par)) {}

template <class M_>
Tvd<M_>::~Tvd() = default;

template <class M_>
auto Tvd<M_>::GetPar() -> Par* {
  return imp->par.get();
}

template <class M_>
void Tvd<M_>::StartStep() {
  imp->StartStep();
}

template <class M_>
void Tvd<M_>::MakeIteration() {
  imp->MakeIteration();
}

template <class M_>
void Tvd<M_>::FinishStep() {
  imp->FinishStep();
}

template <class M_>
auto Tvd<M_>::GetField(Layers l) const -> const FieldCell<Scal>& {
  return imp->fcu_.Get(l);
}

} // namespace solver
