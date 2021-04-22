// Created by Petr Karnakov on 20.11.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <cmath>
#include <sstream>
#include <stdexcept>

#include "convdiffvg.h"
#include "solver/approx_eb.h"
#include "util/convdiff.h"
#include "util/metrics.h"

template <class M_, class CD_>
struct ConvDiffVectGeneric<M_, CD_>::Imp {
  using Owner = ConvDiffVectGeneric<M_, CD_>;

  Imp(Owner* owner, const Args& args)
      : owner_(owner)
      , par(owner_->GetPar())
      , m(owner_->m)
      , eb(owner_->eb)
      , mebc_(args.mebc)
      , dr_(0, m.GetEdim()) {
    for (auto d : dr_) {
      UpdateDerivedCond(d);

      // source term
      vfcs_[d] = GetComponent(*owner_->fcs_, d);

      // initial velocity
      auto fcu = GetComponent(args.fcvel, d);
      fcu.SetName(std::string("velocity_") + "xyz"[d]);

      // solver
      const ConvDiffArgs<EB> cdargs{
          fcu,
          vmebc_[d],
          owner_->fcr_,
          owner_->ffd_,
          &(vfcs_[d]),
          owner_->ffv_,
          owner_->GetTime(),
          owner_->GetTimeStep(),
          args.linsolver,
          par};
      vs_[d] = std::make_shared<CD>(m, eb, cdargs);
    }
    CopyToVect(Step::time_curr, fcvel_);
    lvel_ = Step::time_curr;
  }
  void UpdateDerivedCond(size_t d) {
    // Face conditions for each velocity component
    vmebc_[d] = GetScalarCond(mebc_, d, m);
  }
  // Copy from solvers to vector field.
  // l: layer in component solvers
  // fcv: vector field
  void CopyToVect(Step l, FieldCell<Vect>& fcv) {
    fcv.Reinit(m);
    for (auto d : dr_) {
      SetComponent(fcv, d, vs_[d]->GetField(l));
    }
  }
  void StartStep() {
    auto sem = m.GetSem("convdiffv-start");
    for (auto d : dr_) {
      if (sem("scal-init")) {
        vs_[d]->SetTimeStep(owner_->GetTimeStep());
        vs_[d]->SetPar(par);
      }
      if (sem.Nested("scal-start")) {
        vs_[d]->StartStep();
      }
    }
    if (sem("tovect")) {
      lvel_ = Step::iter_curr;
      owner_->ClearIter();
    }
  }
  void Assemble(const FieldCell<Vect>& fcw, const FieldFaceb<Scal>& ffv) {
    auto sem = m.GetSem("asm");
    if (sem("copy")) {
      for (auto d : dr_) {
        vfcs_[d] = GetComponent(*owner_->fcs_, d);
        vfct_[d] = GetComponent(fcw, d);
        vfct_[d].SetName(std::string("velocity_") + "xyz"[d]);
      }
    }
    for (auto d : dr_) {
      if (sem.Nested("one")) {
        vs_[d]->Assemble(vfct_[d], ffv);
      }
    }
  }
  void MakeIteration() {
    auto sem = m.GetSem("convdiffv-iter");
    if (sem("bc-source")) {
      for (auto d : dr_) {
        UpdateDerivedCond(d);
        vfcs_[d] = GetComponent(*owner_->fcs_, d);
      }
    }

    for (auto d : dr_) {
      if (sem.Nested("scal-iter")) {
        vs_[d]->MakeIteration();
      }
    }

    if (sem("tovect")) {
      CopyToVect(Step::iter_curr, fcvel_);
      lvel_ = Step::iter_curr;
      owner_->IncIter();
    }
  }
  void FinishStep() {
    auto sem = m.GetSem("convdiffv-finish");

    for (auto d : dr_) {
      if (sem.Nested("scal-finish")) {
        vs_[d]->FinishStep();
      }
    }
    if (sem("tovect")) {
      lvel_ = Step::time_curr;
      owner_->IncTime();
    }
  }
  double GetError() const {
    Scal r = 0.; // result
    for (auto d : dr_) {
      r = std::max<Scal>(r, vs_[d]->GetError());
    }
    return r;
  }
  const FieldCell<Vect>& GetVelocity(Step l) const {
    if (l == lvel_) {
      return fcvel_;
    }
    fassert(
        false, "GetVelocity: requested layer '" + GetName(l) + "' but '" +
                   GetName(lvel_) + "' is loaded");
  }
  void CorrectVelocity(Step l, const FieldCell<Vect>& fc) {
    auto sem = m.GetSem("corr");
    for (auto d : dr_) {
      if (sem.Nested("scal-corr")) {
        vs_[d]->CorrectField(l, GetComponent(fc, d));
      }
    }
    if (sem("tovect")) {
      CopyToVect(l, fcvel_);
      lvel_ = l;
    }
  }

  Owner* owner_;
  const Par& par;
  M& m; // mesh
  const EB& eb;

  FieldCell<Vect> fcvel_;
  Step lvel_; // current level loaded in fcvel_
  const MapEmbed<BCond<Vect>>& mebc_; // vect face cond
  GRange<size_t> dr_; // effective dimension range

  template <class T>
  using Array = std::array<T, dim>;

  // Scalar components
  Array<MapEmbed<BCond<Scal>>> vmebc_; // boundary conditions
  Array<std::shared_ptr<CD>> vs_; // solver
  Array<FieldCell<Scal>> vfcs_; // force
  Array<FieldCell<Scal>> vfct_; // tmp
};

template <class M_, class CD_>
ConvDiffVectGeneric<M_, CD_>::ConvDiffVectGeneric(
    M& m_, const EB& eb_, const Args& args)
    : Base(
          args.t, args.dt, m_, eb_, args.par, args.fcr, args.ffd, args.fcs,
          args.ffv)
    , imp(new Imp(this, args)) {}

template <class M_, class CD_>
ConvDiffVectGeneric<M_, CD_>::~ConvDiffVectGeneric() = default;

template <class M_, class CD_>
void ConvDiffVectGeneric<M_, CD_>::Assemble(
    const FieldCell<Vect>& fcw, const FieldFaceb<Scal>& ffv) {
  imp->Assemble(fcw, ffv);
}

template <class M_, class CD_>
void ConvDiffVectGeneric<M_, CD_>::CorrectVelocity(
    Step l, const FieldCell<Vect>& u) {
  imp->CorrectVelocity(l, u);
}

template <class M_, class CD_>
auto ConvDiffVectGeneric<M_, CD_>::GetDiag(size_t d) const -> FieldCell<Scal> {
  return imp->vs_[d]->GetDiag();
}

template <class M_, class CD_>
auto ConvDiffVectGeneric<M_, CD_>::GetConst(size_t d) const -> FieldCell<Scal> {
  return imp->vs_[d]->GetConst();
}

template <class M_, class CD_>
auto ConvDiffVectGeneric<M_, CD_>::GetVelocityCond(size_t d)
    -> MapEmbed<BCond<Scal>>& {
  return imp->vmebc_[d];
}

template <class M_, class CD_>
void ConvDiffVectGeneric<M_, CD_>::StartStep() {
  imp->StartStep();
}

template <class M_, class CD_>
void ConvDiffVectGeneric<M_, CD_>::MakeIteration() {
  imp->MakeIteration();
}

template <class M_, class CD_>
void ConvDiffVectGeneric<M_, CD_>::FinishStep() {
  imp->FinishStep();
}

template <class M_, class CD_>
double ConvDiffVectGeneric<M_, CD_>::GetError() const {
  return imp->GetError();
}

template <class M_, class CD_>
auto ConvDiffVectGeneric<M_, CD_>::GetVelocity(Step l) const
    -> const FieldCell<Vect>& {
  return imp->GetVelocity(l);
}
