#pragma once

#include <cmath>
#include <sstream>
#include <stdexcept>

#include "util/metrics.h"
#include "util/convdiff.h"

namespace solver {

template <class M_, class CD_>
struct ConvDiffVectImp<M_, CD_>::Imp {
  using Owner = ConvDiffVectImp<M_, CD_>;

  Imp(
      Owner* owner, const FieldCell<Vect>& fcvel, 
      const MapFace<std::shared_ptr<CondFace>>& mfc, 
      const MapCell<std::shared_ptr<CondCell>>& mcc, 
      std::shared_ptr<Par> par)
      : owner_(owner), par(par), m(owner_->m)
      , mfc_(mfc), mcc_(mcc), dr_(0, m.GetEdim())
  {
    for (auto d : dr_) {
      // Face conditions for each velocity component
      vmfc_[d] = GetScalarCond(mfc_, d, m);

      // Cell conditions for each velocity component
      for (auto it : mcc) {
        IdxCell c = it.GetIdx();
        CondCell* cb = it.GetValue().get(); // cond base

        if (auto cd = dynamic_cast<CondCellVal<Vect>*>(cb)) {
          // TODO: revise with CondCellValComp
          vmcc_[d][c] = std::make_shared<
              CondCellValFixed<Scal>>(cd->GetValue()[d]);
        } else {
          throw std::runtime_error("convdiffvi: unknown cell condition");
        }
      }

      // Components of source
      for (auto d : dr_) {
        vfcs_[d] = GetComponent(*owner_->fcs_, d);
      }

      // Initialize solver
      vs_[d] = std::make_shared<CD>(
          m, GetComponent(fcvel, d), vmfc_[d], vmcc_[d],
          owner_->fcr_, owner_->ffd_, &(vfcs_[d]), owner_->ffv_, 
          owner_->GetTime(), owner_->GetTimeStep(), par);
    }
    CopyToVect(Layers::time_curr, fcvel_);
    lvel_ = Layers::time_curr;
  }
  // Copy from solvers to vector field.
  // l: layer in component solvers
  // fcv: vector field
  void CopyToVect(Layers l, FieldCell<Vect>& fcv) {
    fcv.Reinit(m);
    for (auto d : dr_) {
      SetComponent(fcv, d, vs_[d]->GetField(l));
    }
  }
  void StartStep() {
    auto sem = m.GetSem("convdiffmulti-start");
    for (auto d : dr_) {
      if (sem("dir-init")) {
        vs_[d]->SetTimeStep(owner_->GetTimeStep());
      }
      if (sem.Nested("dir-start")) {
        vs_[d]->StartStep();
      }
    }
    if (sem("tovect")) {
      lvel_ = Layers::iter_curr;
      owner_->ClearIter();
    }
  }
  void Assemble(const FieldCell<Vect>& fcw, const FieldFace<Scal>& ffv) {
    auto sem = m.GetSem("asm");
    if (sem("copy")) {
      for (auto d : dr_) {
        vfcs_[d] = GetComponent(*owner_->fcs_, d);
        vfct_[d] = GetComponent(fcw, d);
      }
    }
    for (auto d : dr_) {
      if (sem.Nested("one")) {
        vs_[d]->Assemble(vfct_[d], ffv);
      }
    }
  }
  void MakeIteration() {
    auto sem = m.GetSem("convdiffmulti-iter");
    if (sem("source")) {
      for (auto d : dr_) {
        vfcs_[d] = GetComponent(*owner_->fcs_, d);
      }
    }

    for (auto d : dr_) {
      if (sem.Nested("dir-iter")) {
        vs_[d]->MakeIteration();
      }
      if (sem("dir-linreport")) {
        if (m.IsRoot() && par->linreport) {
          std::cout 
              << "v" << ("xyz"[d]) << ":" 
              << " res=" << m.GetResidual()
              << " iter=" << m.GetIter()
              << std::endl;
        }
      }
    }

    if (sem("tovect")) {
      CopyToVect(Layers::iter_curr, fcvel_);
      lvel_ = Layers::iter_curr;
      owner_->IncIter();
    }
  }
  void FinishStep() {
    auto sem = m.GetSem("convdiffmulti-finish");

    for (auto d : dr_) {
      if (sem.Nested("dir-finish")) {
        vs_[d]->FinishStep();
      }
    }
    if (sem("tovect")) {
      lvel_ = Layers::time_curr;
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
  const FieldCell<Vect>& GetVelocity(Layers l) const {
    if (l == lvel_) {
      return fcvel_;
    }
    throw std::runtime_error(
        "GetVelocity: requested layer '" + 
        GetName(l) + "' but '" + 
        GetName(lvel_) + "' is loaded");
  }
  void CorrectVelocity(Layers l, const FieldCell<Vect>& fc) {
    auto sem = m.GetSem("corr");
    for (auto d : dr_) {
      if (sem.Nested("dir-corr")) {
        vs_[d]->CorrectField(l, GetComponent(fc, d));
      }
    }
    if (sem("tovect")) {
      CopyToVect(l, fcvel_);
      lvel_ = l;
    }
  }

  Owner* owner_;
  std::shared_ptr<Par> par;
  M& m; // mesh

  FieldCell<Vect> fcvel_;
  Layers lvel_; // current level loaded in fcvel_
  MapFace<std::shared_ptr<CondFace>> mfc_; // vect face cond
  MapCell<std::shared_ptr<CondCell>> mcc_; // vect cell cond
  GRange<size_t> dr_;  // effective dimension range

  template <class T>
  using VectGeneric = std::array<T, dim>;

  // Scalar components
  VectGeneric<MapFace<std::shared_ptr<CondFace>>> vmfc_; // face cond
  VectGeneric<MapCell<std::shared_ptr<CondCell>>> vmcc_; // cell cond
  VectGeneric<std::shared_ptr<CD>> vs_; // solver
  VectGeneric<FieldCell<Scal>> vfcs_; // force
  VectGeneric<FieldCell<Scal>> vfct_; // tmp
};


template <class M_, class CD_>
ConvDiffVectImp<M_, CD_>::ConvDiffVectImp(
    M& m, const FieldCell<Vect>& fcvel, 
    const MapFace<std::shared_ptr<CondFace>>& mfc, 
    const MapCell<std::shared_ptr<CondCell>>& mcc, 
    const FieldCell<Scal>* fcr, const FieldFace<Scal>* ffd,
    const FieldCell<Vect>* fcs, const FieldFace<Scal>* ffv,
    double t, double dt, std::shared_ptr<Par> par)
    : ConvDiffVect<M>(t, dt, m, fcr, ffd, fcs, ffv)
    , imp(new Imp(this, fcvel, mfc, mcc, par))
{}

template <class M_, class CD_>
ConvDiffVectImp<M_, CD_>::~ConvDiffVectImp() = default;

template <class M_, class CD_>
auto ConvDiffVectImp<M_, CD_>::GetPar() -> Par* {
  return imp->par.get();
}

template <class M_, class CD_>
void ConvDiffVectImp<M_, CD_>::Assemble(const FieldCell<Vect>& fcw, 
                                        const FieldFace<Scal>& ffv) {
  imp->Assemble(fcw, ffv);
}

template <class M_, class CD_>
void ConvDiffVectImp<M_, CD_>::CorrectVelocity(
    Layers l, const FieldCell<Vect>& u) {
  imp->CorrectVelocity(l, u);
}

template <class M_, class CD_>
auto ConvDiffVectImp<M_, CD_>::GetVelocityEquations(size_t d) const 
    -> const FieldCell<Expr>& {
  return imp->vs_[d]->GetEquations();
}

template <class M_, class CD_>
auto ConvDiffVectImp<M_, CD_>::GetVelocityCond(size_t d) 
    -> MapFace<std::shared_ptr<CondFace>>& {
  return imp->vmfc_[d];
}

template <class M_, class CD_>
void ConvDiffVectImp<M_, CD_>::StartStep() {
  imp->StartStep();
}

template <class M_, class CD_>
void ConvDiffVectImp<M_, CD_>::MakeIteration() {
  imp->MakeIteration();
}

template <class M_, class CD_>
void ConvDiffVectImp<M_, CD_>::FinishStep() {
  imp->FinishStep();
}

template <class M_, class CD_>
double ConvDiffVectImp<M_, CD_>::GetError() const {
  return imp->GetError();
}

template <class M_, class CD_>
auto ConvDiffVectImp<M_, CD_>::GetVelocity(Layers l) const 
    -> const FieldCell<Vect>& {
  return imp->GetVelocity(l);
}

} // namespace solver

