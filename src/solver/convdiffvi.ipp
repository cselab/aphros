#pragma once

#include <cmath>
#include <sstream>
#include <stdexcept>

#include "convdiffvi.h"
#include "util/metrics.h"

namespace solver {

template <class M_>
struct ConvDiffVectImp<M_>::Imp {
  using Owner = ConvDiffVectImp<M_>;

  Imp(
      Owner* owner, const FieldCell<Vect>& fcvel, 
      const MapFace<std::shared_ptr<CondFace>>& mfc, 
      const MapCell<std::shared_ptr<CondCell>>& mcc, 
      std::shared_ptr<Par> par)
      : owner_(owner), par(par), m(owner_->m)
      , mfc_(mfc), mcc_(mcc), dr_(0, dim)
  {
    for (auto d : dr_) {
      // Face conditions for each velocity component
      // (copied from given vector conditions)
      for (auto it = mfc_.cbegin(); it != mfc_.cend(); ++it) {
        IdxFace f = it->GetIdx();
        CondFace* cb = it->GetValue().get();
        if (auto p = dynamic_cast<CondFaceVal<Vect>*>(cb)) {
          vmfc_[d][f] = std::make_shared<CondFaceValComp<Vect>>(p, d);
        } else if (auto p = dynamic_cast<CondFaceGrad<Vect>*>(cb)) {
          vmfc_[d][f] = std::make_shared<CondFaceGradComp<Vect>>(p, d);
        } else if (auto p = dynamic_cast<CondFaceReflect*>(cb)) {
          auto nci = cb->GetNci();
          // XXX: adhoc for cartesian grid
          if (d == m.GetNormal(f).abs().argmax()) { 
            // normal, zero value
            vmfc_[d][f] = std::make_shared<CondFaceValFixed<Scal>>(0., nci);
          } else { 
            // tangential, zero gradient
            vmfc_[d][f] = std::make_shared<CondFaceGradFixed<Scal>>(0., nci);
          }
        } else {
          throw std::runtime_error("convdiffvi: Unknown face condition type");
        }
      }

      // Components of source
      for (auto d : dr_) {
        vfcs_[d] = GetComponent(*owner_->fcs_, d);
      }

      // Initialize solver
      vs_[d] = std::make_shared<CD>(
          m, GetComponent(fcvel, d), vmfc_[d],
          MapCell<std::shared_ptr<CondCell>>(), /*TODO *** Cell cond */
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
  GRange<size_t> dr_;  // dimension range

  template <class T>
  using VectGeneric = std::array<T, dim>;

  // Scalar components
  VectGeneric<MapFace<std::shared_ptr<CondFace>>> vmfc_; // face cond
  // TODO *** Cell cond
  VectGeneric<std::shared_ptr<CD>> vs_; // solver
  VectGeneric<FieldCell<Scal>> vfcs_; // force
  VectGeneric<FieldCell<Scal>> vfct_; // tmp
};


template <class M_>
ConvDiffVectImp<M_>::ConvDiffVectImp(
    M& m, const FieldCell<Vect>& fcvel, 
    const MapFace<std::shared_ptr<CondFace>>& mfc, 
    const MapCell<std::shared_ptr<CondCell>>& mcc, 
    const FieldCell<Scal>* fcr, const FieldFace<Scal>* ffd,
    const FieldCell<Vect>* fcs, const FieldFace<Scal>* ffv,
    double t, double dt, std::shared_ptr<Par> par)
    : ConvDiffVect<M>(t, dt, m, fcr, ffd, fcs, ffv)
    , imp(new Imp(this, fcvel, mfc, mcc, par))
{}

template <class M_>
ConvDiffVectImp<M_>::~ConvDiffVectImp() = default;

template <class M_>
auto ConvDiffVectImp<M_>::GetPar() -> Par* {
  return imp->par.get();
}

template <class M_>
void ConvDiffVectImp<M_>::Assemble(const FieldCell<Vect>& fcw, 
                                   const FieldFace<Scal>& ffv) {
  imp->Assemble(fcw, ffv);
}

template <class M_>
void ConvDiffVectImp<M_>::CorrectVelocity(Layers l, const FieldCell<Vect>& u) {
  imp->CorrectVelocity(l, u);
}

template <class M_>
auto ConvDiffVectImp<M_>::GetVelocityEquations(size_t d) const 
    -> const FieldCell<Expr>& {
  return imp->vs_[d]->GetEquations();
}

template <class M_>
auto ConvDiffVectImp<M_>::GetVelocityCond(size_t d) 
    -> MapFace<std::shared_ptr<CondFace>>& {
  return imp->vmfc_[d];
}

template <class M_>
auto ConvDiffVectImp<M_>::GetSolver(size_t d) -> CD& {
  return *imp->vs_[d];
}

template <class M_>
void ConvDiffVectImp<M_>::StartStep() {
  imp->StartStep();
}

template <class M_>
void ConvDiffVectImp<M_>::MakeIteration() {
  imp->MakeIteration();
}

template <class M_>
void ConvDiffVectImp<M_>::FinishStep() {
  imp->FinishStep();
}

template <class M_>
double ConvDiffVectImp<M_>::GetError() const {
  return imp->GetError();
}

template <class M_>
auto ConvDiffVectImp<M_>::GetVelocity(Layers l) const 
    -> const FieldCell<Vect>& {
  return imp->GetVelocity(l);
}

} // namespace solver

