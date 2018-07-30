#pragma once

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <memory>

#include "util/metrics.h"
#include "convdiffv.h"
#include "convdiffi.h"

namespace solver {

template <class M_>
class ConvectionDiffusionImplicit : public ConvectionDiffusion<M_> {
  using M = M_;
  using P = ConvectionDiffusion<M_>;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using CD = ConvDiffScalImp<M>;

  static constexpr size_t dim = M::dim;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;
  template <class T>
  using VectGeneric = std::array<T, dim>;

  M& m;
  using P::fcr_;
  using P::ffd_;
  using P::fcs_;
  using P::ffv_;
  FieldCell<Vect> fcvel_;
  Layers lvel_; // current level loaded in fcvel_
  MapFace<std::shared_ptr<CondFace>> mfc_; // vect face cond
  MapCell<std::shared_ptr<CondCell>> mcc_; // vect cell cond
  GRange<size_t> dr_;  // dimension range

  // Scalar components
  VectGeneric<MapFace<std::shared_ptr<CondFace>>> vmfc_; // face cond
  // TODO *** Cell cond
  VectGeneric<std::shared_ptr<CD>> vs_; // solver
  VectGeneric<FieldCell<Scal>> vfcs_; // force
  VectGeneric<FieldCell<Scal>> vfct_; // tmp

 public:
  using Par = typename CD::Par;
  std::shared_ptr<Par> par;
  Par* GetPar() { return par.get(); }
  // Copy from solvers to vector field.
  // l: layer in component solvers
  // fcv: vector field
  void CopyToVect(Layers l, FieldCell<Vect>& fcv) {
    fcv.Reinit(m);
    for (auto d : dr_) {
      SetComponent(fcv, d, vs_[d]->GetField(l));
    }
  }
  ConvectionDiffusionImplicit(
      M& m, // mesh
      const FieldCell<Vect>& fcvel, // initial velocity
      const MapFace<std::shared_ptr<CondFace>>& mfc, // face conditions
      const MapCell<std::shared_ptr<CondCell>>& mcc, // cell conditions
      const FieldCell<Scal>* fcr, // density
      const FieldFace<Scal>* ffd, // dynamic viscosity
      const FieldCell<Vect>* fcs, // force
      const FieldFace<Scal>* ffv, // volume flux
      double t, double dt, std::shared_ptr<Par> par)
      : ConvectionDiffusion<M>(t, dt, fcr, ffd, fcs, ffv)
      , m(m) , mfc_(mfc) , mcc_(mcc), dr_(0, dim), par(par)
  {
    for (auto d : dr_) {
      // Face conditions for each velocity component
      // (copied from given vector conditions)
      for (auto it = mfc_.cbegin(); it != mfc_.cend(); ++it) {
        IdxFace f = it->GetIdx();
        if (auto p = dynamic_cast<CondFaceVal<Vect>*>(mfc_[f].get())) {
          vmfc_[d][f] = std::make_shared<CondFaceValComp<Vect>>(p, d);
        } else {
          throw std::runtime_error("Unknown face condition type");
        }
      }

      // Initialize solver
      vs_[d] = std::make_shared<CD>(
          m, GetComponent(fcvel, d), vmfc_[d],
          MapCell<std::shared_ptr<CondCell>>(), /*TODO *** Cell cond */
          fcr_, ffd_, &(vfcs_[d]), ffv, t, dt, par);
    }
    CopyToVect(Layers::time_curr, fcvel_);
    lvel_ = Layers::time_curr;
  }
  void StartStep() override {
    auto sem = m.GetSem("convdiffmulti-start");
    for (auto d : dr_) {
      if (sem("dir-init")) {
        vs_[d]->SetTimeStep(this->GetTimeStep());
      }
      if (sem.Nested("dir-start")) {
        vs_[d]->StartStep();
      }
    }
    if (sem("tovect")) {
      lvel_ = Layers::iter_curr;
      this->ClearIter();
    }
  }
  void Assemble(const FieldCell<Vect>& fcw,
                const FieldFace<Scal>& ffv) {
    auto sem = m.GetSem("asm");
    if (sem("copy")) {
      for (auto d : dr_) {
        vfcs_[d] = GetComponent(*fcs_, d);
        vfct_[d] = GetComponent(fcw, d);
      }
    }
    for (auto d : dr_) {
      if (sem.Nested("one")) {
        vs_[d]->Assemble(vfct_[d], ffv);
      }
    }
  }
  void MakeIteration() override {
    auto sem = m.GetSem("convdiffmulti-iter");
    if (sem("source")) {
      for (auto d : dr_) {
        vfcs_[d] = GetComponent(*fcs_, d);
      }
    }

    for (auto d : dr_) {
      if (sem.Nested("dir-iter")) {
        vs_[d]->MakeIteration();
      }
    }

    if (sem("tovect")) {
      CopyToVect(Layers::iter_curr, fcvel_);
      lvel_ = Layers::iter_curr;
      this->IncIter();
    }
  }
  void FinishStep() override {
    auto sem = m.GetSem("convdiffmulti-finish");

    for (auto d : dr_) {
      if (sem.Nested("dir-finish")) {
        vs_[d]->FinishStep();
      }
    }
    if (sem("tovect")) {
      lvel_ = Layers::time_curr;
      this->IncTime();
    }
  }
  double GetError() const override {
    Scal r = 0.; // result
    for (auto d : dr_) {
      r = std::max<Scal>(r, vs_[d]->GetError());
    }
    return r;
  }
  const FieldCell<Vect>& GetVelocity(Layers l) const override {
    if (l == lvel_) {
      return fcvel_;
    }
    throw std::runtime_error(
        "GetVelocity: requested layer '" + 
        GetName(l) + "' but '" + 
        GetName(lvel_) + "' is loaded");
  }
  const FieldCell<Vect>& GetVelocity() const override {
    return GetVelocity(Layers::time_curr);
  }
  void CorrectVelocity(Layers l, const FieldCell<Vect>& fc) override {
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
  const FieldCell<Expr>& GetVelocityEquations(size_t d) const override {
    return vs_[d]->GetEquations();
  }
  MapFace<std::shared_ptr<CondFace>>& GetVelocityCond(size_t d) {
    return vmfc_[d];
  }
  CD& GetSolver(size_t d) {
    return *vs_[d];
  }
};

} // namespace solver

