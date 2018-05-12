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
  using CD = ConvectionDiffusionScalarImplicit<M>;

  static constexpr size_t dim = M::dim;
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;
  template <class T>
  using VectGeneric = std::array<T, dim>;

  M& m;
  using P::fcr_;
  using P::ffd_;
  using P::fcs_;
  using P::ffv_;
  LayersData<FieldCell<Vect>> fcvel_;
  MapFace<std::shared_ptr<ConditionFace>> mfc_; // vect face cond
  MapCell<std::shared_ptr<ConditionCell>> mcc_; // vect cell cond
  GRange<size_t> dr_;  // dimension range

  // Scalar components
  VectGeneric<MapFace<std::shared_ptr<ConditionFace>>> vmfc_; // face cond
  // TODO *** Cell cond
  VectGeneric<std::shared_ptr<CD>> vs_; // solver
  VectGeneric<FieldCell<Scal>> vfcs_; // force

 public:
  using Par = typename CD::Par;
  std::shared_ptr<Par> par;
  Par* GetPar() { return par.get(); }
  // Copy velocity components from solvers
  void CopyToVect(Layers l) {
    fcvel_.Get(l).Reinit(m);
    for (auto d : dr_) {
      SetComponent(fcvel_.Get(l), d, vs_[d]->GetField(l));
    }
  }
  ConvectionDiffusionImplicit(
      M& m, // mesh
      const FieldCell<Vect>& fcvel, // initial velocity
      const MapFace<std::shared_ptr<ConditionFace>>& mfc, // face conditions
      const MapCell<std::shared_ptr<ConditionCell>>& mcc, // cell conditions
      const FieldCell<Scal>* fcr, // density
      const FieldFace<Scal>* ffd, // dynamic viscosity
      const FieldCell<Vect>* fcs, // force
      const FieldFace<Scal>* ffv, // volume flux
      double t, double dt, std::shared_ptr<Par> par)
      : ConvectionDiffusion<M>(t, dt, fcr, ffd, fcs, ffv)
      , m(m) , mfc_(mfc) , mcc_(mcc), par(par), dr_(0, dim)
  {
    for (auto d : dr_) {
      // Face conditions for each velocity component
      // (copied from given vector conditions)
      for (auto it = mfc_.cbegin(); it != mfc_.cend(); ++it) {
        IdxFace f = it->GetIdx();
        if (auto p = dynamic_cast<ConditionFaceValue<Vect>*>(mfc_[f].get())) {
          vmfc_[d][f] =
              std::make_shared<ConditionFaceValueExtractComponent<Vect>>(p, d);
        } else {
          throw std::runtime_error("Unknown face condition type");
        }
      }

      // Initialize solver
      vs_[d] = std::make_shared<CD>(
          m, GetComponent(fcvel, d), vmfc_[d],
          MapCell<std::shared_ptr<ConditionCell>>(), /*TODO *** Cell cond */
          fcr_, ffd_, &(vfcs_[d]), ffv, t, dt, par);
    }
    CopyToVect(Layers::time_curr);
    CopyToVect(Layers::time_prev);
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
      CopyToVect(Layers::iter_curr);
      this->ClearIter();
    }
  }
  void MakeIteration() override {
    auto sem = m.GetSem("convdiffmulti-iter");
    for (auto d : dr_) {
      if (sem("dir-get")) {
        vfcs_[d] = GetComponent(*fcs_, d);
      }
    }

    for (auto d : dr_) {
      if (sem.Nested("dir-iter")) {
        vs_[d]->MakeIteration();
      }
    }

    if (sem("tovect")) {
      CopyToVect(Layers::iter_prev);
      CopyToVect(Layers::iter_curr);
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
      CopyToVect(Layers::time_prev);
      CopyToVect(Layers::time_curr);
      this->IncTime();
    }
  }
  double GetError() const override {
    if (this->GetIter() == 0) {
      return 1.;
    }
    return CalcDiff(fcvel_.iter_curr, fcvel_.iter_prev, m);
  }
  const FieldCell<Vect>& GetVelocity() const override {
    return fcvel_.time_curr;
  }
  const FieldCell<Vect>& GetVelocity(Layers layer) const override {
    return fcvel_.Get(layer);
  }
  void CorrectVelocity(Layers l, const FieldCell<Vect>& fc) override {
    auto sem = m.GetSem("corr");
    for (auto d : dr_) {
      if (sem.Nested("dir-corr")) {
        vs_[d]->CorrectField(l, GetComponent(fc, d));
      }
    }
    if (sem("tovect")) {
      CopyToVect(l);
    }
  }
  const FieldCell<Expr>& GetVelocityEquations(size_t d) const override {
    return vs_[d]->GetEquations();
  }
  MapFace<std::shared_ptr<ConditionFace>>& GetVelocityCond(size_t d) {
    return vmfc_[d];
  }
  const CD& GetSolver(size_t d) {
    return *vs_[d];
  }
};

} // namespace solver

