#pragma once

#include <memory>

#include "advection.h"

namespace solver {

template <class M_>
class Tvd : public AdvectionSolver<M_> {
 public:
  using M = M_;
  using P = AdvectionSolver<M>;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

 public:
  struct Par {
    Scal sharp = 0.;
    Scal sharpo = 0.;
    Scal sharp_max = 1.;
    bool split = false;
  };
  // Constructor
  Tvd(M& m, const FieldCell<Scal>& fcu,
      const MapFace<std::shared_ptr<CondFace>>& mfc,
      const FieldFace<Scal>* ffv, const FieldCell<Scal>* fcs,
      double t, double dt, std::shared_ptr<Par> par);
  // Parameters
  Par* GetPar();
  // ...
  void StartStep() override;
  // ...
  void MakeIteration() override;
  // ...
  void FinishStep() override;
  // ...
  const FieldCell<Scal>& GetField(Layers l) const override;
  // ...
  using P::GetField;
  // ...
  const FieldCell<Scal>& GetCurv() const override;

 private:
  struct Imp; // implementation
  std::unique_ptr<Imp> imp;

  using P::m;
};

} // namespace solver
