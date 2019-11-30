#pragma once

#include <memory>

#include "advection.h"
#include "dump/dumper.h"
#include "partstrmeshm.h"
#include "vof.h"

namespace solver {

template <class M_>
class Vofm final : public AdvectionSolver<M_> {
 public:
  using M = M_;
  using P = AdvectionSolver<M>;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;

  using Par = typename Vof<M>::Par;

  // Constructor
  // fcu: initial volume fraction
  // fccl: initial color
  // mfc: boundary conditions for volume fraction
  // ffv: pointer to mixture flux
  // fcs: poitner to volume sources
  // t,dt: initial time and timestep
  // par: parameters
  Vofm(M& m, const FieldCell<Scal>& fcu, const FieldCell<Scal>& fccl,
      const MapCondFaceAdvection<Scal>& mfc,
      const FieldFace<Scal>* ffv, const FieldCell<Scal>* fcs,
      double t, double dt, std::shared_ptr<Par> par);
  ~Vofm();
  // Parameters
  Par* GetPar();
  // ...
  void StartStep() override;
  // ...
  void MakeIteration() override;
  // ...
  void FinishStep() override;
  // ...
  void PostStep() override;
  // Volume fraction
  const FieldCell<Scal>& GetField(Layers l) const override;
  const FieldCell<Scal>& GetField(Layers l, size_t i) const;
  const FieldCell<Scal>& GetField(size_t i) const;
  // ...
  using P::GetField;
  // Plane constant
  const FieldCell<Scal>& GetAlpha(size_t i) const;
  // Normal to interface
  const FieldCell<Vect>& GetNormal(size_t i) const;
  // Number of layers
  size_t GetNumLayers() const;
  // Color from one layer
  const FieldCell<Scal>& GetColor(size_t i) const;
  // Color from all layers combined
  const FieldCell<Scal>& GetColor() const;
  // Image vector, number of passes through periodic boundaries
  MIdx GetImage(size_t l, IdxCell c) const;
  // Default curvature 
  const FieldCell<Scal>& GetCurv() const override;
  // Default curvature 
  const FieldCell<Scal>& GetCurv(size_t i) const;
  static constexpr Scal kClNone = -1.; // no color

 private:
  struct Imp; // implementation
  std::unique_ptr<Imp> imp;
};

} // namespace solver
