// Created by Petr Karnakov on 01.09.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <memory>

#include "advection.h"
#include "partstrmeshm.h"
#include "vof.h"

template <class M_>
class Vofm final : public AdvectionSolver<M_> {
 public:
  using M = M_;
  using P = AdvectionSolver<M>;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;
  using Par = typename Vof<M>::Par;

  static constexpr Scal kClNone = -1.; // no color

  // Constructor
  // fcu: initial volume fraction
  // fccl: initial color
  // mfc: boundary conditions for volume fraction
  // ffv: pointer to mixture flux
  // fcs: poitner to volume sources
  // t,dt: initial time and timestep
  // par: parameters
  Vofm(
      M& m, const FieldCell<Scal>& fcu, const FieldCell<Scal>& fccl,
      const MapCondFaceAdvection<Scal>& mfc, const FieldFace<Scal>* ffv,
      const FieldCell<Scal>* fcs, double t, double dt, Par par);
  Vofm(
      M& m, const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<const FieldCell<Scal>*>& fccl,
      const MapCondFaceAdvection<Scal>& mfc, const FieldFace<Scal>* ffv,
      const FieldCell<Scal>* fcs, double t, double dt, Par par);
  ~Vofm();
  const Par& GetPar() const;
  void SetPar(Par);
  void StartStep() override;
  void MakeIteration() override;
  void FinishStep() override;
  void PostStep() override;
  // Volume fraction
  const FieldCell<Scal>& GetField(Step l) const override;
  const FieldCell<Scal>& GetField(Step l, size_t i) const;
  Multi<const FieldCell<Scal>*> GetFieldM() const;
  // ...
  using P::GetField;
  // Plane constant
  Multi<const FieldCell<Scal>*> GetAlpha() const;
  // Normal to interface
  Multi<const FieldCell<Vect>*> GetNormal() const;
  // Number of layers
  size_t GetNumLayers() const;
  // Interface mask
  Multi<const FieldCell<bool>*> GetMask() const;
  // Colors
  Multi<const FieldCell<Scal>*> GetColor() const;
  // Colors combined
  const FieldCell<Scal>& GetColorSum() const;
  // Image vector, number of passes through periodic boundaries
  MIdx GetImage(size_t l, IdxCell c) const;
  void DumpInterface(std::string filename) const override;
  void DumpInterfaceMarch(std::string filename) const override;

 private:
  struct Imp; // implementation
  std::unique_ptr<Imp> imp;
};
