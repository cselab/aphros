// Created by Petr Karnakov on 23.01.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <memory>

#include "advection.h"
#include "embed.h"
#include "partstrmeshm.h"
#include "vof.h"

template <class M_>
class VofEmbed final : public AdvectionSolver<M_> {
 public:
  using M = M_;
  using P = AdvectionSolver<M>;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;
  using Par = typename Vof<M>::Par;

  // Constructor
  VofEmbed(
      M& m, const Embed<M>& eb, const FieldCell<Scal>& fcu,
      const FieldCell<Scal>& fccl, const MapCondFaceAdvection<Scal>& mfc,
      const FieldFace<Scal>* ffv, const FieldCell<Scal>* fcs, double t,
      double dt, Par par);
  ~VofEmbed();
  const Embed<M>& GetEmbed() const;
  const Par& GetPar() const;
  void SetPar(Par);
  void StartStep() override;
  void MakeIteration() override;
  void FinishStep() override;
  void PostStep() override;
  // Volume fraction
  const FieldCell<Scal>& GetField(Step l) const override;
  using P::GetField;
  // Interface mask
  const FieldCell<bool>& GetMask() const;
  // Plane constant
  const FieldCell<Scal>& GetAlpha() const;
  // Normal to interface
  const FieldCell<Vect>& GetNormal() const;
  // Color
  const FieldCell<Scal>& GetColor() const;
  static constexpr Scal kClNone = -1; // no color
  // Image vector, number of passes through periodic boundaries
  MIdx GetImage(IdxCell c) const;
  void DumpInterface(std::string filename) const override;
  void DumpInterfaceMarch(std::string filename) const override;

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};
