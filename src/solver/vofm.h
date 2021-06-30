// Created by Petr Karnakov on 01.09.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <memory>

#include "advection.h"
#include "vof.h"

template <class EB_>
class Vofm final : public AdvectionSolver<typename EB_::M> {
 public:
  using EB = EB_;
  using M = typename EB::M;
  using Base = AdvectionSolver<M>;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;
  using Par = typename Vof<M>::Par;
  using Plic = generic::Plic<Vect>;

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
      M& m, const EB& eb, const FieldCell<Scal>& fcu,
      const FieldCell<Scal>& fccl, const MapEmbed<BCondAdvection<Scal>>& mfc,
      const FieldEmbed<Scal>* fev, const FieldCell<Scal>* fcs, double t,
      double dt, Par par);
  Vofm(
      M& m, const EB& eb, const Multi<const FieldCell<Scal>*>& fcu,
      const Multi<const FieldCell<Scal>*>& fccl,
      const MapEmbed<BCondAdvection<Scal>>& mfc, const FieldEmbed<Scal>* fev,
      const FieldCell<Scal>* fcs, double t, double dt, Par par);
  ~Vofm();
  const EB& GetEmbed() const;
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
  using Base::GetField;
  using Base::m;
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
  Multi<const FieldCell<MIdx>*> GetImage() const;
  // Volume fraction, plane constant, normal, color.
  Plic GetPlic() const override;
  void DumpInterface(
      std::string filename,
      std::vector<Multi<const FieldCell<Scal>*>> extra_fields,
      std::vector<std::string> extra_names) const override;
  void DumpInterfaceMarch(std::string filename) const override;
  // Saves current state to directory, created if needed.
  void SaveState(std::string dirpath) const;
  // Loads current state from directory.
  void LoadState(std::string dirpath);
  // Adds a function that modifies the fields at the next iteration
  // and after which is removed.
  void AddModifier(std::function<void(
                       const Multi<FieldCell<Scal>*>& fcu,
                       const Multi<FieldCell<Scal>*>& fccl,
                       GRange<size_t> layers, const EB&)>);
  // Applies sharpening to current field at time time_curr
  void Sharpen();

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};
