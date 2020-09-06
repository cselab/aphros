// Created by Petr Karnakov on 27.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <memory>

#include "advection.h"
#include "multi.h"

namespace generic {

// Piecewise Linear Interface Characterization.
// Describes the multilayer interface.
template <class Scal>
struct Plic {
  using Vect = generic::Vect<Scal, 3>;
  GRange<size_t> layers;
  Multi<const FieldCell<Scal>*> vfcu; // volume fraction
  Multi<const FieldCell<Scal>*> vfca; // plane constant
  Multi<const FieldCell<Vect>*> vfcn; // normal
  Multi<const FieldCell<bool>*> vfci; // interface mask
  Multi<const FieldCell<Scal>*> vfccl; // color
  const MapEmbed<BCondAdvection<Scal>>& me_adv; // boundary conditions
};

} // namespace generic

template <class EB_>
class Vof final : public AdvectionSolver<typename EB_::M> {
 public:
  using EB = EB_;
  using M = typename EB::M;
  using P = AdvectionSolver<M>;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;
  using Plic = generic::Plic<Scal>;

  struct Par {
    size_t dim = 3; // dimension (dim=2 assumes zero velocity in z)
    Scal clipth = 1e-6; // vf clipping threshold
    bool recolor_unionfind = true;
    bool recolor_reduce = true;
    bool recolor_grid = true;
    Scal clfixed = -1; // if >= 0, value for color at point clfixed_x
    Vect clfixed_x = Vect(1e10);
    bool cloverride = false; // XXX adhoc if clear1<1, override color with 0
    bool sharpen = false;
    Scal sharpen_cfl = 0.5;
    size_t layers = 4;
    Scal avgnorm0 = 1; // original normal with sum(u)<avgnorm0
    Scal avgnorm1 = 1; // overriden normal with sum(u)>=acgnorm1
    Scal coalth = 1.5;
    int verb = 0;
    bool bcc_reflectpoly = true; // reflection for DumpPolyMarch
    Scal dumppolymarch_fill = -1; // fill cells outside
    bool vtkbin = false;
    bool vtkmerge = true;
    Scal vtkiso = 0.5;
    enum class Scheme { plain, aulisa, weymouth };
    Scheme scheme = Scheme::weymouth;
    // Enables extrapolation to halo or cut cells,
    // required for the contact angle model.
    // Extrapolate volume fraction to excluded cells before computing normals
    // Extrapolate volume fraction, plane constant and normals to cut cells
    // after the advection step.
    // If embedded boundaries are enabled, supports only periodic condtitions.
    bool extrapolate_boundaries = false;
  };

  // Constructor
  Vof(M& m, const EB& eb, const FieldCell<Scal>& fcu,
      const FieldCell<Scal>& fccl, const MapEmbed<BCondAdvection<Scal>>& mfc,
      const FieldEmbed<Scal>* fev, const FieldCell<Scal>* fcs, double t,
      double dt, Par par);
  ~Vof();
  const EB& GetEmbed() const;
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
  // Volume fraction, plane constant, normal, color.
  Plic GetPlic() const;
  static constexpr Scal kClNone = -1; // no color
  // Image vector, number of passes through periodic boundaries
  MIdx GetImage(IdxCell c) const;
  void DumpInterface(std::string filename) const override;
  void DumpInterfaceMarch(std::string filename) const override;
  // Adds a function that modifies the fields at the next iteration
  // and after which is removed.
  void AddModifier(
      std::function<
          void(FieldCell<Scal>& fcu, FieldCell<Scal>& fccl, const EB&)>);

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};
