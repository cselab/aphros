// Created by Petr Karnakov on 23.01.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <memory>

#include "advection.h"
#include "partstrmeshm.h"

template <class M_>
class VofEmbed final : public AdvectionSolver<M_> {
 public:
  using M = M_;
  using P = AdvectionSolver<M>;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;

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
  };

  // Constructor
  VofEmbed(M& m, const FieldCell<Scal>& fcu, const FieldCell<Scal>& fccl,
      const MapCondFaceAdvection<Scal>& mfc, const FieldFace<Scal>* ffv,
      const FieldCell<Scal>* fcs, double t, double dt, Par par);
  ~VofEmbed();
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
  struct Imp; // implementation
  std::unique_ptr<Imp> imp;
};
