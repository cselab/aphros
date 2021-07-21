// Created by Petr Karnakov on 27.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <memory>

#include "advection.h"

template <class M>
class Labeling;

template <class M>
struct VofPar {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  size_t dim = 3; // dimension (dim=2 assumes zero velocity in z)
  Scal clipth = 1e-10; // vf clipping threshold
  Scal filterth = 0; // orphan filtering threshold
  bool recolor = true; // run connected component labeling on every step
  bool recolor_unionfind = true; // use union-find algorithm
  bool recolor_reduce = true; // reduce set of colors to integers
  bool recolor_grid = true; // use grid heuristic
  Scal clfixed = -1; // if >= 0, value for color at point clfixed_x
  Vect clfixed_x = Vect(1e10);
  bool cloverride = false; // XXX adhoc if clear1<1, override color with 0
  bool sharpen = false;
  Scal sharpen_cfl = 0.5;
  size_t layers = 4;
  Scal avgnorm0 = 1; // original normal with sum(u)<avgnorm0
  Scal avgnorm1 = 1; // overriden normal with sum(u)>=acgnorm1
  Scal coalth = 1e10;
  int verb = 0;
  bool bcc_reflectpoly = true; // reflection for DumpPolyMarch
  Scal dumppolymarch_fill = -1; // fill cells outside
  bool vtkbin = true;
  bool vtkmerge = true;
  bool vtkpoly = true; // dump vtk polygins instead of lines
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
  Labeling<M>* labeling = nullptr; // Pointer to implementation of
                                   // connected component labeling.
                                   // Defaults to Recolor().
};

template <class EB_>
class Vof final : public AdvectionSolver<typename EB_::M> {
 public:
  using EB = EB_;
  using M = typename EB::M;
  using P = AdvectionSolver<M>;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;
  using Plic = generic::Plic<Vect>;
  using Par = VofPar<M>;

  // Constructor
  Vof(M& m, const EB& eb, const FieldCell<Scal>& fcu,
      const FieldCell<Scal>& fccl, const MapEmbed<BCondAdvection<Scal>>& mfc,
      const FieldEmbed<Scal>* fev, const FieldCell<Scal>* fcs, double t,
      double dt, Par par);
  ~Vof();
  const EB& GetEmbed() const;
  const Par& GetPar() const;
  void SetPar(VofPar<M>);
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
  Plic GetPlic() const override;
  static constexpr Scal kClNone = -1; // no color
  // Image vector, number of passes through periodic boundaries
  const FieldCell<MIdx>& GetImage() const;
  void DumpInterface(
      std::string filename,
      std::vector<Multi<const FieldCell<Scal>*>> extra_fields,
      std::vector<std::string> extra_names) const override;
  void DumpInterfaceMarch(std::string filename) const override;
  // Adds a function that modifies the fields at the next iteration
  // and after which is removed.
  void AddModifier(
      std::function<
          void(FieldCell<Scal>& fcu, FieldCell<Scal>& fccl, const EB&)>);
  // Applies sharpening to field iter_curr
  void Sharpen();

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};
