// Created by Petr Karnakov on 28.08.2019
// Copyright 2019 ETH Zurich

#pragma once

#include <memory>

#include "cond.h"
#include "dump/dumper.h"
#include "geom/mesh.h"
#include "geom/range.h"
#include "multi.h"
#include "partstr.h"
#include "vof.h"

template <class M_>
class PartStrMeshM {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using PS = PartStr<Scal>;
  using Plic = generic::Plic<Vect>;

  struct Par {
    typename PartStr<Scal>::Par ps;
    size_t dump_fr = -1; // maximum number of iteration steps to dump
    size_t ns = 3; // number of strings per cell
    Scal tol = 1e-5; // tolerance
    size_t itermax = 20;
    int verb = 0; // debug output
    size_t dim = 3;
    bool vtkbin = true; // write binary vtk in DumpPartInter
    bool vtkmerge = true; // merge close points in DumpPartInter
  };

  PartStrMeshM(M& m, Par par, const GRange<size_t>& layers);
  ~PartStrMeshM();

  // Computes curvature with particles.
  // vfca: plane constant
  // vfcn: normal
  // vfci: interface mask (1: contains interface)
  // vfccl: color
  template <class EB>
  void Part(const Plic& plic, const EB& eb);
  // Dump particles to csv.
  // vfca: plane constant
  // vfcn: normal
  // n: frame index
  // t: time
  void DumpParticles(
      const Multi<const FieldCell<Scal>*>& vfca,
      const Multi<const FieldCell<Vect>*>& vfcn, size_t id, Scal t) const;
  void DumpPartInter(
      const Multi<const FieldCell<Scal>*>& vfca,
      const Multi<const FieldCell<Vect>*>& vfcn, size_t id, Scal t) const;
  void DumpParticles(
      const FieldCell<Scal>& fca, const FieldCell<Vect>& fcn, size_t id,
      Scal t) const;
  void DumpPartInter(
      const FieldCell<Scal>& fca, const FieldCell<Vect>& fcn, size_t id,
      Scal t) const;
  // Returns curvature field from last call of Part()
  Multi<const FieldCell<Scal>*> GetCurv();

 private:
  struct Imp; // implementation
  // std::unique_ptr<Imp> imp;
  std::unique_ptr<Imp> imp;
};
