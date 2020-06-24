// Created by Petr Karnakov on 24.06.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <memory>

#include "cond.h"
#include "geom/mesh.h"
#include "multi.h"
#include "solver/embed.h"

namespace generic {

// Piecewise Linear Interface Characterization.
// Describes the multilayer interface.
template <class Scal>
struct TracerView {
  using Vect = generic::Vect<Scal, 3>;
  GRange<size_t> layers;
  Multi<const FieldCell<Scal>*> vfcvf; // field
  const MapEmbed<BCond<Scal>>& mebc; // boundary conditions
};

} // namespace generic

template <class EB_>
class Tracer {
 public:
  using EB = EB_;
  using M = typename EB::M;
  using P = AdvectionSolver<M>;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using MIdx = typename M::MIdx;
  using TracerView = generic::TracerView<Scal>;

  struct Conf {
    size_t layers;
    Multi<Scal> density;
    Multi<Scal> viscosity;
  };

  // Constructor
  // vfcvf: initial volume fraction
  Tracer(
      M& m, const EB& eb, const Multi<const FieldCell<Scal>*>& vfcvf,
      const MapEmbed<BCond<Scal>>& mebc, Scal t, Conf conf);
  ~Vof();
  const EB& GetEmbed() const;
  const Conf& GetConf() const;
  void SetConf(Conf);
  // dt: time step
  // fe_flux: mixture volume flux
  // vfc_src: volume source, or nullptr
  void Step(
      Scal dt, const FieldEmbed<Scal>& fe_flux,
      const Multi<FieldCell<Scal>*>& vfc_src);
  const FieldCell<Scal>& GetVolumeFraction() const;
  const FieldCell<Scal>& GetDensity() const;
  const FieldCell<Scal>& GetViscosity() const;
  // Returns view with pointers to fields.
  TracerView GetView() const;

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};
