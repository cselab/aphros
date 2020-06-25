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
  Multi<const FieldCell<Scal>*> vfcu; // volume fraction
  const MapEmbed<BCond<Scal>>& mebc; // boundary conditions
};

} // namespace generic

template <class M_>
class TracerInterface {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using TracerView = generic::TracerView<Scal>;

  struct Conf {
    size_t layers;
    Multi<Scal> density;
    Multi<Scal> viscosity;
  };

  virtual ~TracerInterface() {}
  virtual const Conf& GetConf() const = 0;
  virtual void SetConf(Conf) = 0;
  // Makes one time step.
  // dt: time step
  // fe_flux: mixture volume flux
  // vfc_src: volume source, or nullptr
  virtual void Step(Scal dt, const FieldEmbed<Scal>& fe_flux) = 0;
  virtual const Multi<FieldCell<Scal>>& GetVolumeFraction() const = 0;
  virtual const FieldCell<Scal>& GetMixtureDensity() const = 0;
  virtual const FieldCell<Scal>& GetMixtureViscosity() const = 0;
  // Returns view with pointers to fields.
  virtual TracerView GetView() const = 0;
};

template <class EB_>
class Tracer : public TracerInterface<typename EB_::M> {
 public:
  using M = typename EB_::M;
  using Base = TracerInterface<M>;
  using EB = EB_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using TracerView = generic::TracerView<Scal>;
  using Conf = typename Base::Conf;
  using UEB = UEmbed<M>;

  // Constructor
  // vfcu: initial volume fraction
  Tracer(
      M& m, const EB& eb, const Multi<const FieldCell<Scal>*>& vfcu,
      const MapEmbed<BCond<Scal>>& mebc, Scal t, Conf conf);
  ~Tracer();
  const Conf& GetConf() const;
  void SetConf(Conf);
  // dt: time step
  // fe_flux: mixture volume flux
  // vfc_src: volume source, or nullptr
  void Step(Scal dt, const FieldEmbed<Scal>& fe_flux) override;
  const Multi<FieldCell<Scal>>& GetVolumeFraction() const override;
  const FieldCell<Scal>& GetMixtureDensity() const override;
  const FieldCell<Scal>& GetMixtureViscosity() const override;
  // Returns view with pointers to fields.
  TracerView GetView() const override;

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};
