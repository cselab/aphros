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
  GRange<size_t> layers;
  Multi<const FieldCell<Scal>*> vfcu; // volume fraction
  Multi<const MapEmbed<BCond<Scal>>*> vmebc; // boundary conditions
};

} // namespace generic

template <class M_>
class TracerInterface {
 public:
  using M = M_;
  static constexpr size_t dim = M::dim;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using TracerView = generic::TracerView<Scal>;

  enum class SlipType { none, stokes, constant };
  struct Slip {
    SlipType type;
    Vect velocity;
  };

  struct Conf {
    size_t layers;
    Multi<Scal> density;
    Multi<Scal> viscosity;
    Multi<Scal> diameter;
    Multi<Slip> slip;
    Multi<Scal> diffusion;
    Vect gravity;
    ConvSc scheme = ConvSc::superbee;
    Multi<const FieldCell<Scal>*> fc_source;
    bool clip = false;
  };

  virtual ~TracerInterface() {}
  virtual const Conf& GetConf() const = 0;
  virtual void SetConf(Conf) = 0;
  // Makes one time step.
  // dt: time step
  // fe_flux: mixture volume flux
  virtual void Step(Scal dt, const FieldEmbed<Scal>& fe_flux) = 0;
  virtual const Multi<FieldCell<Scal>>& GetVolumeFraction() const = 0;
  virtual void SetVolumeFraction(const Multi<FieldCell<Scal>>&) = 0;
  virtual const FieldCell<Scal>& GetMixtureDensity() const = 0;
  virtual const FieldCell<Scal>& GetMixtureViscosity() const = 0;
  virtual const FieldCell<Scal>& GetMixtureVolumeFraction() const = 0;
  // Returns view with pointers to fields.
  virtual TracerView GetView() const = 0;
  virtual Scal GetTime() const = 0;
  virtual Multi<MapEmbed<BCond<Scal>>>& GetBCondMutable() = 0;
};

template <class EB_>
class Tracer : public TracerInterface<typename EB_::M> {
 public:
  using M = typename EB_::M;
  static constexpr size_t dim = M::dim;
  using Base = TracerInterface<M>;
  using EB = EB_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using TracerView = generic::TracerView<Scal>;
  using Conf = typename Base::Conf;
  using UEB = UEmbed<M>;
  template <class T>
  using FieldFaceb = typename EmbedTraits<EB>::template FieldFaceb<T>;
  using SlipType = typename Base::SlipType;

  // Constructor
  // vfcu: initial volume fraction
  Tracer(
      M& m, const EB& eb, const Multi<const FieldCell<Scal>*>& vfcu,
      const Multi<const MapEmbed<BCond<Scal>>*>& vmebc, Scal t, Conf conf);
  ~Tracer();
  const Conf& GetConf() const override;
  void SetConf(Conf) override;
  // dt: time step
  // fe_flux: mixture volume flux
  void Step(Scal dt, const FieldEmbed<Scal>& fe_flux) override;
  const Multi<FieldCell<Scal>>& GetVolumeFraction() const override;
  void SetVolumeFraction(const Multi<FieldCell<Scal>>&) override;
  const FieldCell<Scal>& GetMixtureDensity() const override;
  const FieldCell<Scal>& GetMixtureViscosity() const override;
  const FieldCell<Scal>& GetMixtureVolumeFraction() const override;
  // Returns view with pointers to fields.
  TracerView GetView() const override;
  Scal GetTime() const override;
  Multi<MapEmbed<BCond<Scal>>>& GetBCondMutable() override;

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};
