// Created by Petr Karnakov on 26.06.2020
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
struct ParticlesView {
  using Vect = generic::Vect<Scal, 3>;
  std::vector<Vect>& x; // positions
  std::vector<Vect>& v; // velocity
  std::vector<Scal>& r; // radius
  std::vector<Scal>& rho; // density
};

} // namespace generic

template <class M_>
class ParticlesInterface {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using ParticlesView = generic::ParticlesView<Scal>;

  struct Conf {
    Scal mixture_density;
    Scal mixture_viscosity;
    Vect gravity;
  };

  virtual ~ParticlesInterface() {}
  virtual const Conf& GetConf() const = 0;
  virtual void SetConf(Conf) = 0;
  // Makes one time step.
  // dt: time step
  // fe_flux: mixture volume flux
  virtual void Step(Scal dt, const FieldEmbed<Scal>& fe_flux) = 0;
  // Returns view with pointers to fields.
  virtual ParticlesView GetView() const = 0;
  virtual Scal GetTime() const = 0;
  virtual void DumpCsv(std::string path) const = 0;
};

template <class EB_>
class Particles : public ParticlesInterface<typename EB_::M> {
 public:
  using M = typename EB_::M;
  using Base = ParticlesInterface<M>;
  using EB = EB_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using ParticlesView = generic::ParticlesView<Scal>;
  using Conf = typename Base::Conf;
  using UEB = UEmbed<M>;

  // Constructor
  // vfcu: initial volume fraction
  Particles(
      M& m, const EB& eb, const ParticlesView& init, Scal t, Conf conf);
  ~Particles();
  const Conf& GetConf() const;
  void SetConf(Conf);
  void Step(Scal dt, const FieldEmbed<Scal>& fe_flux) override;
  ParticlesView GetView() const override;
  Scal GetTime() const;
  void DumpCsv(std::string path) const override;

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};
