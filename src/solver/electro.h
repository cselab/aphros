// Created by Petr Karnakov on 27.09.2020
// Copyright 2020 ETH Zurich

#pragma once

#include <memory>

#include "cond.h"
#include "geom/mesh.h"
#include "solver/embed.h"

template <class M_>
class ElectroInterface {
 public:
  using M = M_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  struct Conf {
    const Vars& var;
    std::shared_ptr<linear::Solver<M>> linsolver;
  };

  struct Stat {
    Scal current = 0;
    Scal potential_min = 0;
    Scal potential_max = 0;
    Scal potential = 0; // Potential difference.
  };

  virtual ~ElectroInterface() {}
  virtual const Conf& GetConf() const = 0;
  // Makes one time step.
  // dt: time step
  // fc_permit: electrical permitivity
  // fc_charge: charge density
  virtual void Step(Scal dt, const FieldCell<Scal>& fc_vf) = 0;
  virtual const FieldCell<Scal>& GetPotential() const = 0;
  virtual const FieldCell<Vect>& GetCurrent() const = 0;
  virtual const FieldEmbed<Scal>& GetFaceCurrent() const = 0;
  virtual const Stat& GetStat() const = 0;
  virtual Scal GetTime() const = 0;
};

template <class EB_>
class Electro : public ElectroInterface<typename EB_::M> {
 public:
  using M = typename EB_::M;
  using Base = ElectroInterface<M>;
  using EB = EB_;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;
  using Conf = typename Base::Conf;
  using Stat = typename Base::Stat;
  using UEB = UEmbed<M>;
  template <class T>
  using FieldFaceb = typename EmbedTraits<EB>::template FieldFaceb<T>;

  // Constructor
  // fc_pot: initial volume fraction
  Electro(
      M& m, const EB& eb, const MapEmbed<BCond<Scal>>& mebc_pot, Scal time,
      const Conf& conf);
  ~Electro();
  Conf& GetConf() const override;
  void Step(Scal dt, const FieldCell<Scal>& fc_vf) override;
  const FieldCell<Scal>& GetPotential() const override;
  const FieldCell<Vect>& GetCurrent() const override;
  const FieldEmbed<Scal>& GetFaceCurrent() const override;
  const Stat& GetStat() const override;
  Scal GetTime() const override;

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};
