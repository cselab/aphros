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

  struct Conf {};

  virtual ~ElectroInterface() {}
  virtual const Conf& GetConf() const = 0;
  virtual void SetConf(Conf) = 0;
  // Makes one time step.
  // dt: time step
  // fc_permit: electrical permitivity
  // fc_charge: charge density
  virtual void Step(
      Scal dt, const FieldCell<Scal>& fc_permit,
      const FieldCell<Scal>& fc_charge) = 0;
  virtual const FieldCell<Scal>& GetPotential() const = 0;
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
  using UEB = UEmbed<M>;
  template <class T>
  using FieldFaceb = typename EmbedTraits<EB>::template FieldFaceb<T>;

  // Constructor
  // fc_pot: initial volume fraction
  Electro(
      M& m, const EB& eb, const FieldCell<Scal>& fc_pot,
      const MapEmbed<BCond<Scal>>& mebc_pot, Scal time, Conf conf);
  ~Electro();
  const Conf& GetConf() const;
  void SetConf(Conf);
  void Step(
      Scal dt, const FieldCell<Scal>& fc_permit,
      const FieldCell<Scal>& fc_charge) override;
  const FieldCell<Scal>& GetPotential() const override;
  Scal GetTime() const;

 private:
  struct Imp;
  std::unique_ptr<Imp> imp;
};
