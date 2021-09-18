// Created by Petr Karnakov on 25.04.2018
// Copyright 2018 ETH Zurich

#pragma once

#include <cmath>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#include "cond.h"
#include "embed.h"
#include "geom/mesh.h"
#include "solver.h"

template <class M_>
class FluidSolver : public UnsteadyIterativeSolver {
 public:
  using M = M_;
  static constexpr size_t dim = M::dim;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

 protected:
  M& m;
  const FieldCell<Scal>* fcr_; // density
  const FieldCell<Scal>* fcd_; // dynamic viscosity
  const FieldCell<Vect>* fcf_; // force
  const FieldEmbed<Scal>* febp_; //  balanced force projections
  const FieldCell<Scal>* fcsv_; // volume source
  const FieldCell<Scal>* fcsm_; // mass source

 public:
  // fcr: density
  // fcd: dynamic viscosity
  // fcf: force
  // febp: projections of balanced force
  // fcsv: volume source
  // fcsm: mass source
  FluidSolver(
      double t, double dt, M& m_, const FieldCell<Scal>* fcr,
      const FieldCell<Scal>* fcd, const FieldCell<Vect>* fcf,
      const FieldEmbed<Scal>* febp, const FieldCell<Scal>* fcsv,
      const FieldCell<Scal>* fcsm)
      : UnsteadyIterativeSolver(t, dt)
      , m(m_)
      , fcr_(fcr)
      , fcd_(fcd)
      , fcf_(fcf)
      , febp_(febp)
      , fcsv_(fcsv)
      , fcsm_(fcsm) {}
  virtual const FieldCell<Vect>& GetVelocity(Step) const = 0;
  virtual const FieldCell<Vect>& GetVelocity() const {
    return GetVelocity(Step::time_curr);
  }
  virtual const FieldCell<Scal>& GetPressure(Step) const = 0;
  virtual const FieldCell<Scal>& GetPressure() const {
    return GetPressure(Step::time_curr);
  }
  virtual const FieldEmbed<Scal>& GetVolumeFlux(Step) const = 0;
  virtual const FieldEmbed<Scal>& GetVolumeFlux() const {
    return GetVolumeFlux(Step::time_curr);
  }
  virtual double GetAutoTimeStep() const {
    return GetTimeStep();
  }
  virtual const MapEmbed<BCond<Vect>>& GetVelocityCond() const = 0;
};

class CondCellFluid : public CondCell {};

namespace fluid_condition {

template <class M>
class GivenPressure : public CondCellFluid {
 public:
  using Scal = typename M::Scal;
  virtual Scal GetPressure() const = 0;
};

template <class M>
class GivenPressureFixed : public GivenPressure<M> {
 public:
  using Scal = typename M::Scal;
  GivenPressureFixed(Scal p) : p_(p) {}
  Scal GetPressure() const override {
    return p_;
  }
  void SetPressure(Scal p) {
    p_ = p;
  }

 private:
  Scal p_;
};

} // namespace fluid_condition

enum class BCondFluidType {
  wall, // no-slip wall [velocity]
  slipwall, // slip wall [velocity]
  inlet, // inlet [velocity]
  inletflux, // inlet only constraining the total flux [velocity]
  inletpressure, // inlet with given pressure [pressure]
  outlet, // outlet [-]
  outletpressure, // outlet with given pressure [pressure]
  symm, // symmetry [-]
};

// Boundary condition for fluid on face or embedded boundary.
template <class Vect>
struct BCondFluid {
  using Scal = typename Vect::value_type;
  BCondFluidType type = BCondFluidType::wall;
  Vect velocity = Vect(0);
  Scal pressure = 0;
  size_t nci = 0; // neighbor cell id on faces and 0 on embedded boundaries
  size_t inletflux_id = 0;
};

template <class M>
MapEmbed<BCond<typename M::Vect>> ConvertBCondFluidToVelocity(
    const MapEmbed<BCondFluid<typename M::Vect>>& me_fluid) {
  using Vect = typename M::Vect;
  MapEmbed<BCond<Vect>> mebc;

  me_fluid.LoopPairs([&](auto p) {
    const auto cf = p.first;
    auto& bcf = p.second;
    auto& bc = mebc[cf];
    bc.nci = bcf.nci;
    switch (bcf.type) {
      case BCondFluidType::wall:
      case BCondFluidType::inlet:
      case BCondFluidType::inletflux:
      case BCondFluidType::inletpressure:
      case BCondFluidType::outlet:
      case BCondFluidType::outletpressure:
        bc.type = BCondType::dirichlet;
        bc.val = bcf.velocity;
        break;
      case BCondFluidType::slipwall:
      case BCondFluidType::symm:
        bc.type = BCondType::mixed;
        bc.val = bcf.velocity;
        break;
    }
  });
  return mebc;
}
