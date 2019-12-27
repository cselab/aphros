#pragma once

#include <cmath>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#include "cond.h"
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
  FieldCell<Scal>* fcr_; // density
  FieldCell<Scal>* fcd_; // dynamic viscosity
  FieldCell<Vect>* fcf_; // force
  FieldFace<Scal>* ffbp_; //  balanced force projections
  FieldCell<Scal>* fcsv_; // volume source
  FieldCell<Scal>* fcsm_; // mass source

 public:
  // fcr: density
  // fcd: dynamic viscosity
  // fcf: force
  // ffbp: projections of balanced force
  // fcsv: volume source
  // fcsm: mass source
  FluidSolver(
      double t, double dt, M& m, FieldCell<Scal>* fcr, FieldCell<Scal>* fcd,
      FieldCell<Vect>* fcf, FieldFace<Scal>* ffbp, FieldCell<Scal>* fcsv,
      FieldCell<Scal>* fcsm)
      : UnsteadyIterativeSolver(t, dt)
      , m(m)
      , fcr_(fcr)
      , fcd_(fcd)
      , fcf_(fcf)
      , ffbp_(ffbp)
      , fcsv_(fcsv)
      , fcsm_(fcsm) {}
  virtual const FieldCell<Vect>& GetVelocity(Layers) const = 0;
  virtual const FieldCell<Vect>& GetVelocity() const {
    return GetVelocity(Layers::time_curr);
  }
  virtual const FieldCell<Scal>& GetPressure(Layers) const = 0;
  virtual const FieldCell<Scal>& GetPressure() const {
    return GetPressure(Layers::time_curr);
  }
  virtual const FieldFace<Scal>& GetVolumeFlux(Layers) const = 0;
  virtual const FieldFace<Scal>& GetVolumeFlux() const {
    return GetVolumeFlux(Layers::time_curr);
  }
  virtual double GetAutoTimeStep() const {
    return GetTimeStep();
  }
  virtual const MapCondFace& GetVelocityCond() const = 0;
};

class CondFaceFluid : public CondFace {
 public:
  CondFaceFluid(size_t nci) : CondFace(nci) {}
};

using MapCondFaceFluid = MapFace<UniquePtr<CondFaceFluid>>;

class CondCellFluid : public CondCell {};

namespace fluid_condition {

template <class M>
class NoSlipWall : public CondFaceFluid {
 public:
  using Vect = typename M::Vect;
  NoSlipWall(size_t nci) : CondFaceFluid(nci) {}
  virtual Vect GetVelocity() const = 0;
};

template <class M>
class NoSlipWallFixed : public NoSlipWall<M> {
 public:
  using Vect = typename M::Vect;
  NoSlipWallFixed(Vect v, size_t nci) : NoSlipWall<M>(nci), v_(v) {}
  Vect GetVelocity() const override {
    return v_;
  }
  void SetVelocity(const Vect& v) {
    v_ = v;
  }

 private:
  Vect v_;
};

template <class M>
class Inlet : public CondFaceFluid {
 public:
  using Vect = typename M::Vect;
  Inlet(size_t nci) : CondFaceFluid(nci) {}
  virtual Vect GetVelocity() const = 0;
  virtual void SetVelocity(Vect) = 0;
};

template <class M>
class InletFixed : public Inlet<M> {
 public:
  using Vect = typename M::Vect;
  InletFixed(Vect v, size_t nci) : Inlet<M>(nci), v_(v) {}
  Vect GetVelocity() const override {
    return v_;
  }
  void SetVelocity(Vect v) override {
    v_ = v;
  }

 private:
  Vect v_;
};

template <class M>
class InletFlux : public Inlet<M> {
 public:
  using Vect = typename M::Vect;
  InletFlux(Vect v, size_t id, size_t nci) : Inlet<M>(nci), v_(v), id_(id) {}
  Vect GetVelocity() const override {
    return v_;
  }
  void SetVelocity(Vect v) override {
    v_ = v;
  }
  size_t GetId() {
    return id_;
  }

 private:
  Vect v_;
  size_t id_;
};

template <class M>
class Outlet : public CondFaceFluid {
 public:
  using Vect = typename M::Vect;
  Outlet(size_t nci) : CondFaceFluid(nci) {}
  virtual Vect GetVelocity() const = 0;
  virtual void SetVelocity(Vect) = 0;
};

template <class M>
class OutletAuto : public Outlet<M> {
 public:
  using Vect = typename M::Vect;
  OutletAuto(size_t nci) : Outlet<M>(nci), v_(Vect::kZero) {}
  Vect GetVelocity() const override {
    return v_;
  }
  void SetVelocity(Vect v) override {
    v_ = v;
  }

 private:
  Vect v_;
};

template <class M>
class SlipWall : public CondFaceFluid {
 public:
  SlipWall(size_t nci) : CondFaceFluid(nci) {}
};

template <class M>
class Symm : public CondFaceFluid {
 public:
  Symm(size_t nci) : CondFaceFluid(nci) {}
};

template <class M>
class ClearBc : public CondFaceFluid {
 public:
  ClearBc(size_t nci) : CondFaceFluid(nci) {}
};

template <class M>
class GivenVelocityAndPressure : public CondCellFluid {
 public:
  using Vect = typename M::Vect;
  using Scal = typename M::Scal;
  virtual Vect GetVelocity() const = 0;
  virtual Scal GetPressure() const = 0;
};

template <class M>
class GivenVelocityAndPressureFixed : public GivenVelocityAndPressure<M> {
 public:
  using Vect = typename M::Vect;
  using Scal = typename M::Scal;
  GivenVelocityAndPressureFixed(Vect v, Scal p) : v_(v), p_(p) {}
  Vect GetVelocity() const override {
    return v_;
  }
  Scal GetPressure() const override {
    return p_;
  }
  void SetVelocity(Vect v) {
    v_ = v;
  }
  void SetPressure(Scal p) {
    p_ = p;
  }

 private:
  Vect v_;
  Scal p_;
};

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

template <class M>
MapCondFace GetVelCond(const M& m, const MapCondFaceFluid& mff) {
  using Vect = typename M::Vect;
  (void)m;
  MapCondFace mf;
  for (auto it : mff) {
    IdxFace f = it.GetIdx();
    auto& cb = it.GetValue();
    size_t nci = cb->GetNci();

    using namespace fluid_condition;
    if (auto cd = cb.Get<NoSlipWall<M>>()) {
      mf[f].Set<CondFaceValFixed<Vect>>(cd->GetVelocity(), nci);
    } else if (auto cd = cb.Get<Inlet<M>>()) {
      mf[f].Set<CondFaceValFixed<Vect>>(cd->GetVelocity(), nci);
    } else if (auto cd = cb.Get<Outlet<M>>()) {
      mf[f].Set<CondFaceValFixed<Vect>>(cd->GetVelocity(), nci);
    } else if (auto cd = cb.Get<SlipWall<M>>()) {
      mf[f].Set<CondFaceReflect>(nci);
    } else if (auto cd = cb.Get<Symm<M>>()) {
      mf[f].Set<CondFaceReflect>(nci);
    } else {
      throw std::runtime_error("GetVelCond: unknown condition");
    }
  }
  return mf;
}

using MapCondFaceFluid = MapCondFaceFluid;
