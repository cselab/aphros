#pragma once

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <memory>
#include <string>

#include "solver.h"

namespace solver {

template <class M_>
class FluidSolver : public UnsteadyIterativeSolver {
  using M = M_;
  static constexpr size_t dim = M::dim;
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

 protected:
  M& m;
  FieldCell<Scal>* fcr_;   // density
  FieldCell<Scal>* fcd_;   // dynamic viscosity
  FieldCell<Vect>* fcf_;   // force
  FieldFace<Scal>* ffbp_; //  balanced force projections
  FieldCell<Scal>* fcsv_;  // volume source
  FieldCell<Scal>* fcsm_;  // mass source

 public:
  // fcr: density
  // fcd: dynamic viscosity
  // fcf: force
  // ffbp: projections of balanced force
  // fcsv: volume source
  // fcsm: mass source
  FluidSolver(double t, double dt, M& m,
              FieldCell<Scal>* fcr, FieldCell<Scal>* fcd, 
              FieldCell<Vect>* fcf, FieldFace<Scal>* ffbp,
              FieldCell<Scal>* fcsv, FieldCell<Scal>* fcsm)
      : UnsteadyIterativeSolver(t, dt)
      , m(m), fcr_(fcr), fcd_(fcd), fcf_(fcf), ffbp_(ffbp)
      , fcsv_(fcsv), fcsm_(fcsm)
  {}
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
  virtual double GetAutoTimeStep() const { return GetTimeStep(); }
};

class CondFaceFluid : public CondFace {
 public:
  CondFaceFluid(size_t nci) : CondFace(nci) {}
};

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
  Vect GetVelocity() const override { return v_; }
  void SetVelocity(const Vect& v) { v_ = v; }

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
  InletFixed(Vect v, size_t nci) : Inlet<M>(nci) , v_(v) {}
  Vect GetVelocity() const override { return v_; }
  void SetVelocity(Vect v) override { v_ = v; }

 private:
  Vect v_;
};

template <class M>
class InletFlux : public Inlet<M> {
 public:
  using Vect = typename M::Vect;
  InletFlux(Vect v, size_t id, size_t nci) : Inlet<M>(nci), v_(v), id_(id) {}
  Vect GetVelocity() const override { return v_; }
  void SetVelocity(Vect v) override { v_ = v; }
  size_t GetId() { return id_; }

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
  OutletAuto(size_t nci) : Outlet<M>(nci) , v_(Vect::kZero) {}
  Vect GetVelocity() const override { return v_; }
  void SetVelocity(Vect v) override { v_ = v; }

 private:
  Vect v_;
};

template <class M>
class SlipWall : public CondFaceFluid {
 public:
  using Vect = typename M::Vect;
  SlipWall(size_t nci) : CondFaceFluid(nci) {}
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
  GivenVelocityAndPressureFixed(Vect v, Scal p) : v_(v) , p_(p) {}
  Vect GetVelocity() const override { return v_; }
  Scal GetPressure() const override { return p_; }
  void SetVelocity(Vect v) { v_ = v; }
  void SetPressure(Scal p) { p_ = p; }

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
  Scal GetPressure() const override { return p_; }
  void SetPressure(Scal p) { p_ = p; }

 private:
  Scal p_;
};

} // namespace fluid_condition


// argstr: argument string
// f: target face 
// nc: target neighbour cell id
template <class M>
std::shared_ptr<CondFaceFluid> Parse(std::string argstr, IdxFace /*f*/,
                                     size_t nc, const M& /*m*/) {
  using namespace fluid_condition;
  using Vect=  typename M::Vect;
  std::stringstream arg(argstr);

  std::string name;
  arg >> name;

  if (name == "wall") {
    // No-slip wall.
    // wall <velocity>
    Vect vel;
    arg >> vel;
    return std::make_shared<NoSlipWallFixed<M>>(vel, nc);
  } else if (name == "inlet") {
    // Fixed velocity inlet.
    // inlet <velocity>
    Vect vel;
    arg >> vel;
    return std::make_shared<InletFixed<M>>(vel, nc);
  } else if (name == "inletflux") {
    // Fixed flux inlet. Flux defined by given velocity is redistributed
    // over all faces with same id.
    // inletflux <velocity> <id>
    Vect vel;
    int id;
    arg >> vel >> id;
    return std::make_shared<InletFlux<M>>(vel, id, nc);
  } else if (name == "outlet") {
    // Outlet. Velocity is extrapolated from neighbour cells and corrected
    // to yield zero total flux over outlet and inlet faces.
    return std::make_shared<OutletAuto<M>>(nc);
  } else if (name == "slipwall") {
    // Zero derivative for both pressure and velocity.
    // TODO: revise, should be inpenetration for velocity
    return std::make_shared<SlipWall<M>>(nc);
  } else {
    throw std::runtime_error("Parse: unknown cond");
  }
}

template <class M>
MapFace<std::shared_ptr<solver::CondFace>> GetVelCond(
    const M& m, const MapFace<std::shared_ptr<solver::CondFaceFluid>>& mff) {
  using Vect = typename M::Vect;
  (void) m;
  MapFace<std::shared_ptr<solver::CondFace>> mf;
  for (auto it : mff) {
    IdxFace i = it.GetIdx();
    solver::CondFaceFluid* cb = it.GetValue().get();
    size_t nci = cb->GetNci();

    using namespace solver;
    using namespace solver::fluid_condition;
    if (auto cd = dynamic_cast<NoSlipWall<M>*>(cb)) {
      mf[i] = std::make_shared<
          CondFaceValFixed<Vect>>(cd->GetVelocity(), nci);
    } else if (auto cd = dynamic_cast<Inlet<M>*>(cb)) {
      mf[i] = std::make_shared<
          CondFaceValFixed<Vect>>(cd->GetVelocity(), nci);
    } else if (auto cd = dynamic_cast<Outlet<M>*>(cb)) {
      mf[i] = std::make_shared<
          CondFaceValFixed<Vect>>(cd->GetVelocity(), nci);
    } else if (auto cd = dynamic_cast<SlipWall<M>*>(cb)) {
      mf[i] = std::make_shared<
          CondFaceReflect>(nci);
    } else {
      throw std::runtime_error("GetVelocityCond: unknown condition");
    }
  }
  return mf;
}


} // namespace solver

