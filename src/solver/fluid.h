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
  using Expr = Expression<Scal, IdxCell, 1 + dim * 2>;

 protected:
  FieldCell<Scal>* fcr_;   // density
  FieldCell<Scal>* fcd_;   // dynamic viscosity
  FieldCell<Vect>* fcf_;   // force
  FieldFace<Scal>* ffbp_; //  balanced force projections
  FieldCell<Scal>* fcsv_;  // volume source
  FieldCell<Scal>* fcsm_;  // mass source

 public:
  FluidSolver(double time, double time_step,
              FieldCell<Scal>* fcr,   // density
              FieldCell<Scal>* fcd,   // dynamic viscosity
              FieldCell<Vect>* fcf,   // force 
              FieldFace<Scal>* ffbp, // balanced force projections 
              FieldCell<Scal>* fcsv,  // volume source
              FieldCell<Scal>* fcsm)  // mass source
      : UnsteadyIterativeSolver(time, time_step)
      , fcr_(fcr) , fcd_(fcd) , fcf_(fcf), ffbp_(ffbp)
      , fcsv_(fcsv) , fcsm_(fcsm)
  {}
  virtual const FieldCell<Vect>& GetVelocity() = 0;
  virtual const FieldCell<Vect>& GetVelocity(Layers layer) = 0;
  virtual const FieldCell<Scal>& GetPressure() = 0;
  virtual const FieldCell<Scal>& GetPressure(Layers layer) = 0;
  virtual const FieldFace<Scal>& GetVolumeFlux() = 0;
  virtual const FieldFace<Scal>& GetVolumeFlux(Layers layer) = 0;
  virtual double GetAutoTimeStep() { return GetTimeStep(); }
};

class ConditionFaceFluid : public ConditionFace {
 public:
  ConditionFaceFluid(size_t nci) : ConditionFace(nci) {}
};

class ConditionCellFluid : public ConditionCell {};

namespace fluid_condition {

template <class M>
class NoSlipWall : public ConditionFaceFluid {
  using Vect = typename M::Vect;
 public:
  NoSlipWall(size_t nci) : ConditionFaceFluid(nci) {}
  virtual Vect GetVelocity() const = 0;
};

template <class M>
class NoSlipWallFixed : public NoSlipWall<M> {
  using Vect = typename M::Vect;
  Vect velocity_;
 public:
  NoSlipWallFixed(Vect velocity, size_t nci)
      : NoSlipWall<M>(nci), velocity_(velocity)
  {}
  Vect GetVelocity() const override {
    return velocity_;
  }
};

template <class M>
class Inlet : public ConditionFaceFluid {
  using Vect = typename M::Vect;
 public:
  Inlet(size_t nci) : ConditionFaceFluid(nci) {}
  virtual Vect GetVelocity() const = 0;
  virtual void SetVelocity(Vect velocity) = 0;
};

template <class M>
class InletFixed : public Inlet<M> {
  using Vect = typename M::Vect;
  Vect velocity_;

 public:
  InletFixed(Vect velocity, size_t nci)
      : Inlet<M>(nci)
      , velocity_(velocity)
  {}
  Vect GetVelocity() const override {
    return velocity_;
  }
  void SetVelocity(Vect velocity) override {
    velocity_ = velocity;
  }
};

template <class M>
class InletFlux : public Inlet<M> {
  using Vect = typename M::Vect;
  Vect vel_;
  size_t id_;

 public:
  InletFlux(Vect vel, size_t id, size_t nci)
      : Inlet<M>(nci)
      , vel_(vel)
      , id_(id)
  {}
  Vect GetVelocity() const override {
    return vel_;
  }
  void SetVelocity(Vect vel) override {
    vel_ = vel;
  }
  size_t GetId() {
    return id_;
  }
};

template <class M>
class Outlet : public ConditionFaceFluid {
  using Vect = typename M::Vect;
 public:
  Outlet(size_t nci) : ConditionFaceFluid(nci) {}
  virtual Vect GetVelocity() const = 0;
  virtual void SetVelocity(Vect velocity) = 0;
};

template <class M>
class OutletAuto : public Outlet<M> {
  using Vect = typename M::Vect;
  Vect velocity_;

 public:
  OutletAuto(size_t nci)
      : Outlet<M>(nci)
      , velocity_(Vect::kZero)
  {}
  Vect GetVelocity() const override {
    return velocity_;
  }
  void SetVelocity(Vect velocity) override {
    velocity_ = velocity;
  }
};

template <class M>
class GivenVelocityAndPressure : public ConditionCellFluid {
  using Vect = typename M::Vect;
  using Scal = typename M::Scal;
 public:
  virtual Vect GetVelocity() const = 0;
  virtual Scal GetPressure() const = 0;
};

template <class M>
class GivenVelocityAndPressureFixed : public GivenVelocityAndPressure<M> {
  using Vect = typename M::Vect;
  using Scal = typename M::Scal;
  Vect velocity_;
  Scal pressure_;

 public:
  GivenVelocityAndPressureFixed(Vect velocity, Scal pressure)
      : velocity_(velocity)
      , pressure_(pressure)
  {}
  Vect GetVelocity() const override {
    return velocity_;
  }
  Scal GetPressure() const override {
    return pressure_;
  }
};

template <class M>
class GivenPressure : public ConditionCellFluid {
  using Scal = typename M::Scal;
 public:
  virtual Scal GetPressure() const = 0;
};

template <class M>
class GivenPressureFixed : public GivenPressure<M> {
  using Scal = typename M::Scal;
  Scal pressure_;

 public:
  GivenPressureFixed(Scal pressure)
      : pressure_(pressure)
  {}
  Scal GetPressure() const override {
    return pressure_;
  }
};

} // namespace fluid_condition


template <class M>
std::shared_ptr<ConditionFaceFluid> Parse(std::string argstr,
                                          IdxFace idxface,
                                          size_t nc, // neighbour cell id
                                          const M& m) {
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
  } else {
    throw std::runtime_error("Parse: Unknown boundary condition type");
  }
}

} // namespace solver

